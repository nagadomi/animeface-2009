#include "nv_core.h"
#include "nv_ml.h"
#include "nv_ip.h"
#include "nv_num.h"
#include "nv_face_detect.h"
#include "nv_face_feature.h"
#include "nv_face_detect_internal.h"

#define NV_FACE_DETECT_MAX_CANDIDATES 8192

int 
nv_face_detect(nv_face_position_t *face_pos, 
			   int maxface,
			   const nv_matrix_t *gray_integral, 
			   const nv_matrix_t *edge_integral, 
			   const nv_rect_t *image_size,
			   const nv_mlp_t *dir_mlp,
			   const nv_mlp_t *detector_mlp,
			   const nv_mlp_t **bagging_mlp, int bagging_mlps,
			   const nv_mlp_t *parts_mlp,
			   float step,
			   float scale_factor,
			   float min_window_size
			   )
{
	nv_candidate candidates[NV_FACE_DETECT_MAX_CANDIDATES] = {0}; // max
	int i, j;
	float scale = min_window_size / 32.0f;
	float xs = 0.0f;
	float ys = 0.0f;
	int ncandidate = 0;
	float frame_avg;
	int avg_count;
	int nface;
#ifdef _OPENMP
	int threads = omp_get_num_procs();
#else
	int threads = 1;
#endif
	nv_matrix_t **haar = (nv_matrix_t **)malloc(sizeof(nv_matrix_t *) * threads);

	for (i = 0; i < threads; ++i) {
		haar[i] = nv_matrix_alloc(NV_FACE_HAARLIKE_DIM, 1);
	}
	while (min(image_size->width, image_size->height) / scale > min_window_size) {
		int window = (int)(32.0f * scale);
		int ye = (int)((image_size->height - (32.0f * scale) - ys) / (step * scale));
		int xe = (int)((image_size->width - (32.0f * scale) - xs) / (step * scale));

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
		for (i = 0; i < ye * xe; ++i) {
#ifdef _OPENMP
			int thread_idx = omp_get_thread_num();
#else
			int thread_idx = 0;
#endif
			int yi = i / xe;
			int xi = i % xe;
			int y = (int)(ys + (yi * step * scale));
			int x = (int)(xs + (xi * step * scale));
			double z1, z;
			int label;
			int ex = NV_ROUND_INT(x + (scale * 32.0f));
			int ey = NV_ROUND_INT(y + (scale * 32.0f));
			float area;

			if (ex >= image_size->width
				|| ey >= image_size->height)
			{
				continue;
			}
			// エッジで枝刈り
			area = NV_MAT3D_V(edge_integral, ey, ex, 0)
				- NV_MAT3D_V(edge_integral, ey, x, 0)
				- (NV_MAT3D_V(edge_integral, y, ex, 0) - NV_MAT3D_V(edge_integral, y, x, 0));
			if (!nv_is_face_edge(window, scale, area)) {
				continue;
			}

			// 特徴量抽出
			nv_face_haarlike(
				NV_NORMALIZE_MAX,
				haar[thread_idx], 0, 
				gray_integral,
				x, y, window, window);

			// 顔方向判定
			label = nv_mlp_predict_label(dir_mlp, haar[thread_idx], 0);
			if (!(label == 0 )) {
				continue; // 0 = -30°～30° 以外だったらはじく
			}

			// 顔判別1
			z1 = nv_mlp_predict_d(detector_mlp, haar[thread_idx], 0, 0);//
			if (z1 > 0.1) {
				if (bagging_mlps == 0) {
					z = z1;
				} else {
					// 顔判別2
					z = nv_mlp_bagging_predict_d(bagging_mlp, bagging_mlps, haar[thread_idx], 0, 0);
				}
				if (z > 0.5) {
					// 顔
#ifdef _OPENMP
#pragma omp critical
#endif
					if (NV_FACE_DETECT_MAX_CANDIDATES > ncandidate) {
						candidates[ncandidate].x1 = (float)x;
						candidates[ncandidate].y1 = (float)y;
						candidates[ncandidate].x2 = (float)ex;
						candidates[ncandidate].y2 = (float)ey;
						candidates[ncandidate].z = z;
						candidates[ncandidate].flag = 1;
						candidates[ncandidate].parts = NULL;
						++ncandidate;
					}
				}
			}
			
		}
		scale *= scale_factor;
	}

	// 重複領域の除去
	qsort(candidates, ncandidate, sizeof(nv_candidate), nv_candidate_cmp);
	for (i = ncandidate-1; i >= 0; --i) {
		if (!candidates[i].flag) {
			continue;
		}
		for (j = i - 1; j >= 0; --j) {
			float intersect = 0.0f;

			if (!candidates[j].flag) {
				continue;
			}

			intersect = nv_rect_intersect(
				candidates[i].x1, candidates[i].y1, candidates[i].x2, candidates[i].y2, 
				candidates[j].x1, candidates[j].y1, candidates[j].x2, candidates[j].y2);
			if (intersect > 0.1f) {
				if (candidates[i].z == candidates[j].z) {
					if ((candidates[i].x2 - candidates[i].x1) > (candidates[j].x2 - candidates[j].x1)) {
						candidates[i].flag = 0;
						break;
					} else {
						candidates[j].flag = 0;
					}
				} else {
					if (candidates[i].z <= candidates[j].z) {
						candidates[i].flag = 0;
						break;
					} else {
						candidates[j].flag = 0;
					}
				}
			}
		}
	}

	// 部品推定
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
	for (i = 0; i < min(ncandidate, maxface); ++i) {
		if (candidates[i].flag) {
			int window = (int)candidates[i].x2 - (int)candidates[i].x1;
#ifdef _OPENMP
			int thread_idx = omp_get_thread_num();
#else
			int thread_idx = 0;
#endif
			nv_face_haarlike(
				NV_NORMALIZE_NORM,
				haar[thread_idx], 0,
				gray_integral,
				(int)candidates[i].x1, (int)candidates[i].y1, window, window);
			candidates[i].parts = nv_matrix_alloc(parts_mlp->output, 1);
			nv_mlp_regression(parts_mlp, haar[thread_idx], 0, candidates[i].parts, 0);
		}
	}

	// 結果作成
	nface = 0;
	for (i = 0; i < ncandidate && i < maxface; ++i) {
		if (candidates[i].flag) {
			float d = (candidates[i].x2 - candidates[i].x1);

			face_pos[nface].likelihood = (float)candidates[i].z;
			face_pos[nface].face.x = NV_ROUND_INT(candidates[i].x1);
			face_pos[nface].face.y = NV_ROUND_INT(candidates[i].y1);
			face_pos[nface].face.width = NV_ROUND_INT((candidates[i].x2 - candidates[i].x1));
			face_pos[nface].face.height = NV_ROUND_INT((candidates[i].y2 - candidates[i].y1));
			face_pos[nface].right_eye.x = NV_ROUND_INT(candidates[i].x1 + d * (NV_MAT_V(candidates[i].parts, 0, 0) - (NV_MAT_V(candidates[i].parts, 0, 2))));
			face_pos[nface].right_eye.y = NV_ROUND_INT(candidates[i].y1 + d * (NV_MAT_V(candidates[i].parts, 0, 1) - (NV_MAT_V(candidates[i].parts, 0, 3))));
			face_pos[nface].right_eye.width = NV_ROUND_INT(d * NV_MAT_V(candidates[i].parts, 0, 2) * 2.0f);
			face_pos[nface].right_eye.height = NV_ROUND_INT(d * NV_MAT_V(candidates[i].parts, 0, 3) * 2.0f);
			face_pos[nface].left_eye.x = NV_ROUND_INT(candidates[i].x1 + d * (NV_MAT_V(candidates[i].parts, 0, 4) - (NV_MAT_V(candidates[i].parts, 0, 6))));
			face_pos[nface].left_eye.y = NV_ROUND_INT(candidates[i].y1 + d * (NV_MAT_V(candidates[i].parts, 0, 5) - (NV_MAT_V(candidates[i].parts, 0, 7))));
			face_pos[nface].left_eye.width = NV_ROUND_INT(d * NV_MAT_V(candidates[i].parts, 0, 6) * 2.0f);
			face_pos[nface].left_eye.height = NV_ROUND_INT(d * NV_MAT_V(candidates[i].parts, 0, 7) * 2.0f);
			face_pos[nface].nose.x = NV_ROUND_INT(candidates[i].x1 + d * (NV_MAT_V(candidates[i].parts, 0, 8) - 1.0f * (1.0f / 32.0f)));
			face_pos[nface].nose.y = NV_ROUND_INT(candidates[i].y1 + d * (NV_MAT_V(candidates[i].parts, 0, 9) - 1.0f * (1.0f / 32.0f)));
			face_pos[nface].nose.width = NV_ROUND_INT(d * (1.0f / 32.0f * 2.0f));
			face_pos[nface].nose.height = NV_ROUND_INT(d * (1.0f / 32.0f * 2.0f));
			face_pos[nface].mouth.x = NV_ROUND_INT(candidates[i].x1 + d * (NV_MAT_V(candidates[i].parts, 0, 10) - (NV_MAT_V(candidates[i].parts, 0, 12))));
			face_pos[nface].mouth.y = NV_ROUND_INT(candidates[i].y1 + d * (NV_MAT_V(candidates[i].parts, 0, 11) - (NV_MAT_V(candidates[i].parts, 0, 13))));
			face_pos[nface].mouth.width = NV_ROUND_INT(d * NV_MAT_V(candidates[i].parts, 0, 12) * 2.0f);
			face_pos[nface].mouth.height = NV_ROUND_INT(d * NV_MAT_V(candidates[i].parts, 0, 13) * 2.0f);
			face_pos[nface].chin.x = NV_ROUND_INT(candidates[i].x1 + d * (NV_MAT_V(candidates[i].parts, 0, 14) - 1.0f * (1.0f / 32.0f)));
			face_pos[nface].chin.y = NV_ROUND_INT(candidates[i].y1 + d * (NV_MAT_V(candidates[i].parts, 0, 15) - 1.0f * (1.0f / 32.0f)));
			face_pos[nface].chin.width = NV_ROUND_INT(d * (1.0f / 32.0f) * 2.0f);
			face_pos[nface].chin.height = NV_ROUND_INT(d * (1.0f / 32.0f) * 2.0f);
			++nface;
		}
		if (candidates[i].parts != NULL) {
			nv_matrix_free(&candidates[i].parts);
		}
	}

	for (i = 0; i < threads; ++i) {
		nv_matrix_free(&haar[i]);
	}
	free(haar);
	haar = NULL;

	return nface;
}
