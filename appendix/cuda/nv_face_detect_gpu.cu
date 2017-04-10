#include "nv_face_detect.h"
#include "nv_face_feature.h"
#include "nv_face_detect_internal.h"
#include <cutil.h>
#include <cuda_runtime.h>

#define NV_INTEGRAL_V(sum, x, y, xw, yh) \
(NV_MAT3D_V((sum), (yh), (xw), 0) \
- NV_MAT3D_V((sum), (yh), (x), 0) \
- (NV_MAT3D_V((sum), (y), (xw), 0) - NV_MAT3D_V((sum), (y), (x), 0))) 

#define NV_FACE_FEATURE_MAX 30000

__global__
void
feature_kernel(nv_matrix_t *feature, 
			   nv_matrix_t *gray_integral,
			   const nv_rect_t *rect, int data_step, int data_m);
__global__
void
dir_kernel(nv_matrix_t *feature, int data_m,
		   const nv_matrix_t *input_w, 
		   const nv_matrix_t *hidden_w, 
		   const nv_matrix_t *input_bias, 
		   const nv_matrix_t *hidden_bias, 
		   nv_matrix_t *input_y, 
		   nv_matrix_t *hidden_y,
		   nv_matrix_t *flags
		   );
__global__
void
face_kernel(nv_matrix_t *feature, 
		   const int *idx, int nidx,
		   const nv_matrix_t *input_w, 
		   const nv_matrix_t *hidden_w, 
		   const nv_matrix_t *input_bias, 
		   const nv_matrix_t *hidden_bias, 
		   nv_matrix_t *input_y, 
		   nv_matrix_t *hidden_y,
		   nv_matrix_t *flags
		   );


__host__
static void nv_alloc_device_windows(nv_rect_t **rects, int *nrect,
									const nv_matrix_t *edge_integral,
									float base, float start_scale, float up_scale, 
									int width, int height)
{
	float scale = start_scale;
	float xs = 0.0f;
	float ys = 0.0f;
	int count = 0;
	nv_rect_t *candidates;

	while (NV_MIN(width, height) / scale > base) {
		int yi, ye;
		int window = (int)(32.0f * scale);
		ye = (height - (base * scale) - ys) / (4.0f * scale);
		for (yi = 0; yi < ye; ++yi) {
			int y = (int)(ys + (yi * 4.0f * scale));
			int xi, xe;
			xe = (width - (32.0f * scale) - xs) / (4.0f * scale);
			for (xi = 0; xi < xe; ++xi) {
				++count;
			}
		}
		scale *= up_scale;
	}
	candidates = (nv_rect_t *)malloc(sizeof(nv_rect_t) * count);
	count = 0;
	scale = start_scale;

	while (NV_MIN(width, height) / scale > base) {
		int yi, ye;
		int window = (int)(32.0f * scale);
		ye = (height - (base * scale) - ys) / (4.0f * scale);
		for (yi = 0; yi < ye; ++yi) {
			int y = (int)(ys + (yi * 4.0f * scale));
			int xi, xe;
			xe = (width - (32.0f * scale) - xs) / (4.0f * scale);
			for (xi = 0; xi < xe; ++xi) {
				int x = (int)(xs + (xi * 4.0f * scale));
				int px = x;
				int py = y;
				int ex = (x + ((scale * 32.0f) + 0.5f));
				int ey = (y + ((scale * 32.0f) + 0.5f));
				float area = NV_MAT3D_V(edge_integral, ey, ex, 0)
					- NV_MAT3D_V(edge_integral, ey, px, 0)
					- (NV_MAT3D_V(edge_integral, py, ex, 0) - NV_MAT3D_V(edge_integral, py, px, 0));
				if (!nv_is_face_edge(window, scale, area)) {
					continue;
				}
				candidates[count].x = px;
				candidates[count].y = py;
				candidates[count].width = ex - px;
				candidates[count].height = ey - py;
				++count;
			}
		}
		scale *= up_scale;
	}

	// cuda malloc
	*rects = candidates;
	*nrect = count;
}

int nv_face_detect_gpu(const nv_mlp_t **mlp, int nmlp,
							const nv_mlp_t *dir_mlp, const nv_mlp_t *parts_mlp,
							const nv_matrix_t *gray_integral, 
							const nv_matrix_t *edge_integral, 
							const nv_rect_t *image_size,
							nv_face_position_t *face_pos, 
							int maxface)
{
	nv_rect_t *rects, *dev_rects;
	int nrect, i, flag_m;
	int data_m, data_m_all, data_step;
	nv_matrix_t *dev_feature;
	nv_matrix_t *dev_integral;
	nv_matrix_t *flags;
	nv_matrix_t *label;
	int threads, blocks;
	nv_mlp_t kdir_mlp, kface_mlp;
	nv_matrix_t *kdir_iy, *kdir_hy, *kflags;
	nv_matrix_t *kface_iy, *kface_hy;
	nv_matrix_t *feature = nv_matrix_alloc(NV_FACE_HAARLIKE_DIM, NV_FACE_FEATURE_MAX);
	float *kflags_v;
	int *kidx, idx[NV_FACE_FEATURE_MAX], nidx;
	int face_count;
	int t=nv_clock();
	nv_matrix_t *mat_t;

	CUDA_SAFE_CALL(cudaMalloc((void **)&kidx, sizeof(int) * NV_FACE_FEATURE_MAX));
	dev_feature = nv_gpu_matrix_alloc(NULL, NV_FACE_FEATURE_MAX, NV_FACE_HAARLIKE_DIM);
	dev_integral = nv_gpu_matrix_copy(gray_integral);

	kdir_mlp.input_w = nv_gpu_matrix_copy(dir_mlp->input_w);
	kdir_mlp.input_bias = nv_gpu_matrix_copy(dir_mlp->input_bias);
	kdir_mlp.hidden_w = nv_gpu_matrix_copy(dir_mlp->hidden_w);
	kdir_mlp.hidden_bias = nv_gpu_matrix_copy(dir_mlp->hidden_bias);
	kdir_iy = nv_gpu_matrix_alloc(NULL, NV_FACE_FEATURE_MAX, dir_mlp->input_w->m);
	kdir_hy = nv_gpu_matrix_alloc(NULL, NV_FACE_FEATURE_MAX, dir_mlp->hidden_w->m);
	kface_mlp.input_w = nv_gpu_matrix_copy(mlp[0]->input_w);
	kface_mlp.input_bias = nv_gpu_matrix_copy(mlp[0]->input_bias);
	kface_mlp.hidden_w = nv_gpu_matrix_copy(mlp[0]->hidden_w);
	kface_mlp.hidden_bias = nv_gpu_matrix_copy(mlp[0]->hidden_bias);
	kface_iy = nv_gpu_matrix_alloc(NULL, NV_FACE_FEATURE_MAX, mlp[0]->input_w->m);
	kface_hy = nv_gpu_matrix_alloc(NULL, NV_FACE_FEATURE_MAX, mlp[0]->hidden_w->m);
	kflags = nv_gpu_matrix_alloc(&kflags_v, 1, NV_FACE_FEATURE_MAX);

	nv_alloc_device_windows(&rects, &nrect,
		edge_integral,
		32.0f, 1.44f, 1.095f,
		image_size->width, image_size->height);

	CUDA_SAFE_CALL(cudaMalloc((void **)&dev_rects, sizeof(nv_rect_t) * nrect));
	CUDA_SAFE_CALL(cudaMemcpy(dev_rects, rects, sizeof(nv_rect_t) * nrect, cudaMemcpyHostToDevice));

	flags = nv_matrix_alloc(1, nrect);
	nv_matrix_zero(flags);

	data_step = 0;
	flag_m = 0;
	data_m_all = nrect;
	while (data_m_all > 0) {
		if (data_m_all > NV_FACE_FEATURE_MAX) {
			data_m = NV_FACE_FEATURE_MAX;
		} else {
			data_m = data_m_all;
		}
		// face feature kernel
		blocks = nv_gpu_block(data_m);
		threads = nv_gpu_thread(data_m);
		t = nv_clock();

		// ì¡í•íäèo
		feature_kernel<<<blocks, threads>>>(
			dev_feature, 
			dev_integral, 
			dev_rects, NV_FACE_FEATURE_MAX * data_step, data_m
		);
		CUT_CHECK_ERROR("feature_kernel() execution failed\n");
		cudaThreadSynchronize();
		//printf("feature extraction kernel: %d\n", nv_clock() -t);
		//t = nv_clock();

		// äÁï˚å¸îªíË 0 = flag 1
		dir_kernel<<<blocks, threads>>>(
			dev_feature, data_m,
			kdir_mlp.input_w,
			kdir_mlp.hidden_w,
			kdir_mlp.input_bias,
			kdir_mlp.hidden_bias,
			kdir_iy, kdir_hy,
			kflags
		);
		CUT_CHECK_ERROR("dir_kernel() execution failed\n");
		cudaThreadSynchronize();
		//printf("direction kernel: %d\n", nv_clock() -t);
		//t = nv_clock();
		CUDA_SAFE_CALL(cudaMemcpy(
			&NV_MAT_V(flags, NV_FACE_FEATURE_MAX * data_step, 0), 
			kflags_v, sizeof(float) * data_m, 
			cudaMemcpyDeviceToHost));
		//printf("flag tran: %d\n", nv_clock() -t);
		//t = nv_clock();

		// äÁîªíË
		nidx = 0;
		for (i = 0; i < data_m; ++i) {
			if (NV_MAT_V(flags, NV_FACE_FEATURE_MAX * data_step + i, 0) == 1.0f) {
				idx[nidx] = i;
				++nidx;
			}
		}
		//printf("flag: %d\n", nv_clock() -t);
		//t = nv_clock();
		CUDA_SAFE_CALL(cudaMemcpy(kidx, idx, sizeof(int) * nidx, cudaMemcpyHostToDevice));
		blocks = nv_gpu_block(nidx);
		threads = nv_gpu_thread(nidx);
		face_kernel<<<blocks, threads>>>(
			dev_feature,
			kidx, nidx,
			kface_mlp.input_w,
			kface_mlp.hidden_w,
			kface_mlp.input_bias,
			kface_mlp.hidden_bias,
			kface_iy, kface_hy,
			kflags
			);
		CUT_CHECK_ERROR("face_kernel() execution failed\n");
		cudaThreadSynchronize();
		//printf("face kernel: %d\n", nv_clock() -t);
		//t = nv_clock();
		CUDA_SAFE_CALL(cudaMemcpy(
			&NV_MAT_V(flags, NV_FACE_FEATURE_MAX * data_step, 0), 
			kflags_v, sizeof(float) * data_m, 
			cudaMemcpyDeviceToHost));
		//printf("flag tran: %d\n", nv_clock() -t);
		data_m_all -= data_m;
		++data_step;
	}
	face_count = 0;
	for (i = 0; i < nrect; ++i) {
		if (NV_MAT_V(flags, i, 0) == 1.0f) {
			++face_count;
		}
	}
	//printf("face: %d (%d)\n", face_count, nv_clock()-t);

	return 0;
}

__device__
static float nv_face_feature_filter2_gpu(int type,
									const nv_matrix_t *sum,
									int px, int py,
									float xscale, float yscale)
{
	int ix, iy;
	int i = 0;
	float p = 0.0f, p1 = 0.0f, p2 = 0.0f;
	int area_all = 0;
	int ystep = (int)(1.0f * yscale + 0.5f);
	int xstep = (int)(1.0f * xscale + 0.5f);

	if (type == 1) {
		// |Å_|
		for (i = 0; i < 7; ++i) {
			int ppx = px + (int)((1.0f + i) * xscale + 0.5f);
			int ppy = py + (int)(i * yscale + 0.5f);
			int eex = px + (int)(8.0f * xscale + 0.5f);
			int eey = py + (int)((i + 1) * yscale + 0.5f);

			//printf("p1: %d, %d, %d, %d\n", 1+i,8,i,i+1);

			p1 += NV_INTEGRAL_V(sum, ppx, ppy, eex, eey);
			area_all += (eex - ppx) * (eey - ppy);	
		}
		for (i = 1; i < 8; ++i) {
			int ppx = px + (int)(0.0f * xscale + 0.5f);
			int ppy = py + (int)(i * yscale + 0.5f);
			int eex = px + (int)(i * xscale + 0.5f);
			int eey = py + (int)((i + 1) * yscale + 0.5f);

			//printf("p2: %d, %d, %d, %d\n", 0,i,i,i+1);
			p2 += NV_INTEGRAL_V(sum, ppx, ppy, eex, eey);
			area_all += (eex - ppx) * (eey - ppy);	
		}
		p = (p1 - p2) / (area_all * 255.0f);
	} else {
		// |/|
		for (i = 0; i < 7; ++i) {
			int ppx = px + (int)(0.0f * xscale + 0.5f);
			int ppy = py + (int)(i * yscale + 0.5f);
			int eex = px + (int)((7.0f - i) * xscale + 0.5f);
			int eey = py + (int)((i + 1) * yscale + 0.5f);

			//printf("p1: %d, %d, %d, %d\n", 0, 7-i, i, i+1);
			
			p1 += NV_INTEGRAL_V(sum, ppx, ppy, eex, eey);
			area_all += (eex - ppx) * (eey - ppy);	
		}
		for (i = 1; i < 8; ++i) {
			int ppx = px + (int)((8.0f - i) * xscale + 0.5f);
			int ppy = py + (int)(i * yscale + 0.5f);
			int eex = px + (int)(8.0f * xscale + 0.5f);
			int eey = py + (int)((i + 1) * yscale + 0.5f);
			//printf("p2: %d, %d, %d, %d\n", 8-i, 8, i, i+1);

			p2 += NV_INTEGRAL_V(sum, ppx, ppy, eex, eey);
			area_all += (eex - ppx) * (eey - ppy);	
		}
		p = (p1 - p2) / (area_all * 255.0f);
	}

	return p;
}


__device__
float nv_face_feature_filter_gpu(const float *filter,
								 const nv_matrix_t *sum,
								 int px, int py,
								 float xscale, float yscale)
{
	int ix, iy;
	int i = 0;
	float p = 0.0f;
	int area_all = 0;
	int ystep = (int)(1.0f * yscale + 0.5f);
	int xstep = (int)(1.0f * xscale + 0.5f);
	return 0.0f;

	for (iy = 0; iy < 8; ++iy) {
		int ppy = py + (int)(iy * yscale + 0.5f);
		int eey = ppy + ystep;
		for (ix = 0; ix < 8; ++ix) {
			int ppx = px + (int)(ix * xscale + 0.5f);
			int eex = ppx + xstep;
			float area = NV_MAT3D_V(sum, eey, eex, 0)
				- NV_MAT3D_V(sum, eey, ppx, 0)
				- (NV_MAT3D_V(sum, ppy, eex, 0) - NV_MAT3D_V(sum, ppy, ppx, 0));

			p += area * filter[i];
			if (filter[i] != 0.0f) {
				area_all += (eey - ppy) * (eex - ppx);
			}
			++i;
		}
	}
	p /= area_all * 255.0f;
	return p;
}

__constant__
const float filter_diagonal1_gpu[] = {
	0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
	-1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
	-1.0f, -1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
	-1.0f, -1.0f, -1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f,
	-1.0f, -1.0f, -1.0f, -1.0f, 0.0f, 1.0f, 1.0f, 1.0f,
	-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 0.0f, 1.0f, 1.0f,
	-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 0.0f, 1.0f,
	-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, 0.0f,
};
__constant__
const float filter_diagonal2_gpu[] = {
	1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
	1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, -1.0f,
	1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, -1.0f, -1.0f,
	1.0f, 1.0f, 1.0f, 1.0f, 0.0f, -1.0f, -1.0f, -1.0f,
	1.0f, 1.0f, 1.0f, 0.0f, -1.0f, -1.0f, -1.0f, -1.0f,
	1.0f, 1.0f, 0.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f,
	1.0f, 0.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f,
	0.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f,
};

__device__
void nv_face_feature_gpu(nv_face_haarlike_normalize_e normalize_type,
						 nv_matrix_t *feature, 
						 int feature_m,
						 const nv_matrix_t *sum,
						 int x, int y, int width, int height,
						 __shared__ float *shd_mem, int shd_idx)
{
	int ix, iy, n;
	int py, px, ey, ex, sy, sx, hy, hx;
	float p1, p2, area;
	float scale_area;
	float ey_ex, ey_px, py_ex, py_px;
	float v, vmax, vmin;
	float xscale = width / 32.0f;
	float yscale = height / 32.0f;
	float ystep = yscale;
	float xstep = xscale;
	int hystep = (32 - 8) / 2 * 8;

	// level1
	for (iy = 0; iy < 32-8; iy += 2) {
		py = y + (int)(ystep * iy + 0.5f);
		ey = py + (int)(8.0f * ystep + 0.5f);
		sy = (int)(4.0f * ystep + 0.5f);
		hy = iy / 2;
		for (ix = 0; ix < 32-8; ix += 2) {
			px = x + (int)(xstep * ix + 0.5f);
			ex = px + (int)(8.0f * xstep + 0.5f);
			sx = (int)(4.0f * xstep + 0.5f);
			hx = ix / 2;
			scale_area = 1.0f / ((ex - px) * (ey - py) * 255.0f);
			ey_ex = NV_MAT3D_V(sum, ey, ex, 0);
			ey_px = NV_MAT3D_V(sum, ey, px, 0);
			py_ex = NV_MAT3D_V(sum, py, ex, 0);
			py_px = NV_MAT3D_V(sum, py, px, 0);

			// ëSÉGÉäÉA
			area = ey_ex - ey_px - (py_ex - py_px);

			// 1
			// [ ]
			// [ ]
			p1 = NV_MAT3D_V(sum, py + sy, ex, 0)
				   - NV_MAT3D_V(sum, py + sy, px, 0)
				   - (py_ex - py_px);
			p2 = area - p1;
			if (p1 > p2) {
				shd_mem[NV_GPU_THREAD_MAX * 0 + shd_idx] = (p1 - p2) * scale_area;
				shd_mem[NV_GPU_THREAD_MAX * 1 + shd_idx] = 0.0f;
			} else {
				shd_mem[NV_GPU_THREAD_MAX * 0 + shd_idx] = 0.0f;
				shd_mem[NV_GPU_THREAD_MAX * 1 + shd_idx] = (p2 - p1) * scale_area;
			}

			// 2
			// [ ][ ]
			p1 = NV_MAT3D_V(sum, ey, px + sx, 0)
				- ey_px
				- (NV_MAT3D_V(sum, py, px + sx, 0) - py_px);
			p2 = area - p1;
			if (p1 > p2) {
				shd_mem[NV_GPU_THREAD_MAX * 2 + shd_idx] = (p1 - p2) * scale_area;
				shd_mem[NV_GPU_THREAD_MAX * 3 + shd_idx] = 0.0f;
			} else {
				shd_mem[NV_GPU_THREAD_MAX * 2 + shd_idx] = 0.0f;
				shd_mem[NV_GPU_THREAD_MAX * 3 + shd_idx] = (p2 - p1) * scale_area;
			}

			// 3
			// |/|
			p1 = nv_face_feature_filter2_gpu(1, sum, px, py, xscale, yscale);
			if (p1 > 0.0f) {
				shd_mem[NV_GPU_THREAD_MAX * 4 + shd_idx] = p1;
				shd_mem[NV_GPU_THREAD_MAX * 5 + shd_idx] = 0.0f;
			} else {
				shd_mem[NV_GPU_THREAD_MAX * 4 + shd_idx] = 0.0f;
				shd_mem[NV_GPU_THREAD_MAX * 5 + shd_idx] = -1.0f * p1;
			}

			// 4
			// |Å_|
			p1 = nv_face_feature_filter2_gpu(2, sum, px, py, xscale, yscale);
			if (p1 > 0.0f) {
				shd_mem[NV_GPU_THREAD_MAX * 6 + shd_idx] = p1;
				shd_mem[NV_GPU_THREAD_MAX * 7 + shd_idx] = 0.0f;
			} else {
				shd_mem[NV_GPU_THREAD_MAX * 6 + shd_idx] = 0.0f;
				shd_mem[NV_GPU_THREAD_MAX * 7 + shd_idx] = -1.0f * p1;
			}
			
			for (n = 0; n < 8; ++n) {
				NV_MAT_VT(feature, feature_m, hy * hystep + hx * 8 + n) = shd_mem[NV_GPU_THREAD_MAX * n + shd_idx];
			}
		}
	}
	
	// ê≥ãKâª
	switch (normalize_type) {
	case NV_NORMALIZE_MAX:
		// ç≈ëÂíl=1.0
		vmax = 0.0f;
		vmin = FLT_MAX;
		for (n = 0; n < feature->m; ++n) {
			if (NV_MAT_VT(feature, feature_m, n) > vmax) {
				vmax = NV_MAT_VT(feature, feature_m, n);
			}
			if (NV_MAT_VT(feature, feature_m, n) != 0.0f
				&& NV_MAT_VT(feature, feature_m, n) < vmin) 
			{
				vmin = NV_MAT_VT(feature, feature_m, n);
			}
		}
		if (vmax != 0.0f && vmax > vmin) {
			v = 1.0f / (vmax - vmin);
			for (n = 0; n < feature->m; ++n) {
				if (NV_MAT_VT(feature, feature_m, n) != 0.0f) {
					NV_MAT_VT(feature, feature_m, n) = (NV_MAT_VT(feature, feature_m, n) - vmin) * v;
				}
			}
		}
		break;
	case NV_NORMALIZE_NORM:
		// ÉxÉNÉgÉãÉmÉãÉÄ=1.0
		v = 0.0f;
		for (n = 0; n < feature->m; ++n) {
			v += NV_MAT_VT(feature, feature_m, n) * NV_MAT_VT(feature, feature_m, n);
		}
		if (v != 0.0) {
			v = 1.0f / sqrtf(v);
			for (n = 0; n < feature->m; ++n) {
				NV_MAT_VT(feature, feature_m, n) *= v;
			}
		}
		break;
	case NV_NORMALIZE_NONE:
	default:
		break;
	}
}

__device__
void nv_mlp_predict_gpu(const nv_matrix_t *input_w, 
						 const nv_matrix_t *hidden_w, 
						 const nv_matrix_t *input_bias, 
						 const nv_matrix_t *hidden_bias, 
						 int y_idx,
						 nv_matrix_t *input_y, 
						 nv_matrix_t *hidden_y, 
						 const nv_matrix_t *x, int xm,
						 float *shd_mem, int shd_idx)
{
	int n, m;

	// ì¸óÕëw
	for (m = 0; m < input_w->m; ++m) {
		float y = NV_MAT_V(input_bias, m, 0);
		for (n = 0; n < input_w->n; ++n) {
			y += NV_MAT_VT(x, xm, n) * NV_MAT_V(input_w, m, n);
		}
		NV_MAT_VT(input_y, y_idx, m) = 1.0f / (1.0f + expf(-y));
	}

	// íÜä‘ëw
	for (m = 0; m < hidden_w->m; ++m) {
		float y = NV_MAT_V(hidden_bias, m, 0);
		for (n = 0; n < hidden_w->n; ++n) {
			y += NV_MAT_VT(input_y, y_idx, n) * NV_MAT_V(hidden_w, m, n);
		}
		NV_MAT_VT(hidden_y, y_idx, m) = 1.0f / (1.0f + expf(-y));
	}
}

__global__
void
feature_kernel(nv_matrix_t *feature, 
			   nv_matrix_t *gray_integral,
			   const nv_rect_t *rect, int data_step, int data_m)
{
	__shared__ float shd_mem[NV_GPU_THREAD_MAX * 8];
	int shd_idx = threadIdx.x;
	int my_m = blockDim.x * blockIdx.x + threadIdx.x;
	if (my_m < data_m) {
		__syncthreads();

		nv_face_feature_gpu(
			NV_NORMALIZE_MAX, 
			feature, my_m, 
			gray_integral,
			rect[data_step + my_m].x, rect[data_step + my_m].y,
			rect[data_step + my_m].width, rect[data_step + my_m].height,
			shd_mem, shd_idx);
	}
}

__global__
void
dir_kernel(nv_matrix_t *feature, int data_m,
		   const nv_matrix_t *input_w, 
		   const nv_matrix_t *hidden_w, 
		   const nv_matrix_t *input_bias, 
		   const nv_matrix_t *hidden_bias, 
		   nv_matrix_t *input_y, 
		   nv_matrix_t *hidden_y,
		   nv_matrix_t *flags
		   )
{
	__shared__ float shd_mem[NV_FACE_HAARLIKE_DIM * 2];
	int my_m = blockDim.x * blockIdx.x + threadIdx.x;
	int shd_idx = threadIdx.x;
	float mp;
	int l, n;

	if (my_m < data_m) {
		__syncthreads();
		nv_mlp_predict_gpu(
			input_w, hidden_w, 
			input_bias, hidden_bias,
			my_m,
			input_y, hidden_y,
			feature, my_m,
			shd_mem, shd_idx);
		l = -1; // nega

		for (n = 0; n < hidden_y->m; ++n) {
			if (NV_MAT_VT(hidden_y, my_m, n) > 0.5f
				&&
				mp < NV_MAT_VT(hidden_y, my_m, n)) 
			{
				mp = NV_MAT_VT(hidden_y, my_m, n);
				l = n;
			}
		}

		if (l == 0) {
			NV_MAT_V(flags, my_m, 0) = 1.0f;
		} else {
			NV_MAT_V(flags, my_m, 0) = 0.0f;
		}
	}
}

__global__
void
face_kernel(nv_matrix_t *feature, 
		   const int *idx, int nidx,
		   const nv_matrix_t *input_w, 
		   const nv_matrix_t *hidden_w, 
		   const nv_matrix_t *input_bias, 
		   const nv_matrix_t *hidden_bias, 
		   nv_matrix_t *input_y, 
		   nv_matrix_t *hidden_y,
		   nv_matrix_t *flags
		   )
{
	__shared__ float shd_mem[NV_FACE_HAARLIKE_DIM * 2];
	int shd_idx = threadIdx.x;
	int my_n = blockDim.x * blockIdx.x + threadIdx.x;
	if (my_n < nidx) {
		int my_m = idx[my_n];

		__syncthreads();
		nv_mlp_predict_gpu(
			input_w, hidden_w, 
			input_bias, hidden_bias,
			my_n,
			input_y, hidden_y,
			feature, my_m,
			shd_mem, shd_idx);

		if (NV_MAT_VT(hidden_y, my_n, 0) > 0.01f) {
			NV_MAT_V(flags, my_m, 0) = 1.0f;
		} else {
			NV_MAT_V(flags, my_m, 0) = 0.0f;
		}
	}
}

