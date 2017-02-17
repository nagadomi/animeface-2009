#include "nv_core.h"
#include "nv_num.h"
#include "nv_ip.h"
#include "nv_face_detect.h"
#include "nv_face_analyze.h"


#define NV_SKIN_COLOR_SAMPLES 1024
#define NV_SKIN_COLOR_CLASS 8

#define NV_HAIR_COLOR_SAMPLES 1024
#define NV_HAIR_COLOR_CLASS 8

#define NV_EYE_COLOR_SAMPLES 256
#define NV_EYE_COLOR_CLASS 5

static int
nv_color_nn(const nv_matrix_t *means, int mean_m,
			const nv_matrix_t *data,
			const nv_matrix_t *labels)
{
	int m, idx = 0;
	float min_dist = FLT_MAX;
	float mean_label = (float)mean_m;

	for (m = 0; m < data->m; ++m) {
		if (NV_MAT_V(labels, m, 0) == mean_label) {
			float dist = nv_euclidean2(means, mean_m, data, m);
			if (dist < min_dist) {
				min_dist = dist;
				idx = m;
			}
		}
	}
	return idx;
}

static void
nv_get_skin_color(nv_color_t *skin_color,
				  nv_cov_t *skin_cov,
				  const nv_face_position_t *face,
				  const nv_matrix_t *img)
{
	nv_matrix_t *sample = nv_matrix_alloc(3, NV_SKIN_COLOR_SAMPLES);
	nv_matrix_t *skin_sample = nv_matrix_alloc(3, NV_SKIN_COLOR_SAMPLES);
	nv_matrix_t *labels = nv_matrix_alloc(1, NV_SKIN_COLOR_SAMPLES);
	nv_matrix_t *means = nv_matrix_alloc(3, NV_SKIN_COLOR_CLASS);
	nv_matrix_t *count = nv_matrix_alloc(1, NV_SKIN_COLOR_CLASS);
	int m, max_label, k, skin_m;
	nv_rect_t skin_rect;

	// 目と口の間
	skin_rect.x = face->right_eye.x;
	skin_rect.y = max(
		face->left_eye.y + face->left_eye.height,
		face->right_eye.y + face->left_eye.height);
	skin_rect.width = face->left_eye.x + face->left_eye.width - face->right_eye.x;
	skin_rect.height = face->mouth.y - skin_rect.y;

	for (m = 0; m < NV_SKIN_COLOR_SAMPLES / 2; ++m) {
		int x = skin_rect.x + (int)(skin_rect.width * nv_rand());
		int y = skin_rect.y + (int)(skin_rect.height * nv_rand());

		nv_color_bgr2euclidean_scalar(sample, m, img, NV_MAT_M(img, y, x));
	}
	// 目と目の間
	skin_rect.x = face->right_eye.x + face->right_eye.width;
	skin_rect.y = min(face->left_eye.y, face->right_eye.y);
	skin_rect.width = face->left_eye.x - skin_rect.x;
	skin_rect.height = face->nose.y - skin_rect.y;

	for (m = NV_SKIN_COLOR_SAMPLES / 2; m < NV_SKIN_COLOR_SAMPLES; ++m) {
		int x = skin_rect.x + (int)(skin_rect.width * nv_rand());
		int y = skin_rect.y + (int)(skin_rect.height * nv_rand());

		nv_color_bgr2euclidean_scalar(sample, m, img, NV_MAT_M(img, y, x));
	}
	// クラスタリング
	k = nv_kmeans(means, count, labels, sample, NV_SKIN_COLOR_CLASS, 0);

	// 肌選択
	max_label = nv_vector_maxsum_m(count);
	skin_color->v[0] = NV_MAT_V(means, max_label, 0);
	skin_color->v[1] = NV_MAT_V(means, max_label, 1);
	skin_color->v[2] = NV_MAT_V(means, max_label, 2);

	// 肌の分散共分散行列,固有値を計算
	skin_m = 0;
	for (m = 0; m < NV_SKIN_COLOR_SAMPLES; ++m) {
		if (NV_MAT_V(labels, m, 0) == (float)max_label) {
			nv_vector_copy(skin_sample, skin_m, sample, m);
			++skin_m;
		}
	}
	nv_matrix_m(skin_sample, skin_m);
	nv_cov_eigen(skin_cov, skin_sample);

	nv_matrix_free(&sample);
	nv_matrix_free(&skin_sample);
	nv_matrix_free(&means);
	nv_matrix_free(&count);
	nv_matrix_free(&labels);
}



static void nv_get_hair_color(nv_color_t *hair_color, 
							  nv_cov_t *hair_cov,
							  const nv_cov_t *skin_cov,
							  const nv_color_t *skin_color,
							  const nv_face_position_t *face,
							  const nv_matrix_t *img)
{
	nv_matrix_t *sample = nv_matrix_alloc(3, NV_HAIR_COLOR_SAMPLES);
	nv_matrix_t *hair_sample = nv_matrix_alloc(3, NV_HAIR_COLOR_SAMPLES);
	nv_matrix_t *labels = nv_matrix_alloc(1, NV_HAIR_COLOR_SAMPLES);
	nv_matrix_t *means = nv_matrix_alloc(3, NV_HAIR_COLOR_CLASS);
	nv_matrix_t *count = nv_matrix_alloc(1, NV_HAIR_COLOR_CLASS);
	nv_matrix_t *skin_likelihood = nv_matrix_alloc(1, NV_HAIR_COLOR_SAMPLES);
	nv_matrix_t *likelihood_means = nv_matrix_alloc(1, NV_HAIR_COLOR_CLASS);
	int m, max_label, skin_label, hair_m;
	nv_rect_t hair_rect;
	int k;

	// 目の上から顔の1/3
	hair_rect.x = face->right_eye.x;
	hair_rect.y = max(0, min(face->left_eye.y, face->right_eye.y) - face->face.height / 3);
	hair_rect.width = face->left_eye.x + face->left_eye.width - hair_rect.x;
	hair_rect.height = min(face->left_eye.y, face->right_eye.y) - hair_rect.y;

	if (hair_rect.y + hair_rect.height >= img->rows) {
		hair_rect.height = img->rows - hair_rect.y;
	}

	for (m = 0; m < sample->m; ++m) {
		int x = hair_rect.x + (int)(hair_rect.width * nv_rand());
		int y = hair_rect.y + (int)(hair_rect.height * nv_rand());

		nv_color_bgr2euclidean_scalar(sample, m, img, NV_MAT_M(img, y, x));
	}

	// 肌除去
	for (m = 0; m < sample->m; ++m) {
		NV_MAT_V(skin_likelihood, m, 0) = nv_gaussian_log_predict(0, skin_cov, sample, m);
	}
	k = nv_kmeans(likelihood_means, count, labels, skin_likelihood, likelihood_means->m, 0);
	skin_label = nv_vector_maxsum_m(likelihood_means);
	hair_m = 0;
	for (m = 0; m < sample->m; ++m) {
		if (NV_MAT_V(labels, m, 0) != (float)skin_label) {
			nv_vector_copy(sample, hair_m, sample, m);
			++hair_m;
		}
	}
	nv_matrix_m(sample, hair_m);

	// クラスタリング
	k = nv_kmeans(means, count, labels, sample, NV_HAIR_COLOR_CLASS, 0);
	// 最大メンバクラス選択

	max_label = nv_vector_maxsum_m(count);
	hair_color->v[0] = NV_MAT_V(means, max_label, 0);
	hair_color->v[1] = NV_MAT_V(means, max_label, 1);
	hair_color->v[2] = NV_MAT_V(means, max_label, 2);

	// 髪の分散共分散行列,固有値を計算
	hair_m = 0;
	for (m = 0; m < sample->m; ++m) {
		if (NV_MAT_V(labels, m, 0) == (float)max_label) {
			nv_vector_copy(hair_sample, hair_m, sample, m);
			++hair_m;
		}
	}
	nv_matrix_m(hair_sample, hair_m);
	nv_cov_eigen(hair_cov, hair_sample);

	nv_matrix_free(&sample);
	nv_matrix_free(&hair_sample);
	nv_matrix_free(&means);
	nv_matrix_free(&count);
	nv_matrix_free(&labels);
	nv_matrix_free(&skin_likelihood);
	nv_matrix_free(&likelihood_means);
}

// [3]countで降順ソート用
static int nv_cmp_means(const void *p1, const void *p2)
{
	float *lhs = (float *)p1;
	float *rhs = (float *)p2;

	if (lhs[3] < rhs[3]) {
		return 1;
	} else if (lhs[3] > rhs[3]) {
		return -1;
	}
	return 0;
}

static void nv_get_eye_color(nv_color_t *eye_colors,
							 nv_matrix_t *sample,
							 const nv_color_t *skin_color,
							 const nv_cov_t *skin_cov,
							 const nv_color_t *hair_color,
							 const nv_cov_t *hair_cov,
							 const nv_rect_t *eye_rect,
							 const nv_matrix_t *img)
{
	int m, c;
	int skin_label, hair_label;
	int color_samples = sample->m;

	nv_cov_t *eye_cov = nv_cov_alloc(skin_cov->n);
	nv_matrix_t *skin_likelihood = nv_matrix_alloc(1, color_samples);
	nv_matrix_t *hair_likelihood = nv_matrix_alloc(1, color_samples);
	nv_matrix_t *labels = nv_matrix_alloc(1, color_samples);
	nv_matrix_t *means = nv_matrix_alloc(3, NV_EYE_COLOR_CLASS);
	nv_matrix_t *likelihood_means = nv_matrix_alloc(1, NV_EYE_COLOR_CLASS);
	nv_matrix_t *count = nv_matrix_alloc(1, NV_EYE_COLOR_CLASS);
	nv_matrix_t *sorted_means = nv_matrix_alloc(3 + 1, NV_EYE_COLOR_CLASS);
	int k;
	int eye_m;
	int x, y;

	// 目の色サンプリング
	m = 0;
	for (y = eye_rect->y; y < eye_rect->y + eye_rect->height; ++y) {
		for (x = eye_rect->x; x < eye_rect->x + eye_rect->width; ++x) {
			nv_color_bgr2euclidean_scalar(sample, m, img, NV_MAT_M(img, y, x));
			++m;
		}
	}

	// 肌と髪を除去

	// 肌除去
	for (m = 0; m < sample->m; ++m) {
		NV_MAT_V(skin_likelihood, m, 0) = nv_gaussian_log_predict(0, skin_cov, sample, m);
	}
	k = nv_kmeans(likelihood_means, count, labels, skin_likelihood, likelihood_means->m, 0);
	skin_label = nv_vector_maxsum_m(likelihood_means);
	eye_m = 0;
	for (m = 0; m < sample->m; ++m) {
		if (NV_MAT_V(labels, m, 0) != (float)skin_label) {
			nv_vector_copy(sample, eye_m, sample, m);
			++eye_m;
		}
	}
	nv_matrix_m(sample, eye_m);

	// 髪除去
	for (m = 0; m < sample->m; ++m) {
		NV_MAT_V(hair_likelihood, m, 0) = nv_gaussian_log_predict(0, hair_cov, sample, m);
	}
	nv_matrix_m(hair_likelihood, sample->m);

	k = nv_kmeans(likelihood_means, count, labels, hair_likelihood, likelihood_means->m, 0);
	hair_label = nv_vector_maxsum_m(likelihood_means);
	eye_m = 0;
	for (m = 0; m < sample->m; ++m) {
		if (NV_MAT_V(labels, m, 0) != (float)hair_label) {
			nv_vector_copy(sample, eye_m, sample, m);
			++eye_m;
		}
	}
	nv_matrix_m(sample, eye_m);

	// 目のクラスタリング
	k = nv_kmeans(means, count, labels, sample, NV_EYE_COLOR_CLASS, 0);

	// カウントでソート
	for (m = 0; m < k; ++m) {
		NV_MAT_V(sorted_means, m, 0) = NV_MAT_V(means, m, 0);
		NV_MAT_V(sorted_means, m, 1) = NV_MAT_V(means, m, 1);
		NV_MAT_V(sorted_means, m, 2) = NV_MAT_V(means, m, 2);
		NV_MAT_V(sorted_means, m, 3) = NV_MAT_V(count, m, 0);
	}
	// ソート (label破壊)
	qsort(sorted_means->v, k, sizeof(float) * sorted_means->step, nv_cmp_means);

	// 色
	c = 0;
	for (m = 0; m < k; ++m) {
		eye_colors[c].v[0] = NV_MAT_V(sorted_means, m, 0);
		eye_colors[c].v[1] = NV_MAT_V(sorted_means, m, 1);
		eye_colors[c].v[2] = NV_MAT_V(sorted_means, m, 2);
		eye_colors[c].v[3] = NV_MAT_V(count, m, 0);
		++c;
		if (c == 4) {
			break;
		}
	}
	nv_matrix_free(&skin_likelihood);
	nv_matrix_free(&hair_likelihood);
	nv_matrix_free(&means);
	nv_matrix_free(&likelihood_means);
	nv_matrix_free(&count);
	nv_matrix_free(&labels);
	nv_matrix_free(&sorted_means);
	nv_cov_free(&eye_cov);
}

static void 
nv_get_eye_colors(nv_color_t *eye_colors, 
				  nv_color_t *left_eye_colors,
				  nv_color_t *right_eye_colors,
				  const nv_color_t *skin_color,
				  const nv_cov_t *skin_cov,
				  const nv_color_t *hair_color,
				  const nv_cov_t *hair_cov,
				  const nv_face_position_t *face,
				  const nv_matrix_t *img)
{
	nv_matrix_t *left_sample = nv_matrix_alloc(3, face->left_eye.height * face->left_eye.width);
	nv_matrix_t *right_sample = nv_matrix_alloc(3, face->right_eye.height * face->right_eye.width);
	int color_samples = left_sample->m + right_sample->m;
	nv_matrix_t *sample = nv_matrix_alloc(3, color_samples);
	nv_matrix_t *labels = nv_matrix_alloc(1, color_samples);
	nv_matrix_t *means = nv_matrix_alloc(3, NV_EYE_COLOR_CLASS);
	nv_matrix_t *count = nv_matrix_alloc(1, NV_EYE_COLOR_CLASS);
	nv_matrix_t *sorted_means = nv_matrix_alloc(3 + 1, NV_EYE_COLOR_CLASS);
	int m, c, epoch;

	// 左目
	nv_get_eye_color(
		left_eye_colors, left_sample,
		skin_color, skin_cov,
		hair_color, hair_cov,
		&face->left_eye, img);
	// 右目
	nv_get_eye_color(
		right_eye_colors, right_sample, 
		skin_color, skin_cov,
		hair_color, hair_cov,
		&face->right_eye, img);

	// サンプル数更新
	color_samples = right_sample->m + left_sample->m;
	nv_matrix_m(sample, color_samples);

	// 両目
	// 目の色サンプリング
	for (m = 0; m < left_sample->m; ++m) {
		NV_MAT_V(sample, m, 0) = NV_MAT_V(left_sample, m, 0);
		NV_MAT_V(sample, m, 1) = NV_MAT_V(left_sample, m, 1);
		NV_MAT_V(sample, m, 2) = NV_MAT_V(left_sample, m, 2);
	}
	for (m = 0; m < right_sample->m; ++m) {
		NV_MAT_V(sample, m + left_sample->m, 0) = NV_MAT_V(right_sample, m, 0);
		NV_MAT_V(sample, m + left_sample->m, 1) = NV_MAT_V(right_sample, m, 1);
		NV_MAT_V(sample, m + left_sample->m, 2) = NV_MAT_V(right_sample, m, 2);
	}
	epoch = nv_kmeans(means, count, labels, sample, NV_EYE_COLOR_CLASS, 0);

	// カウントでソート
	for (m = 0; m < NV_EYE_COLOR_CLASS; ++m) {
		NV_MAT_V(sorted_means, m, 0) = NV_MAT_V(means, m, 0);
		NV_MAT_V(sorted_means, m, 1) = NV_MAT_V(means, m, 1);
		NV_MAT_V(sorted_means, m, 2) = NV_MAT_V(means, m, 2);
		NV_MAT_V(sorted_means, m, 3) = NV_MAT_V(count, m, 0);
	}
	// カウントでソート, !labels破壊!
	qsort(sorted_means->v, NV_EYE_COLOR_CLASS, sizeof(float) * sorted_means->step, nv_cmp_means);

	// 色
	c = 0;
	for (m = 0; m < NV_EYE_COLOR_CLASS; ++m) {
		eye_colors[c].v[0] = NV_MAT_V(sorted_means, m, 0);
		eye_colors[c].v[1] = NV_MAT_V(sorted_means, m, 1);
		eye_colors[c].v[2] = NV_MAT_V(sorted_means, m, 2);
		++c;
		if (c == 4) {
			break;
		}
	}

	nv_matrix_free(&sample);
	nv_matrix_free(&left_sample);
	nv_matrix_free(&right_sample);
	nv_matrix_free(&means);
	nv_matrix_free(&sorted_means);
	nv_matrix_free(&count);
	nv_matrix_free(&labels);
}


float nv_eye_ratio(const nv_face_position_t *face)
{
	float left = (float)face->left_eye.width / face->left_eye.height;
	float right = (float)face->right_eye.width / face->right_eye.height;

	return max(left, right);
}

float nv_face_ratio(const nv_face_position_t *face)
{
	float left = (float)face->left_eye.width / face->left_eye.height;
	float right = (float)face->right_eye.width / face->right_eye.height;
	float eye_height = (float)((left > right) ? face->left_eye.height: face->right_eye.height);
	float face_height = (float)max(face->left_eye.y, face->right_eye.y) - face->chin.y;

	return face_height / eye_height;
}

static void nv_ec2bgr(const nv_color_t *ec, nv_color_t *bgr)
{
	nv_matrix_t *tmp = nv_matrix_alloc(3, 2);
	NV_MAT_V(tmp, 0, 0) = ec->v[0];
	NV_MAT_V(tmp, 0, 1) = ec->v[1];
	NV_MAT_V(tmp, 0, 2) = ec->v[2];
	nv_color_euclidean2bgr_scalar(tmp, 1, tmp, 0);
	bgr->v[0] = NV_MAT_V(tmp, 1, 0);
	bgr->v[1] = NV_MAT_V(tmp, 1, 1);
	bgr->v[2] = NV_MAT_V(tmp, 1, 2);

	nv_matrix_free(&tmp);
}

void 
nv_face_analyze(nv_face_feature_t *feature,
				const nv_face_position_t *face,
				const nv_matrix_t *img)
{
	int i;
	nv_cov_t *skin_cov = nv_cov_alloc(3);
	nv_cov_t *hair_cov = nv_cov_alloc(3);

	nv_get_skin_color(&feature->skin_ec, skin_cov, face, img);
	nv_get_hair_color(&feature->hair_ec, hair_cov, skin_cov, &feature->skin_ec, face, img);
	nv_get_eye_colors(
		feature->eye_ec,
		feature->left_eye_ec,
		feature->right_eye_ec,
		&feature->skin_ec, skin_cov,
		&feature->hair_ec, hair_cov,
		face, img);
	feature->eye_ratio = nv_eye_ratio(face);
	feature->face_ratio = nv_face_ratio(face);

	// 色空間変換
	nv_ec2bgr(&feature->skin_ec, &feature->skin_bgr);
	nv_ec2bgr(&feature->hair_ec, &feature->hair_bgr);
	for (i = 0; i < 4; ++i) {
		nv_ec2bgr(&feature->eye_ec[i], &feature->eye_bgr[i]);
		nv_ec2bgr(&feature->left_eye_ec[i], &feature->left_eye_bgr[i]);
		nv_ec2bgr(&feature->right_eye_ec[i], &feature->right_eye_bgr[i]);
	}

	nv_cov_free(&skin_cov);
	nv_cov_free(&hair_cov);
}

