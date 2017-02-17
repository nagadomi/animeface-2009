#include "nv_core.h"
#include "nv_num_distance.h"
#include "nv_num_vector.h"
#if NV_ENABLE_SSE2
#include <emmintrin.h>
#include <xmmintrin.h>
#endif

float nv_vector_dot(const nv_matrix_t *vec1, int m1,
					const nv_matrix_t *vec2, int m2)
{
	assert(vec1->n == vec2->n);
#if NV_ENABLE_SSE2
	{
		NV_ALIGNED(float, mm[4], 16);
		__m128 w, x, u;
		div_t lp = div(vec1->n, 4);
		int n, pk_lp = lp.quot * 4, nm_lp = lp.rem;
		float dp = 0.0f;

		u = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f);
		for (n = 0; n < pk_lp; n += 4) {
			w = _mm_loadu_ps(&NV_MAT_V(vec1, m1, n));
			x = _mm_loadu_ps(&NV_MAT_V(vec2, m2, n));
			x = _mm_mul_ps(x, w);
			u = _mm_add_ps(u, x);
		}
		_mm_store_ps(mm, u);
		dp = mm[0] + mm[1] + mm[2] + mm[3];
  
		for (n = pk_lp; n < vec1->n; ++n) {
			dp += NV_MAT_V(vec1, m1, n) * NV_MAT_V(vec2, m2, n);
		}
  
		return dp;
	}
#else
	{
		int n;
		float dp = 0.0f;
		for (n = 0; n < vec1->n; ++n) {
			dp += NV_MAT_V(vec1, m1, n) * NV_MAT_V(vec2, m2, n);
		}
		return dp;
	}
#endif
}

float nv_vector_norm(const nv_matrix_t *vec, int vec_m)
{
	int n;
	float norm = 0.0F;

	for (n = 0; n < vec->n; ++n) {
		norm += NV_MAT_V(vec, vec_m, n) * NV_MAT_V(vec, vec_m, n);
	}

	return sqrtf(norm);
}

int nv_vector_nn(const nv_matrix_t *mat,
				 const nv_matrix_t *vec, int vec_m)
{
	int m, min_index = 0;
	float min_dist = FLT_MAX;

	assert(mat->n == vec->n);

	for (m = 0; m < mat->m; ++m) {
		float dist = nv_euclidean2(mat, m, vec, vec_m);
		if (dist < min_dist) {
			min_dist = dist;
			min_index = m;
		}
	}

	return min_index;
}

int nv_vector_max_n(const nv_matrix_t *v, int m)
{
	int n, max_n = -1;
	float v_max = -FLT_MAX;

	for (n = 0; n < v->n; ++n) {
		if (NV_MAT_V(v, m, n) > v_max) {
			max_n = n;
			v_max = NV_MAT_V(v, m, n);
		}
	}
	return max_n;
}

int nv_vector_min_n(const nv_matrix_t *v, int m)
{
	int n, min_n = -1;
	float v_min = FLT_MAX;

	for (n = 0; n < v->n; ++n) {
		if (NV_MAT_V(v, m, n) < v_min) {
			min_n = n;
			v_min = NV_MAT_V(v, m, n);
		}
	}
	return min_n;
}

int nv_vector_maxnorm_m(const nv_matrix_t *v)
{
	int m, max_m = -1;
	float v_max = -FLT_MAX;

	for (m = 0; m < v->m; ++m) {
		float norm = nv_vector_norm(v, m);
		if (norm > v_max) {
			max_m = m;
			v_max = norm;
		}
	}
	return max_m;
}

int nv_vector_minnorm_m(const nv_matrix_t *v)
{
	int m, min_m = -1;
	float v_min = FLT_MAX;

	for (m = 0; m < v->m; ++m) {
		float norm = nv_vector_norm(v, m);
		if (norm < v_min) {
			min_m = m;
			v_min = norm;
		}
	}
	return min_m;
}

float nv_vector_sum(const nv_matrix_t *v, int m)
{
	int n;
	float sum = 0.0f;
	for (n = 0; n < v->n; ++n) {
		sum += NV_MAT_V(v, m, n);
	}
	return sum;
}


int nv_vector_maxsum_m(const nv_matrix_t *v)
{
	int m, max_m = -1;
	float v_max = -FLT_MAX;

	for (m = 0; m < v->m; ++m) {
		float sum = nv_vector_sum(v, m);
		if (sum > v_max) {
			max_m = m;
			v_max = sum;
		}
	}
	return max_m;
}

int nv_vector_minsum_m(const nv_matrix_t *v)
{
	int m, min_m = -1;
	float v_min = FLT_MAX;

	for (m = 0; m < v->m; ++m) {
		float sum = nv_vector_sum(v, m);
		if (sum < v_min) {
			min_m = m;
			v_min = sum;
		}
	}
	return min_m;
}

void nv_vector_mean(nv_matrix_t *mean, int mean_m, const nv_matrix_t *mat)
{
	float factor = 1.0f / mat->m;
	int m;

	assert(mean->n == mat->n);

	nv_vector_zero(mean, mean_m);
	for (m = 0; m < mat->m; ++m) {
		int n;
		for (n = 0; n < mat->n; ++n) {
			NV_MAT_V(mean, mean_m, n) += factor * NV_MAT_V(mat, m, n);
		}
	}
}