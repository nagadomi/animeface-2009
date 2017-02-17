#include "nv_core.h"
#include "nv_num.h"
#include "nv_num_distance.h"

// ユークリッド距離
float nv_euclidean(const nv_matrix_t *vec1, int m1, const nv_matrix_t *vec2, int m2)
{
	return sqrtf(nv_euclidean2(vec1, m1, vec2, m2));
}

// ユークリッド距離^2
float nv_euclidean2(const nv_matrix_t *vec1, int m1, const nv_matrix_t *vec2, int m2)
{
	int n;
	float dist = 0.0f;

	assert(vec1->n == vec2->n);

	for (n = 0; n < vec1->n; ++n) {
		dist += (NV_MAT_V(vec1, m1, n) - NV_MAT_V(vec2, m2, n))
			* (NV_MAT_V(vec1, m1, n) - NV_MAT_V(vec2, m2, n));
	}
	return dist;
}

// マハラノビス距離
float nv_mahalanobis(const nv_cov_t *cov, const nv_matrix_t *x, int xm)
{
	int n;
	nv_matrix_t *y = nv_matrix_alloc(x->n, 1);
	nv_matrix_t *x2 = nv_matrix_alloc(x->n, 1);
	float distance;
	float delta2 = 0.0f;

	nv_matrix_zero(y);
	nv_matrix_zero(x2);
	for (n = 0; n < x2->n; ++n) {
		NV_MAT_V(x2, 0, n) = NV_MAT_V(x, xm, n) - NV_MAT_V(cov->u, 0, n);
	}
	nv_gemv(y, 0, NV_MAT_TR, cov->eigen_vec, x2, xm);
	for (n = 0; n < x->n; ++n) {
		float ev = NV_MAT_V(cov->eigen_val, n, 0);
		float xv = NV_MAT_V(y, 0, n);
		delta2 += (xv * xv) / ev;
	}

	distance = sqrtf(delta2);
	nv_matrix_free(&x2);
	nv_matrix_free(&y);

	return distance;
}
