#if !NV_XS

#include "nv_core.h"
#include "nv_num.h"
#include "nv_ml_gaussian.h"

// ƒKƒEƒX•ª•z

float nv_gaussian_predict(const nv_cov_t *cov, const nv_matrix_t *x, int xm)
{
	int n;
	nv_matrix_t *y = nv_matrix_alloc(x->n, 1);
	nv_matrix_t *x2 = nv_matrix_alloc(x->n, 1);
	float p = 1.0f;
	float d = (float)x->n;
	float delta2 = 0.0f;
	float lambda = 1.0f;

	nv_matrix_zero(y);
	nv_matrix_zero(x2);
	for (n = 0; n < x2->n; ++n) {
		NV_MAT_V(x2, 0, n) = NV_MAT_V(x, xm, n) - NV_MAT_V(cov->u, 0, n);
	}
	nv_gemv(y, 0, NV_MAT_TR, cov->eigen_vec, x2, 0);
	for (n = 0; n < x->n; ++n) {
		float ev = NV_MAT_V(cov->eigen_val, n, 0);
		float xv = NV_MAT_V(y, 0, n);
		if (ev > 0.0f) {
			delta2 += (xv * xv) / ev;
			lambda *= sqrtf(ev);
		}
	}
	p = (1.0f / powf(2.0f * NV_PI, d / 2.0f)) * (1.0f / lambda) * expf(-0.5f * delta2);

	nv_matrix_free(&x2);
	nv_matrix_free(&y);

	return p;
}

float nv_gaussian_log_predict(int npca, const nv_cov_t *cov, const nv_matrix_t *x, int xm)
{
	int n;
	nv_matrix_t *y = nv_matrix_alloc(x->n, 1);
	nv_matrix_t *x2 = nv_matrix_alloc(x->n, 1);
	float p = 0.0f;
	float delta2 = 0.0f;
	float lambda = 1.0f;
	float d = (float)x->n;

	assert(npca <= x->n);

	nv_matrix_zero(y);
	nv_matrix_zero(x2);
	for (n = 0; n < x2->n; ++n) {
		NV_MAT_V(x2, 0, n) = NV_MAT_V(x, xm, n) - NV_MAT_V(cov->u, 0, n);
	}
	nv_gemv(y, 0, NV_MAT_TR, cov->eigen_vec, x2, 0);
	for (n = x->n - 1; n >= (x->n - (x->n - npca)); --n) {
	//for (n = 0; n < npca; ++n) {
		float ev = NV_MAT_V(cov->eigen_val, n, 0);
		float xv = NV_MAT_V(y, 0, n);
		if (ev > 0.0f) {
			p += logf(1/sqrtf(2.0f * NV_PI * ev)) - (xv * xv)/(2.0f * ev);
		}
	}
	nv_matrix_free(&x2);
	nv_matrix_free(&y);

	return p;
}
#endif
