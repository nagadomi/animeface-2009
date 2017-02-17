#include "nv_core.h"
#include "nv_num_matrix.h"
#include "nv_num_vector.h"
#include "nv_num_lapack.h"

// 転置
nv_matrix_t *nv_matrix_tr(const nv_matrix_t *mat)
{
	int m, n;
	nv_matrix_t *tr = nv_matrix_alloc(mat->m, mat->n);

	for (m = 0; m < mat->m; ++m) {
		for (n = 0; n < mat->n; ++n) {
			NV_MAT_V(tr, n, m) = NV_MAT_V(mat, m, n);		
		}
	}

	return tr;
}

nv_matrix_t *nv_matrix3d_tr(const nv_matrix_t *mat)
{
	int row, col, n;
	nv_matrix_t *tr = nv_matrix3d_alloc(mat->n, mat->cols, mat->rows);

	for (row = 0; row < mat->rows; ++row) {
		for (col = 0; col < mat->cols; ++col) {
			for (n = 0; n < mat->n; ++n) {
				NV_MAT3D_V(tr, col, row, n) = NV_MAT3D_V(mat, row, col, n);
			}
		}
	}

	return tr;
}

// y = (A, x)
void nv_gemv(nv_matrix_t *y, int ym,
			 nv_matrix_tr_t a_tr,
			 const nv_matrix_t *a,
			 const nv_matrix_t *x,
			 int xm)
{
	float alpha = 1.0f;
	float beta = 0.0f;
	integer o = 1;
	integer a_step = (integer)a->step;
	integer a_m = (integer)a->m;
	integer a_n = (integer)a->n;

	nv_matrix_zero(y);
	sgemv_(a_tr == NV_MAT_TR ? "T":"N",
		&a_m, &a_n, &alpha, a->v, &a_step,
		&NV_MAT_V(x, xm, 0), &o, &beta, &NV_MAT_V(y, ym, 0), &o);
}

//連立一次方程式
int nv_gesv(nv_matrix_t *x, int xm,
			const nv_matrix_t *a, // NxN Matrix
			const nv_matrix_t *b, // N-Vector
			int bm)
{
	nv_matrix_t *t_a = nv_matrix_alloc(a->n, a->m);
	nv_matrix_t *t_b = nv_matrix_alloc(b->n, 1);
	integer *ipiv = (integer *)malloc(sizeof(integer) * b->n);
	integer t_a_n = (integer)t_a->n;
	integer t_a_step = (integer)t_a->step;
	integer t_b_step = (integer)t_b->step;
	integer nrhs = 1;
	integer info = 0;

	assert(a->n == a->m);
	assert(a->n == b->n);

	memset(ipiv, 0, sizeof(integer) * b->n);
	nv_matrix_copy(t_a, 0, a, 0, a->m);
	nv_vector_copy(t_b, 0, b, 0);

	sgesv_(&t_a_n, &nrhs, t_a->v, &t_a_step, ipiv, t_b->v, &t_b_step, &info);

	nv_vector_copy(x, xm, t_b, 0);

	nv_matrix_free(&t_a);
	nv_matrix_free(&t_b);
	free(ipiv);

	return info;
}

//連立一次方程式 最小二乗解
int nv_gels(nv_matrix_t *x, int xm,
			const nv_matrix_t *a, // NxN Matrix
			const nv_matrix_t *b, // N-Vector
			int bm)
{
	nv_matrix_t *t_a = nv_matrix_alloc(a->n, a->m);
	nv_matrix_t *t_b = nv_matrix_alloc(b->n, 1);
	integer *ipiv = (integer *)malloc(sizeof(integer) * b->n);
	integer nil = -1;
	integer nb = 10, spec = 1, lwork;
	real *work;
	integer t_a_m = (integer)t_a->m;
	integer t_a_n = (integer)t_a->n;
	integer t_a_step = (integer)t_a->step;
	integer t_b_step = (integer)t_b->step;
	integer nrhs = 1;
	integer info = 0;

	assert(a->n == a->m);
	assert(a->n == b->n);

	memset(ipiv, 0, sizeof(integer) * b->n);
	nv_matrix_copy(t_a, 0, a, 0, a->m);
	nv_vector_copy(t_b, 0, b, 0);

	//ilaenv_(&spec, "SSYTRD", "U", &nb, &nil, &nil, &nil);
	lwork = (nb + 2) * b->n;
	work = (real *)malloc(sizeof(real) * lwork);
	sgels_("N", &t_a_m, &t_a_n, &nrhs, t_a->v, &t_a_step, t_b->v, &t_b_step, work, &lwork, &info);

	nv_vector_copy(x, xm, t_b, 0);

	nv_matrix_free(&t_a);
	nv_matrix_free(&t_b);
	free(ipiv);
	free(work);

	return info;
}

//連立一次方程式 特異値分解+近似
int nv_gelss(nv_matrix_t *x,
			 nv_matrix_t *s, // NxN Matrix
			 const nv_matrix_t *a, // NxN Matrix
			 const nv_matrix_t *b)
{
	nv_matrix_t *t_a = nv_matrix_alloc(a->n, a->m);
	nv_matrix_t *t_b = nv_matrix_alloc(b->n, b->m);
	nv_matrix_t *t_sv = nv_matrix_alloc(a->n, 1);
	real epsilon = FLT_EPSILON;
	integer rank = a->n;
	integer lwork = 3 * b->n + 2 * b->n;
	real *work = (real *)malloc(sizeof(real) * lwork);
	integer t_a_m = (integer)t_a->m;
	integer t_b_m = (integer)t_b->m;
	integer t_a_n = (integer)t_a->n;
	integer t_a_step = (integer)t_a->step;
	integer t_b_step = (integer)t_b->step;
	integer nrhs = 1;
	integer info = 0;

	assert(a->n == a->m);
	assert(a->n == b->n);
	assert(x->n == b->n);
	assert(x->m == b->m);

	nv_matrix_copy(t_a, 0, a, 0, a->m);
	nv_matrix_copy(t_b, 0, b, 0, b->m);
	nv_matrix_zero(t_sv);

	sgelss_(&t_a_m, &t_a_n, &t_b_m, t_a->v, &t_a->step, t_b->v, &t_b_step, t_sv->v, &epsilon, &rank, work, &lwork, &info);
	nv_matrix_copy(x, 0, t_b, 0, x->m);
	if (s != NULL) {
		nv_matrix_copy(s, 0, t_sv, 0, s->m);
	}
	nv_matrix_free(&t_a);
	nv_matrix_free(&t_b);
	nv_matrix_free(&t_sv);
	free(work);

	return info;
}
