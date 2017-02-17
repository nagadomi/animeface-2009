#include "nv_core.h"
#include "nv_num_cov.h"
#include "nv_num_matrix.h"
#include "nv_num_eigen.h"

// ï™éUã§ï™éU

nv_cov_t *nv_cov_alloc(int n)
{
	nv_cov_t *cov = (nv_cov_t *)malloc(sizeof(nv_cov_t));

	cov->n = n;
	cov->u = nv_matrix_alloc(n, 1);
	cov->s = nv_matrix_alloc(n, 1);
	cov->cov = nv_matrix_alloc(n, n);
	cov->eigen_val = nv_matrix_alloc(1, n);
	cov->eigen_vec = nv_matrix_alloc(n, n);
	cov->data_m = 0;

	return cov;
}

void nv_cov_free(nv_cov_t **cov)
{
	nv_matrix_free(&(*cov)->cov);
	nv_matrix_free(&(*cov)->eigen_val);
	nv_matrix_free(&(*cov)->eigen_vec);
	nv_matrix_free(&(*cov)->s);
	nv_matrix_free(&(*cov)->u);

	free(*cov);
	*cov = NULL;
}

void nv_cov_eigen(nv_cov_t *cov, const nv_matrix_t *data)
{
	int info;

	nv_cov(cov->cov, cov->u, cov->s, data);
	info = nv_eigen_dm(cov->eigen_vec, cov->eigen_val, cov->cov);
	cov->data_m = data->m;

	assert(info == 0);
}

void nv_cov(nv_matrix_t *cov,
			nv_matrix_t *u,
			nv_matrix_t *s,
			const nv_matrix_t *data)
{
	int m, n, i;
	int alloc_u = 0;
	int alloc_s = 0;
	float factor = 1.0f / data->m;

	if (u == NULL) {
		u = nv_matrix_alloc(cov->n, 1);
		alloc_u =1;
	}
	if (s == NULL) {
		s = nv_matrix_alloc(cov->n, 1);
		alloc_s =1;
	}
	assert(cov->n == data->n && cov->n == cov->m
		&& u->n == cov->n
		&& s->n == cov->n);

	// ïΩãœ
	nv_matrix_zero(u);
	for (m = 0; m < data->m; ++m) {
		for (n = 0; n < data->n; ++n) {
			NV_MAT_V(u, 0, n) += NV_MAT_V(data, m, n) * factor;
		}
	}

	// è„éOäp ï™éUã§ï™éUçsóÒ
	nv_matrix_zero(cov);
	nv_matrix_zero(s);
	for (n = 0; n < cov->n; ++n) {
		for (m = n; m < cov->m; ++m) {
			float v = 0.0f;
			for (i = 0; i < data->m; ++i) {
				float a = NV_MAT_V(data, i, n) - NV_MAT_V(u, 0, n);
				float b = NV_MAT_V(data, i, m) - NV_MAT_V(u, 0, m);
				v += a * b * factor;
			}
			NV_MAT_V(cov, m, n) = v;
		}
	}

	if (alloc_u) {
		nv_matrix_free(&u);
	}
	if (alloc_s) {
		nv_matrix_free(&s);
	}
}
