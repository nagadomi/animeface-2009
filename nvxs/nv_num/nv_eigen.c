#include "nv_core.h"
#include "nv_num_eigen.h"
#include "nv_num_lapack.h"

#define NV_SYMBOL_N "N"
#define NV_SYMBOL_V "V"
#define NV_LAPACK_SPEC_FAST 1


// 対角行列の固有値,固有ベクトルを求める
int nv_eigen_dm(nv_matrix_t *eigen_vec, 
				nv_matrix_t *eigen_val,
				const nv_matrix_t *dmat)
{
	integer spec = NV_LAPACK_SPEC_FAST;
	integer nil = -1;
	integer nb = 10;
	integer lwork;
	integer vec_step = (integer)eigen_vec->step;
	integer vec_n = (integer)eigen_vec->n;
	real *work;
	integer info = 0;
	//int m, n;

	assert(eigen_vec->n == dmat->n
		&& eigen_vec->m == dmat->m
		&& eigen_val->n == 1
		&& eigen_val->m == eigen_vec->m);

	nv_matrix_zero(eigen_val);
	nv_matrix_copy(eigen_vec, 0, dmat, 0, dmat->m);
	//ilaenv_(&spec, "SSYTRD", "U", &nb, &nil, &nil, &nil);
	lwork = (nb + 2) * eigen_vec->n;
	work = (real *)malloc(sizeof(real) * lwork);

	ssyev_("V", "U", &vec_n, eigen_vec->v, &vec_step,
		eigen_val->v, work, &lwork, 
		&info);

	free(work);

	return info;
}
