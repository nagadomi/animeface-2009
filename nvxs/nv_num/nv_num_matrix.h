#ifndef __NV_NUM_MATRIX_H
#define __NV_NUM_MATRIX_H

#include "nv_core.h"

#ifdef __cplusplus
extern "C" {
#endif

nv_matrix_t *nv_matrix_tr(const nv_matrix_t *mat);
nv_matrix_t *nv_matrix3d_tr(const nv_matrix_t *mat);


// blas
typedef enum {
	NV_MAT_TR,
	NV_MAT_NOTR
} nv_matrix_tr_t;

void nv_gemv(nv_matrix_t *y, int ym,
			 nv_matrix_tr_t a_tran,
			 const nv_matrix_t *a,
			 const nv_matrix_t *x,
			 int xm);

int nv_gesv(nv_matrix_t *x, int xm,
			const nv_matrix_t *a, // NxN Matrix
			const nv_matrix_t *b, // N-Vector
			int bm);

int nv_gels(nv_matrix_t *x, int xm,
			const nv_matrix_t *a, // NxN Matrix
			const nv_matrix_t *b, // N-Vector
			int bm);

int nv_gelss(nv_matrix_t *x,
			 nv_matrix_t *s, // NxN Matrix
			 const nv_matrix_t *a, // NxN Matrix
			 const nv_matrix_t *b); // N-Matrix


#ifdef __cplusplus
}
#endif

#endif
