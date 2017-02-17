#ifndef __NV_NUM_EIGEN_H
#define __NV_NUM_EIGEN_H
#include "nv_core.h"
#ifdef __cplusplus
extern "C" {
#endif


int nv_eigen_dm(nv_matrix_t *eigen_vec, 
				nv_matrix_t *eigen_val,
				const nv_matrix_t *dmat);


#ifdef __cplusplus
}
#endif
#endif
