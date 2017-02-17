#ifndef __NV_IP_LAPLACIAN_H
#define __NV_IP_LAPLACIAN_H

#include "nv_core.h"

#ifdef __cplusplus
extern "C"{ 
#endif

void nv_laplacian1(nv_matrix_t *edge, const nv_matrix_t *gray, float level);
void nv_laplacian2(nv_matrix_t *edge, const nv_matrix_t *gray, float level);
void nv_laplacian3(nv_matrix_t *edge, const nv_matrix_t *gray, float level);
void nv_laplacian(nv_matrix_t *edge, const float kernel[3][3], const nv_matrix_t *gray, float level);

#ifdef __cplusplus
}
#endif

#endif
