#ifndef __NV_IP_GRAY_H
#define __NV_IP_GRAY_H
#include "nv_core.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NV_GRAY_B_RATE 0.11F
#define NV_GRAY_G_RATE 0.59F
#define NV_GRAY_R_RATE 0.30F

void nv_gray(nv_matrix_t *gray, const nv_matrix_t *bgr);
void nv_gray_cpu(nv_matrix_t *gray, const nv_matrix_t *bgr);


#ifdef __cplusplus
}
#endif
#endif
