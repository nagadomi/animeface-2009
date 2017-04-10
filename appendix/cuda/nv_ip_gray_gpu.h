#ifndef NV_IP_GRAY_GPU_H
#define NV_IP_GRAY_GPU_H
#include "nv_core.h"
#include "nv_ip_gray.h"

#ifdef __cplusplus
extern "C" {
#endif

void nv_gray_gpu(nv_matrix_t *gray, const nv_matrix_t *bgr);


#ifdef __cplusplus
}
#endif
#endif
