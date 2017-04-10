#ifndef NV_GPU_H
#define NV_GPU_H

#include "nv_core_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NV_GPU_THRESHOLD 256
#define NV_GPU_THREAD_MAX 192
#define NV_GPU_PIN_MALLOC 0 // use cudaMallocHost

int nv_gpu_init(void);
int nv_gpu_available(void);
int nv_gpu_block(int n);
int nv_gpu_thread(int n);
int nv_gpu_optz_block();
int nv_gpu_optz_thread();

nv_matrix_t *nv_gpu_matrix_copy(const nv_matrix_t *hostmat);
nv_matrix_t *nv_gpu_matrix_alloc(float **v, int n, int m);
void nv_gpu_matrix_free(nv_matrix_t *devmat);
#ifdef __cplusplus
}
#endif

#endif

