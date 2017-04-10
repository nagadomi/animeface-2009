#include <cutil.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include "nv_core.h"
#include "nv_core_gpu.h"

static int nv_gpu_is_available = 0;
static int nv_gpu_sm_count = 0;
static int nv_gpu_thread_max = 0;


int nv_gpu_init(void)
{
#if __DEVICE_EMULATION__
	nv_gpu_is_available = 1;
	nv_gpu_sm_count = 16;
	return 0;
#else
	int count = 0;
	int i = 0;

	cudaGetDeviceCount(&count);
	if(count == 0) {
		return -1;
	}
	
	for(i = 0; i < count; i++) {
		cudaDeviceProp prop;
		if(cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
			if(prop.major >= 1) {
				nv_gpu_sm_count = prop.multiProcessorCount;
				nv_gpu_thread_max = prop.maxThreadsPerBlock;

				// NV_MAX
				if (nv_gpu_thread_max > NV_GPU_THREAD_MAX) {
					nv_gpu_thread_max = NV_GPU_THREAD_MAX;
				}

				break;
			}
		}
	}
	if(i == count) {
		return -1;
	}
	cudaSetDevice(i);
	cublasStatus err = cublasInit();

	nv_gpu_is_available = 1;

	return 0;
#endif
}


int nv_gpu_available(void)
{
	return nv_gpu_is_available;
}


int nv_gpu_block(int n)
{
	if (n < nv_gpu_sm_count) {
		return 1;
	}
	if (n < nv_gpu_sm_count * 32) {
		return n / 32 + (n % 32 != 0 ? 1:0);
	}
	if (n < nv_gpu_sm_count * nv_gpu_thread_max) {
		return nv_gpu_sm_count;
	}
	return n / nv_gpu_thread_max + (n % nv_gpu_thread_max != 0 ? 1:0);
}

int nv_gpu_thread(int n)
{
	if (n < nv_gpu_sm_count) {
		return n;
	}
	if (n < nv_gpu_sm_count * 32) {
		return 32;
	}

	if (n < nv_gpu_sm_count * nv_gpu_thread_max) {
		return n / nv_gpu_sm_count 
			+ (32 - (nv_gpu_sm_count >= 32 ? 0:n % nv_gpu_sm_count));
	}
	return nv_gpu_thread_max;
}

int nv_gpu_optz_block()
{
	return nv_gpu_sm_count;
}

int nv_gpu_optz_thread()
{
	return nv_gpu_thread_max > 32 ? 32: nv_gpu_thread_max;
}

nv_matrix_t *nv_gpu_matrix_copy(const nv_matrix_t *mat)
{
	nv_matrix_t *dev_mat;
	float *dev_v;

	CUDA_SAFE_CALL(cudaMalloc((void **)&dev_mat, sizeof(nv_matrix_t)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&dev_v, mat->list * mat->list_step * sizeof(float)));
	CUDA_SAFE_CALL(cudaMemcpy(dev_mat, mat, sizeof(nv_matrix_t), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dev_v, mat->v, mat->list * mat->list_step * sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(&dev_mat->v, &dev_v, sizeof(float *), cudaMemcpyHostToDevice));

	return dev_mat;
}

nv_matrix_t *nv_gpu_matrix_alloc(float **v, int n, int m)
{
	nv_matrix_t *mat = nv_matrix_alloc(n, m);
	nv_matrix_t *dev_mat;
	float *dev_v;

	CUDA_SAFE_CALL(cudaMalloc((void **)&dev_mat, sizeof(nv_matrix_t)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&dev_v, mat->list * mat->list_step * sizeof(float)));
	CUDA_SAFE_CALL(cudaMemcpy(dev_mat, mat, sizeof(nv_matrix_t), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(&dev_mat->v, &dev_v, sizeof(float *), cudaMemcpyHostToDevice));

	nv_matrix_free(&mat);
	if (v != NULL) {
		*v = dev_v;
	}


	return dev_mat;
}

void nv_gpu_matrix_free(nv_matrix_t *dev_mat)
{
	float *v;
	CUDA_SAFE_CALL(cudaMemcpy(&v, &dev_mat->v, sizeof(float *), cudaMemcpyDeviceToHost));
	cudaFree(v);
	cudaFree(dev_mat);
}
