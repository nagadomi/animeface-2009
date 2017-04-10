#include <cuda_runtime.h>
#include <cutil.h>
#include "nv_core.h"
#include "nv_ip_gray_gpu.h"

static __global__ 
void nv_gray_kernel(float *gray, const float *bgr, int gray_m);


void nv_gray_gpu(nv_matrix_t *gray, const nv_matrix_t *bgr)
{
	int blocks = nv_gpu_block(gray->m);
	int threads = nv_gpu_thread(gray->m);
	int gray_size = sizeof(float) * gray->m * gray->step;
	int bgr_size = sizeof(float) * bgr->m * bgr->step;
	float *dev_gray;
	float *dev_bgr;

	cudaMalloc((void **)&dev_gray, gray_size);
	cudaMalloc((void **)&dev_bgr, bgr_size);
	cudaMemcpy(dev_gray, gray->v, gray_size, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_bgr, bgr->v, bgr_size, cudaMemcpyHostToDevice);
	nv_gray_kernel<<<blocks, threads, sizeof(float) * 3 * NV_GPU_THREAD_MAX>>>(
		dev_gray, dev_bgr, gray->m
	);
	cudaMemcpy(gray->v, dev_gray, gray_size, cudaMemcpyDeviceToHost);
	CUT_CHECK_ERROR("nv_gray_kernel() execution failed\n");

	cudaFree(dev_gray);
	cudaFree(dev_bgr);
}


__global__ 
void nv_gray_kernel(float *gray, const float *bgr, int gray_m)
{
	extern __shared__ float shared_mem[];
	const int my_m = blockDim.x * blockIdx.x + threadIdx.x;
	const int my_m3 = my_m * 3;
	const int shd_idx = threadIdx.x;
	float g;

	if (my_m >= gray_m) {
		return;
	}
	__syncthreads();
	shared_mem[shd_idx + 0] = bgr[my_m3 + 0];
	shared_mem[shd_idx + NV_GPU_THREAD_MAX] = bgr[my_m3 + 1];
	shared_mem[shd_idx + NV_GPU_THREAD_MAX * 2] = bgr[my_m3 + 2];
	g = NV_GRAY_B_RATE * shared_mem[shd_idx + 0];
	g += NV_GRAY_G_RATE * shared_mem[shd_idx + NV_GPU_THREAD_MAX];
	g += NV_GRAY_R_RATE * shared_mem[shd_idx + NV_GPU_THREAD_MAX * 2];
	gray[my_m] = g;
}
