#include "nv_core.h"
#include "nv_ip_gray.h"
#if NV_ENABLE_CUDA
#include "nv_ip_gray_gpu.h"
#endif

void nv_gray(nv_matrix_t *gray, const nv_matrix_t *bgr)
{
#if NV_ENABLE_CUDA
	if (nv_gpu_available()) {
		// 遅い　ピクセルの制限あり
		nv_gray_gpu(gray, bgr);
	} else {
		nv_gray_cpu(gray, bgr);
	}
#else
	nv_gray_cpu(gray, bgr);
#endif
}


