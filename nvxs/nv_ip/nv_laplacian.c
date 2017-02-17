#include "nv_core.h"
#include "nv_ip_laplacian.h"

#define NV_LAPLACIAN_KERNEL_SIZE 3

void nv_laplacian1(nv_matrix_t *edge, const nv_matrix_t *gray, float level)
{
	static const float kernel[3][3] = {
		{ 0.0f, 1.0f, 0.0f },
		{ 1.0f,-4.0f, 1.0f },
		{ 0.0f, 1.0f, 0.0f }
	};
	nv_laplacian(edge, kernel, gray, level);
}

void nv_laplacian2(nv_matrix_t *edge, const nv_matrix_t *gray, float level)
{
	static const float kernel[3][3] = {
		{ 1.0f, 1.0f, 1.0f },
		{ 1.0f,-8.0f, 1.0f },
		{ 1.0f, 1.0f, 1.0f }
	};
	nv_laplacian(edge, kernel, gray, level);
}

void nv_laplacian3(nv_matrix_t *edge, const nv_matrix_t *gray, float level)
{
	static const float kernel[3][3] = {
		{ -1.0f, 2.0f, -1.0f },
		{ 2.0f,-4.0f, 2.0f },
		{ -1.0f, 2.0f, -1.0f }
	};
	nv_laplacian(edge, kernel, gray, level);
}

void nv_laplacian(nv_matrix_t *edge,
				  const float kernel[NV_LAPLACIAN_KERNEL_SIZE][NV_LAPLACIAN_KERNEL_SIZE],
				  const nv_matrix_t *gray,
				  float level)
{
	int row;
	int kernel_offset = NV_LAPLACIAN_KERNEL_SIZE / 2;
	float max_v = -FLT_MAX;
#ifdef _OPENMP
	int i;
	int threads = omp_get_num_threads();
	float *thread_max_v = (float *)malloc(sizeof(float) * threads);
	for (i = 0; i < threads; ++i) {
		thread_max_v[i] = -FLT_MAX;
	}
#endif

	nv_matrix_zero(edge);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
	for (row = kernel_offset; row < gray->rows - kernel_offset; ++row) {
		int col;
#ifdef _OPENMP
		int thread_index = omp_get_thread_num();
#endif
		for (col = kernel_offset; col < gray->cols - kernel_offset; ++col) {
			int krow, kcol;
			float v = 0.0f;

			for (krow = 0; krow < NV_LAPLACIAN_KERNEL_SIZE; ++krow) {
				for (kcol = 0; kcol < NV_LAPLACIAN_KERNEL_SIZE; ++kcol) {
					v += NV_MAT3D_V(gray, row + krow - kernel_offset, col + kcol - kernel_offset, 0) 
						 * kernel[krow][kcol];
				}
			}
			v = max(v, 0.0f);
			if (level != 0.0f) {
				NV_MAT3D_V(edge, row, col, 0) = min(v * level, 255.0f);
			} else {
#ifdef _OPENMP
				thread_max_v[thread_index] = max(v, thread_max_v[thread_index]);
#else
				max_v = max(v, max_v);
#endif
				NV_MAT3D_V(edge, row, col, 0) = v;
			}
		}
	}

#ifdef _OPENMP
	if (level == 0.0f) {
		for (i = 0; i < threads; ++i) {
			max_v = max(max_v, thread_max_v[i]);
		}
	}
#endif

	if (level == 0.0f) {
		float th = max_v / 8.0f;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (row = kernel_offset; row < gray->rows - kernel_offset; ++row) {
			int col;
			for (col = kernel_offset; col < gray->cols - kernel_offset; ++col) {
				if (NV_MAT3D_V(edge, row, col, 0) > th) {
					NV_MAT3D_V(edge, row, col, 0) = 255.0f;
				} else {
					NV_MAT3D_V(edge, row, col, 0) = 0.0f;
				}
			}
		}
	}
#ifdef _OPENMP
	free(thread_max_v);
#endif
}
