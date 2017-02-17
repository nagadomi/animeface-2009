#include "nv_core.h"
#include "nv_ip_gaussian.h"

void nv_gaussian5x5(nv_matrix_t *dest, int dch, const nv_matrix_t *src, int sch)
{
	int row;
	int kernel_offset = 2;
	static const float kernel[5][5] = {
		{ 1.831564E-002f, 8.208500E-002f, 1.353353E-001f, 8.208500E-002f, 1.831564E-002f },
		{ 8.208500E-002f, 3.678795E-001f, 6.065307E-001f, 3.678795E-001f, 8.208500E-002f },
		{ 1.353353E-001f, 6.065307E-001f, 1.000000E+000f, 6.065307E-001f, 1.353353E-001f },
		{ 8.208500E-002f, 3.678795E-001f, 6.065307E-001f, 3.678795E-001f, 8.208500E-002f },
		{ 1.831564E-002f, 8.208500E-002f, 1.353353E-001f, 8.208500E-002f, 1.831564E-002f }
	};
	static const float scale = 1.621028E-001f;
	nv_matrix_copy(dest, 0, src, 0, dest->m);
/*
	if (scale == 0.0f) {
		float sum = 0.0f;
		int krow, kcol;
		for (krow = 0; krow < 5; ++krow) {
			for (kcol = 0; kcol < 5; ++kcol) {
				float gy = krow - 2.0f;
				float gx = kcol - 2.0f;
				float gaussian = expf(-(gx * gx) / 2.0f) * expf(-(gy * gy) / 2.0f);
				kernel[krow][kcol] = gaussian;
				sum += gaussian;
			}
			printf("\n");
		}
		scale = 1.0f / sum;
	}
*/

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (row = kernel_offset; row < src->rows - kernel_offset; ++row) {
		int col;
		for (col = kernel_offset; col < src->cols - kernel_offset; ++col) {
			int krow, kcol;
			float v = 0.0f;

			for (krow = 0; krow < 5; ++krow) {
				for (kcol = 0; kcol < 5; ++kcol) {
					v += NV_MAT3D_V(src, row + krow - kernel_offset, col + kcol - kernel_offset, sch) 
						 * kernel[krow][kcol];
				}
			}
			NV_MAT3D_V(dest, row, col, dch) = min(v * scale, 255.0f);
		}
	}
}