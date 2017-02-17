#include "nv_core.h"
#include "nv_ip_integral.h"

// Integral Image

void nv_integral(nv_matrix_t *integral, const nv_matrix_t *img, int channel)
{
	int row, col;
	int erow = img->rows + 1;
	int ecol = img->cols + 1;

	assert(
		integral->rows - 1 == img->rows 
		&& integral->cols - 1 == img->cols
	);

	nv_matrix_zero(integral);

	for (row = 1; row < erow; ++row) {
		float col_sum = NV_MAT3D_V(img, row - 1, 0, channel);
		for (col = 1; col < ecol; ++col) {
			float col_val = NV_MAT3D_V(img, row - 1, col - 1, channel);
			NV_MAT3D_V(integral, row, col, 0) = NV_MAT3D_V(integral, row - 1, col, 0) + col_sum + col_val;
			col_sum += col_val;
		}
	}
}