#include "nv_core.h"
#include "nv_ip_gray.h"


void nv_gray_cpu(nv_matrix_t *gray, const nv_matrix_t *bgr)
{
	int m;

	assert(bgr->n == 3 && gray->n == 1 && gray->m == bgr->m);

#ifdef _OPENMP
#pragma omp parallel for 
#endif
	for (m = 0; m < gray->m; ++m) {
		// bgr->v 0-255
		float g = NV_GRAY_B_RATE * NV_MAT_V(bgr, m, NV_CH_B);
		g += NV_GRAY_G_RATE * NV_MAT_V(bgr, m, NV_CH_G);
		g += NV_GRAY_R_RATE * NV_MAT_V(bgr, m, NV_CH_R);
		NV_MAT_V(gray, m, 0) = g;
	}
}
