#include "nv_core.h"
#include "nv_ip_euclidean_color.h"

// ユークリッド空間で比較する用の色空間
// 中身は変わる可能性あり

// 現在
// C1 = (R + G + B) / 3
// C2 = (R + (255 - B)) / 2
// C3 = (R + 2 * (255 - G) + B) / 4
// B  = (4 * C3 - 6 * C2 + 6 * C1 + 255) / 6
// G  = (-4 * C3 + 3 * C1 + 510) / 3
// R  = (4 * C3 + 6 * C2 + 6 * C1 - 1275) / 6

static float v_0_255(float v)
{
	if (v > 255.0f) {
		return 255.0f;
	}
	if (v < 0.0f) {
		return 0.0f;
	}
	return v;
}

void nv_color_euclidean2bgr_scalar(nv_matrix_t *bgr, int bgr_m, const nv_matrix_t *ec, int ec_m)
{
	float c1 = NV_MAT_V(ec, ec_m, 0);
	float c2 = NV_MAT_V(ec, ec_m, 1);
	float c3 = NV_MAT_V(ec, ec_m, 2);

	assert(ec->n == bgr->n && ec->n == 3);

	NV_MAT_V(bgr, bgr_m, NV_CH_B) = floorf((4.0f * c3 - 6.0f * c2 + 6.0f * c1 + 255.0f) / 6.0f);
	NV_MAT_V(bgr, bgr_m, NV_CH_G) = floorf((-4.0f * c3 + 3.0f * c1 + 510.0f) / 3.0f);
	NV_MAT_V(bgr, bgr_m, NV_CH_R) = floorf((4.0f * c3 + 6.0f * c2 + 6.0f * c1 - 1275.0f) / 6.0f);
	NV_MAT_V(bgr, bgr_m, 0) = v_0_255(NV_MAT_V(bgr, bgr_m, 0));
	NV_MAT_V(bgr, bgr_m, 1) = v_0_255(NV_MAT_V(bgr, bgr_m, 1));
	NV_MAT_V(bgr, bgr_m, 2) = v_0_255(NV_MAT_V(bgr, bgr_m, 2));
}


void nv_color_euclidean2bgr(nv_matrix_t *bgr, const nv_matrix_t *ec)
{
	int m;

	assert(ec->n == bgr->n && ec->n == 3);
	assert(ec->m == bgr->m);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (m = 0; m < bgr->m; ++m) {
		nv_color_euclidean2bgr_scalar(bgr, m, ec, m);
	}
}

void nv_color_bgr2euclidean_scalar(nv_matrix_t *ec, int ec_m, const nv_matrix_t *bgr, int bgr_m)
{
	assert(ec->n == bgr->n && ec->n == 3);
	NV_MAT_V(ec, ec_m, 0) = floorf((NV_MAT_V(bgr, bgr_m, NV_CH_R) + NV_MAT_V(bgr, bgr_m, NV_CH_G) + NV_MAT_V(bgr, bgr_m, NV_CH_B)) / 3.0f);
	NV_MAT_V(ec, ec_m, 1) = floorf((NV_MAT_V(bgr, bgr_m, NV_CH_R) + (255.0f - NV_MAT_V(bgr, bgr_m, NV_CH_B))) / 2.0f);
	NV_MAT_V(ec, ec_m, 2) = floorf((NV_MAT_V(bgr, bgr_m, NV_CH_R) + 2.0f * (255.0f - NV_MAT_V(bgr, bgr_m, NV_CH_G)) + NV_MAT_V(bgr, bgr_m, NV_CH_B)) / 4.0f);
}

void nv_color_bgr2euclidean(nv_matrix_t *ec, const nv_matrix_t *bgr)
{
	int m;

	assert(ec->n == bgr->n && ec->n == 3);
	assert(ec->m == bgr->m);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (m = 0; m < ec->m; ++m) {
		nv_color_bgr2euclidean_scalar(ec, m, bgr, m);
	}
}
