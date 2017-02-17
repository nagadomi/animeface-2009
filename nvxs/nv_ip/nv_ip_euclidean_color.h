#ifndef __NV_IP_EUCLIDEAN_COLOR_H
#define __NV_IP_EUCLIDEAN_COLOR_H

#ifdef __cplusplus
extern "C" {
#endif
#include "nv_core.h"

void nv_color_bgr2euclidean(nv_matrix_t *ec, const nv_matrix_t *bgr);
void nv_color_bgr2euclidean_scalar(nv_matrix_t *ec, int ec_m, const nv_matrix_t *bgr, int bgr_m);

void nv_color_euclidean2bgr_scalar(nv_matrix_t *bgr, int bgr_m, const nv_matrix_t *ec, int ec_m);
void nv_color_euclidean2bgr(nv_matrix_t *bgr, const nv_matrix_t *ec);

#ifdef __cplusplus
}
#endif

#endif
