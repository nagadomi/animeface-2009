#ifndef __NV_ML_GAUSSIAN_H
#define __NV_ML_GAUSSIAN_H

#include "nv_core.h"
#include "nv_num.h"

#ifdef __cplusplus
extern "C" {
#endif

float nv_gaussian_predict(const nv_cov_t *cov, const nv_matrix_t *x, int m);
float nv_gaussian_log_predict(int npca, const nv_cov_t *cov, const nv_matrix_t *x, int xm);
#ifdef __cplusplus
}
#endif

#endif
