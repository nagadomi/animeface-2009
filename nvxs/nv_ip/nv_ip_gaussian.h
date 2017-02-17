#ifndef __NV_IP_GAUSSIAN_H
#define __NV_IP_GAUSSIAN_H
#include "nv_core.h"

#ifdef __cplusplus
extern "C" {
#endif

void nv_gaussian5x5(nv_matrix_t *dest, int dch, const nv_matrix_t *src, int sch);

#ifdef __cplusplus
}
#endif

#endif
