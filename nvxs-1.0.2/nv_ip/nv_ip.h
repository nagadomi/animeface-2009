#ifndef __NV_IP_H
#define __NV_IP_H

#include "nv_core.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "nv_ip_gray.h"
#if NV_ENABLE_CUDA
#include "nv_ip_gray_gpu.h"
#endif
#include "nv_ip_integral.h"
#include "nv_ip_laplacian.h"
#include "nv_ip_gaussian.h"
#include "nv_ip_euclidean_color.h"

#ifdef __cplusplus
}
#endif

#endif
