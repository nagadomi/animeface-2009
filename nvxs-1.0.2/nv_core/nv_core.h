#ifndef __NNCV_CORE_H
#define __NNCV_CORE_H
#include "nv_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "nv_portable.h"
#include "nv_core_matrix.h"
#include "nv_core_util.h"

#if NV_ENABLE_CUDA
#include "nv_core_gpu.h"
#endif

#ifdef __cplusplus
}
#endif
#endif
