#ifndef __NV_FACE_FEATURE_H
#define __NV_FACE_FEATURE_H
#include "nv_core.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NV_FACE_HAARLIKE_DIM 1152


#define NV_INTEGRAL_V(sum, x, y, xw, yh) \
(NV_MAT3D_V((sum), (yh), (xw), 0) \
- NV_MAT3D_V((sum), (yh), (x), 0) \
- (NV_MAT3D_V((sum), (y), (xw), 0) - NV_MAT3D_V((sum), (y), (x), 0))) 


typedef enum {
	NV_NORMALIZE_NONE,
	NV_NORMALIZE_MAX,
	NV_NORMALIZE_NORM
} nv_face_haarlike_normalize_e;

void nv_face_haarlike(nv_face_haarlike_normalize_e type,
					  nv_matrix_t *feature, 
					  int feature_m,
					  const nv_matrix_t *sum,
					  int x, int y, int width, int height);



#ifdef __cplusplus
}
#endif

#endif
