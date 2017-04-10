#ifndef NV_FACE_DETECTION_GPU_H
#define NV_FACE_DETECTION_GPU_H
#include "nv_core.h"
#include "nv_ml.h"
#include "nv_face_detect.h"

#ifdef __cplusplus
extern "C" {
#endif

int nv_face_detect_gpu(const nv_mlp_t **mlp, int nmlp,
					   const nv_mlp_t *dir_mlp, const nv_mlp_t *parts_mlp,
					   const nv_matrix_t *gray_integral, 
					   const nv_matrix_t *edge_integral, 
					   const nv_rect_t *image_size,
					   nv_face_position_t *face_pos, 
					   int maxface);

#ifdef __cplusplus
}
#endif

#endif