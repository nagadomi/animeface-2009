#ifndef __NV_FACE_DETECTION_H
#define __NV_FACE_DETECTION_H
#include "nv_core.h"
#include "nv_ml.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	float likelihood;
	nv_rect_t face;
	nv_rect_t left_eye;
	nv_rect_t right_eye;
	nv_rect_t nose;
	nv_rect_t mouth;
	nv_rect_t chin;
} nv_face_position_t;

int nv_face_detect(nv_face_position_t *face_pos, 
				   int maxface,
				   const nv_matrix_t *gray_integral, 
				   const nv_matrix_t *edge_integral, 
				   const nv_rect_t *image_size,
				   const nv_mlp_t *dir_mlp,
				   const nv_mlp_t *detector_mlp,
				   const nv_mlp_t **bagging_mlp, int bagging_mlps,
				   const nv_mlp_t *parts_mlp,
				   float step,
				   float scale_factor,
				   float min_window_size
				   );

#if NV_ENABLE_CUDA
#include "nv_face_detect_gpu.h"
#endif

#ifdef __cplusplus
}
#endif

#endif
