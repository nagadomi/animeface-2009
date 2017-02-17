#ifndef __NV_NUM_VECTOR_H
#define __NV_NUM_VECTOR_H

#include "nv_core.h"
#include "nv_num_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

float nv_vector_dot(const nv_matrix_t *vec1, int m1, const nv_matrix_t *vec2, int m2);
float nv_vector_norm(const nv_matrix_t *v, int m);
int nv_vector_nn(const nv_matrix_t *mat, const nv_matrix_t *vec, int m);
int nv_vector_max_n(const nv_matrix_t *v, int m);
int nv_vector_min_n(const nv_matrix_t *v, int m);
int nv_vector_max_m(const nv_matrix_t *v);
int nv_vector_min_m(const nv_matrix_t *v);
void nv_vector_mean(nv_matrix_t *mean, int mean_m, const nv_matrix_t *mat);

int nv_vector_maxnorm_m(const nv_matrix_t *v);
int nv_vector_maxsum_m(const nv_matrix_t *v);
float nv_vector_sum(const nv_matrix_t *v, int vm);
#ifdef __cplusplus
}
#endif
#endif
