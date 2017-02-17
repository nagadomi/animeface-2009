#ifndef __NV_ML_KMEANS
#define __NV_ML_KMEANS

#ifdef __cplusplus
extern "C" {
#endif

int nv_kmeans(nv_matrix_t *means,  // k
			  nv_matrix_t *count,  // k
			  nv_matrix_t *labels, // data->m
			  const nv_matrix_t *data,
			  const int k,
			  const int max_epoch);


#ifdef __cplusplus
}
#endif


#endif
