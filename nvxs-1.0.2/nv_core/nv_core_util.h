#ifndef __NV_CORE_UTIL_H
#define __NV_CORE_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#define NV_PI 3.1415926f

typedef struct {
	int x;
	int y;
	int width;
	int height;
} nv_rect_t;

typedef struct {
	float v[4];
} nv_color_t;

void nv_srand_time(void);
void nv_srand(unsigned long seed);
float nv_rand(void);
float nv_gaussian_rand(float average, float variance);
void nv_shaffule_idx(int *a, int start, int end);
void nv_vector_shaffule(nv_matrix_t *mat);
void nv_vector_shaffule_pair(nv_matrix_t *mat1, nv_matrix_t *mat2);


// swap
#define nv_swap(type, a, b) \
{ \
	type nv_swap_temp = a; \
	a = b; \
	b = nv_swap_temp; \
}


#ifdef __cplusplus
}
#endif


#endif
