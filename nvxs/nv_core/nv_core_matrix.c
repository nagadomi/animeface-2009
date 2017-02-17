#include "nv_core.h"
#include "nv_core_matrix.h"

#if (NV_ENABLE_CUDA && NV_GPU_PIN_MALLOC)
#include <cuda_runtime.h>
#endif


static void *nv_malloc(unsigned long n)
{
	void *mem;
#if (NV_ENABLE_CUDA && NV_GPU_PIN_MALLOC)
	if (nv_gpu_available()) {
		cudaMallocHost(&mem, n);
	} else {
		mem = malloc(n);
	}
#else
	mem = malloc(n);
#endif
	return mem;
}

static void *nv_realloc(void *oldmem, unsigned long n)
{
	void *mem;
	mem = realloc(oldmem, n);
	return mem;
}

nv_matrix_t *nv_matrix_alloc(int n, int m)
{
	void *mem;
	
	int step = n * sizeof(float);//(n + 4 - (n & 3)) * sizeof(float); // SSE2
	int mem_size = step * m + sizeof(nv_matrix_t) + 0x10;
	nv_matrix_t *matrix = (nv_matrix_t *)nv_malloc(mem_size);

	if (matrix == NULL) {
		return NULL;
	}
	mem = ((char *)matrix) + sizeof(nv_matrix_t);
	matrix->v = (float *)(((char *)mem) + 0x10 - ((unsigned int)mem & 0xf));

	matrix->list = 1;

	matrix->n = n;
	matrix->m = m;
	matrix->cols = m;
	matrix->rows = 1;
	matrix->step = step / sizeof(float);
	matrix->alias = 0;
	matrix->list_step = matrix->step * matrix->m;

	return matrix;
}

nv_matrix_t *nv_matrix_realloc(nv_matrix_t *oldmat, int new_m)
{
	unsigned int offset = (unsigned int)(((char *)oldmat->v) - ((char *)oldmat));
	int mem_size = oldmat->step * new_m * sizeof(float) + sizeof(nv_matrix_t) + 0x10;
	nv_matrix_t *matrix = (nv_matrix_t *)nv_realloc(oldmat, mem_size);

	if (matrix == NULL) {
		return NULL;
	}
	matrix->m = new_m;
	matrix->cols = new_m;
	matrix->list_step = matrix->step * matrix->m;
	matrix->v = (float *)(((char *)matrix) + offset);

	return matrix;
}

nv_matrix_t *nv_matrix_list_alloc(int n, int m, int list)
{
	void *mem;

	int step = n * sizeof(float);//(n + 4 - (n & 3)) * sizeof(float); // SSE2
	int mem_size = sizeof(nv_matrix_t) + (step * m + 0x20) * list;
	nv_matrix_t *matrix = (nv_matrix_t *)nv_malloc(mem_size);

	if (matrix == NULL) {
		return NULL;
	}
	mem = ((char *)matrix) + sizeof(nv_matrix_t);
	matrix->v = (float *)(((char *)mem) + 0x10 - ((unsigned int)mem & 0xf));

	matrix->list = list;
	matrix->n = n;
	matrix->m = m;
	matrix->cols = m;
	matrix->rows = 1;
	matrix->step = step / sizeof(float);
	matrix->alias = 0;
	matrix->list_step = matrix->step * matrix->m;

	return matrix;
}

nv_matrix_t *nv_matrix3d_alloc(int n, int rows, int cols)
{
	nv_matrix_t *mat = nv_matrix_alloc(n, rows * cols);
	mat->rows = rows;
	mat->cols = cols;

	return mat;
}

nv_matrix_t *nv_matrix3d_list_alloc(int n, int rows, int cols, int list)
{
	nv_matrix_t *mat = nv_matrix_list_alloc(n, rows * cols, list);
	mat->rows = rows;
	mat->cols = cols;

	return mat;
}

nv_matrix_t *nv_matrix_alias(const nv_matrix_t *parent, int sn, int sm, int n, int m)
{
	nv_matrix_t *matrix = (nv_matrix_t *)malloc(sizeof(nv_matrix_t));
	matrix->list = 1;
	matrix->n = n;
	matrix->m = m;
	matrix->rows = 1;
	matrix->cols = m;
	matrix->v = &parent->v[sm * parent->step + sn];
	matrix->step = parent->step;
	matrix->list_step = matrix->step * matrix->m;
	matrix->alias = 1;

	return matrix;
}

nv_matrix_t *nv_matrix_list_get(const nv_matrix_t *parent, int list)
{
	nv_matrix_t *matrix = (nv_matrix_t *)malloc(sizeof(nv_matrix_t));
	matrix->list = 1;
	matrix->n = parent->n;
	matrix->m = parent->m;
	matrix->rows = parent->rows;
	matrix->cols = parent->cols;
	matrix->v = &NV_MAT_LIST_V(parent, list, 0, 0);
	matrix->step = parent->step;
	matrix->list_step = parent->list_step;
	matrix->alias = 1;

	return matrix;
}

void nv_matrix_zero(nv_matrix_t *mat)
{
	if (mat->list > 1) {
		memset(mat->v, 0, mat->list_step * mat->list * sizeof(float));
	} else {
		memset(mat->v, 0, mat->m * mat->step * sizeof(float));
	}
}

void nv_vector_zero(nv_matrix_t *mat, int m)
{
	memset(&NV_MAT_V(mat, m, 0), 0, mat->step * sizeof(float));
}

void nv_matrix_copy(nv_matrix_t *dest, int dm, const nv_matrix_t *src, int sm, int count_m)
{
	assert(dest->n == src->n);
	memcpy(&NV_MAT_V(dest, dm, 0), &NV_MAT_V(src, sm, 0), dest->step * count_m * sizeof(float));
}

void nv_matrix_free(nv_matrix_t **matrix)
{
	if (*matrix != NULL) {
		free(*matrix);
		*matrix = NULL;
	}
}

void nv_matrix_print(FILE *out, const nv_matrix_t *mat)
{
	int m, n;
	for (n = 0; n < mat->n; ++n) {
		for (m = 0; m < mat->m; ++m) {
			if (m != 0) {
				fprintf(out, " ");
			}
			fprintf(out, "%10E", NV_MAT_V(mat, m, n));
		}
		fprintf(out, "\n");
	}
}

void nv_matrix_dump_c(FILE *out, const nv_matrix_t *mat, const char *name, int static_variable)
{
	int i, m, start_flag = 1;

	fprintf(out, "static float %s_v[%d] = {\n", name, mat->m * mat->step);
	for (m = 0; m < mat->m; ++m) {
		for (i = 0; i < mat->step; ++i) {
			if (!start_flag) {
				fprintf(out, ",");
				if (i != 0 && i % 5 == 0) {
					fprintf(out, "\n");
				}
			} else {
				start_flag = 0;
			}
			fprintf(out, "%15Ef", NV_MAT_V(mat, m, i));
		}
		fprintf(out, "\n");
	}
	fprintf(out, "};\n", name);
	fprintf(out, "%snv_matrix_t %s = {\n %d, %d, %d, %d, %d, %d, %d, %d, %s_v\n};\n",
		static_variable ? "static ":"",
		name, mat->list, mat->list_step, mat->n, mat->m, mat->rows, mat->cols, mat->step, 0,
		name);
	fflush(out);
}

void nv_vector_copy(nv_matrix_t *dest, int dm, const nv_matrix_t *src, int sm)
{
	assert(dest->n == src->n);

	memcpy(&NV_MAT_V(dest, dm, 0), &NV_MAT_V(src, sm, 0), dest->step * sizeof(float));
}

void nv_matrix_m(nv_matrix_t *mat, int m)
{
	if (mat->rows == 1) {
		mat->cols = m;
	} else {
		int diff = mat->m - m;
		if (diff % mat->cols) {
			mat->rows -= (mat->m - m) / mat->cols; 
		} else {
			mat->rows -= (mat->m - m) / mat->cols;
			if (diff > 0) {
				--mat->rows;
			} else {
				++mat->rows;
			}
		}
	}
	mat->m = m;
}