#include "nv_core.h"
#include "nv_core_matrix.h"
#include "nv_core_util.h"
#include <time.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

#define IDUM_N 2048
static long idums[IDUM_N] = { 8273 };

void nv_srand_time(void)
{
	unsigned long t = (unsigned long)time(NULL);
	int i;
	for (i = 0; i < IDUM_N; ++i) {
		idums[i] = t + i * 3;
	}
}

void nv_srand(unsigned long seed)
{
	int i;
	for (i = 0; i < IDUM_N; ++i) {
		idums[i] = seed + i * 3;
	}
}

float nv_rand(void)
{
	long k;
	float ans;
#ifdef _OPENMP
	int thread_idx = omp_get_thread_num();
#else
	int thread_idx = 0;
#endif
	long idum = idums[thread_idx];
	idum ^= MASK;
	k = idum / IQ;
	idum = IA * (idum - k * IQ) - IR * k;
	if (idum < 0) { 
		idum += IM; 
	}
	ans = (float)(AM * idum);
	idum ^= MASK;
	idums[thread_idx] = idum;

	return ans;
}

float nv_gaussian_rand(float average, float variance)
{
	float a = nv_rand();
	float b = nv_rand();
	return sqrtf(-2.0f * logf(a)) * sinf(2.0f * NV_PI * b) * variance + average;
}

void nv_shaffule_idx(int *a, int start, int end)
{
	int i, tmp;
	int n = end - start;
	int rnd;

	for (i = 0; i < n; ++i) {
		a[i] = i + start;
	}
	for (i = 1; i < n; ++i) {
		rnd = (int)(nv_rand() * i);
		if (rnd < n) {
			tmp = a[i];
			a[i] = a[rnd];
			a[rnd] = tmp;
		}
	}
}

void nv_vector_shaffule(nv_matrix_t *mat)
{
	int i, rnd;
	int n = mat->m;
	nv_matrix_t *tmp = nv_matrix_alloc(mat->n, 1);

	for (i = 1; i < n; ++i) {
		rnd = (int)(nv_rand() * i);
		if (rnd < n) {
			nv_vector_copy(tmp, 0, mat, i);
			nv_vector_copy(mat, i, mat, rnd);
			nv_vector_copy(mat, rnd, tmp, 0);
		}
	}

	nv_matrix_free(&tmp);
}

void nv_vector_shaffule_pair(nv_matrix_t *mat1, nv_matrix_t *mat2)
{
	int i, rnd;
	int n = mat1->m;
	nv_matrix_t *tmp1 = nv_matrix_alloc(mat1->n, 1);
	nv_matrix_t *tmp2 = nv_matrix_alloc(mat2->n, 1);

	for (i = 1; i < n; ++i) {
		rnd = (int)(nv_rand() * i);
		if (rnd < n) {
			// mat1
			nv_vector_copy(tmp1, 0, mat1, i);
			nv_vector_copy(mat1, i, mat1, rnd);
			nv_vector_copy(mat1, rnd, tmp1, 0);
			// mat2
			nv_vector_copy(tmp2, 0, mat2, i);
			nv_vector_copy(mat2, i, mat2, rnd);
			nv_vector_copy(mat2, rnd, tmp2, 0);
		}
	}
	nv_matrix_free(&tmp1);
	nv_matrix_free(&tmp2);
}
