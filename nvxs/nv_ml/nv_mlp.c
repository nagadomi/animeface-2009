#include "nv_core.h"
#include "nv_num.h"
#include "nv_ml_mlp.h"
#include <stdlib.h>


// ‘½‘wƒp[ƒZƒvƒgƒƒ“
// 2 Layer


float nv_mlp_sigmoid(float a)
{
	return 1.0f / (1.0f + expf(-a));
}

double nv_mlp_sigmoid_d(double a)
{
	return 1.0 / (1.0 + exp(-a));
}

nv_mlp_t *nv_mlp_alloc(int input, int hidden, int k)
{
	nv_mlp_t *mlp = (nv_mlp_t *)malloc(sizeof(nv_mlp_t));

	mlp->input = input;
	mlp->hidden = hidden;
	mlp->output = k;

	mlp->input_w = nv_matrix_alloc(input, hidden);
	mlp->hidden_w = nv_matrix_alloc(mlp->input_w->m, k);
	mlp->input_bias = nv_matrix_alloc(1, hidden);
	mlp->hidden_bias = nv_matrix_alloc(1, k);

	return mlp;
}

void nv_mlp_free(nv_mlp_t **mlp)
{
	nv_matrix_free(&(*mlp)->input_w);
	nv_matrix_free(&(*mlp)->input_bias);
	nv_matrix_free(&(*mlp)->hidden_w);
	nv_matrix_free(&(*mlp)->hidden_bias);
	free(*mlp);
	*mlp = NULL;
}

void nv_mlp_dump_c(FILE *out, const nv_mlp_t *mlp, const char *name, int static_variable)
{
	char var_name[4][1024];

	snprintf(var_name[0], sizeof(var_name[0]), "%s_input_w", name);
	nv_matrix_dump_c(out, mlp->input_w, var_name[0], 1);
	snprintf(var_name[1], sizeof(var_name[1]), "%s_hidden_w", name);
	nv_matrix_dump_c(out, mlp->hidden_w, var_name[1], 1);
	snprintf(var_name[2], sizeof(var_name[2]), "%s_input_bias", name);
	nv_matrix_dump_c(out, mlp->input_bias, var_name[2], 1);
	snprintf(var_name[3], sizeof(var_name[3]), "%s_hidden_bias", name);
	nv_matrix_dump_c(out, mlp->hidden_bias, var_name[3], 1);

	fprintf(out, "%snv_mlp_t %s = {\n %d, %d, %d, &%s, &%s, &%s, &%s\n};\n",
		static_variable ? "static ":"",
		name, mlp->input, mlp->hidden, mlp->output,
		var_name[0],  var_name[1],  var_name[2],  var_name[3]);
	fflush(out);
}

// ƒNƒ‰ƒX•ª—Ş

int nv_mlp_predict_label(const nv_mlp_t *mlp, const nv_matrix_t *x, int xm)
{
	int m, n;
	int l = -1;
	float mp = -FLT_MAX, y;

	nv_matrix_t *input_y = nv_matrix_alloc(mlp->input_w->m, 1);
	nv_matrix_t *output_y = nv_matrix_alloc(mlp->output, 1);

	/* ‡“`”d */
	// “ü—Í‘w

#ifdef _OPENMP
#pragma omp parallel for private(y) 
#endif
	for (m = 0; m < mlp->hidden; ++m) {
		y = NV_MAT_V(mlp->input_bias, m, 0);
		y += nv_vector_dot(x, xm, mlp->input_w, m);
		NV_MAT_V(input_y, 0, m) = nv_mlp_sigmoid(y);
	}

	// ‰B‚ê‘w
	// o—Í‘w
	for (m = 0; m < mlp->output; ++m) {
		y = NV_MAT_V(mlp->hidden_bias, m, 0);
		y += nv_vector_dot(input_y, 0, mlp->hidden_w, m);
		NV_MAT_V(output_y, 0, m) = nv_mlp_sigmoid(y);
	}

	l = -1; // nega
	if (output_y->n == 1) {
		if (NV_MAT_V(output_y, 0, 0) > 0.5f) {
			l = 0;
		} else {
			l = 1;
		}
	} else {
		for (n = 0; n < output_y->n; ++n) {
			if (mp < NV_MAT_V(output_y, 0, n)
				&& NV_MAT_V(output_y, 0, n) > 0.5f
				) 
			{
				mp = NV_MAT_V(output_y, 0, n);
				l = n;
			}
		}
	}

	nv_matrix_free(&input_y);
	nv_matrix_free(&output_y);

	return l;
}

double nv_mlp_predict_d(const nv_mlp_t *mlp,
						const nv_matrix_t *x, int xm, int cls)
{
	int m;
	int l = -1;
	float mp = -FLT_MAX, y;
	nv_matrix_t *input_y = nv_matrix_alloc(mlp->input_w->m, 1);
	nv_matrix_t *output_z = nv_matrix_alloc(mlp->output, 1);
	nv_matrix_t *output_y = nv_matrix_alloc(mlp->output, 1);
	double p;

	// “ü—Í‘w
#ifdef _OPENMP
// omp_set_nested(false);
#pragma omp parallel for private(y)
#endif
	for (m = 0; m < mlp->hidden; ++m) {
		y = NV_MAT_V(mlp->input_bias, m, 0);
		y += nv_vector_dot(x, xm, mlp->input_w, m);
		NV_MAT_V(input_y, 0, m) = nv_mlp_sigmoid(y);
	}

	// ‰B‚ê‘w
	// o—Í‘w
	for (m = 0; m < mlp->output; ++m) {
		y = NV_MAT_V(mlp->hidden_bias, m, 0);;
		y += nv_vector_dot(input_y, 0, mlp->hidden_w, m);
		NV_MAT_V(output_z, 0, m) = y;
		NV_MAT_V(output_y, 0, m) = y;
	}
	p = (double)NV_MAT_V(output_z, 0, cls);
	p = nv_mlp_sigmoid_d(p);

	nv_matrix_free(&input_y);
	nv_matrix_free(&output_y);
	nv_matrix_free(&output_z);

	return p;
}


float nv_mlp_predict(const nv_mlp_t *mlp,
					 const nv_matrix_t *x, int xm, int cls)
{
	int m;
	int l = -1;
	float mp = -FLT_MAX, y;
	nv_matrix_t *input_y = nv_matrix_alloc(mlp->input_w->m, 1);
	nv_matrix_t *output_y = nv_matrix_alloc(mlp->output, 1);
	float p;

	// “ü—Í‘w
#ifdef _OPENMP
// omp_set_nested(false);
#pragma omp parallel for private(y)
#endif
	for (m = 0; m < mlp->hidden; ++m) {
		y = NV_MAT_V(mlp->input_bias, m, 0);
		y += nv_vector_dot(x, xm, mlp->input_w, m);
		NV_MAT_V(input_y, 0, m) = nv_mlp_sigmoid(y);
	}

	// ‰B‚ê‘w
	// o—Í‘w
	for (m = 0; m < mlp->output; ++m) {
		y = NV_MAT_V(mlp->hidden_bias, m, 0);;
		y += nv_vector_dot(input_y, 0, mlp->hidden_w, m);
		NV_MAT_V(output_y, 0, m) = nv_mlp_sigmoid(y);
	}
	p = NV_MAT_V(output_y, 0, cls);

	nv_matrix_free(&input_y);
	nv_matrix_free(&output_y);

	return p;
}

double nv_mlp_bagging_predict_d(const nv_mlp_t **mlp, int nmlp, 
							   const nv_matrix_t *x, int xm, int cls)
{
	double p = 0.0f;
	double factor = 1.0 / nmlp;
	int i;
	
	for (i = 0; i < nmlp; ++i) {
		p += factor * nv_mlp_predict_d(mlp[i], x, xm, cls);
	}

	return p;
}

float nv_mlp_bagging_predict(const nv_mlp_t **mlp, int nmlp, 
							 const nv_matrix_t *x, int xm, int cls)
{
	float p = 0.0f;
	float factor = 1.0f / nmlp;
	int i;
	
	for (i = 0; i < nmlp; ++i) {
		p += factor * nv_mlp_predict(mlp[i], x, xm, cls);
	}

	return p;
}

// ”ñüŒ`d‰ñ‹A

void nv_mlp_regression(const nv_mlp_t *mlp, 
					   const nv_matrix_t *x, int xm, nv_matrix_t *out, int om)
{
	int m;
	int l = -1;
	float mp = -FLT_MAX, y;
	nv_matrix_t *input_y = nv_matrix_alloc(mlp->input_w->m, 1);
	nv_matrix_t *hidden_y = nv_matrix_alloc(mlp->hidden_w->m, 1);

	// “ü—Í‘w
#ifdef _OPENMP
#pragma omp parallel for private(y)
#endif
	for (m = 0; m < mlp->input_w->m; ++m) {
		y = NV_MAT_V(mlp->input_bias, m, 0);
		y += nv_vector_dot(x, xm, mlp->input_w, m);
		NV_MAT_V(input_y, 0, m) = nv_mlp_sigmoid(y);
	}

	// ‰B‚ê‘w
	for (m = 0; m < mlp->hidden_w->m; ++m) {
		y = NV_MAT_V(mlp->hidden_bias, m, 0);
		y += nv_vector_dot(input_y, 0, mlp->hidden_w, m);
		NV_MAT_V(hidden_y, 0, m) = y;
	}

	// o—Í‘w
	nv_vector_copy(out, om, hidden_y, 0);

	nv_matrix_free(&input_y);
	nv_matrix_free(&hidden_y);
}

