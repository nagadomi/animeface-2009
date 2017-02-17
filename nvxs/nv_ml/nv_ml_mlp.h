#ifndef __NV_ML_MLP_H
#define __NV_ML_MLP_H

#include "nv_core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
	int input;
	int hidden;
	int output;
	nv_matrix_t *input_w;
	nv_matrix_t *hidden_w;
	nv_matrix_t *input_bias;
	nv_matrix_t *hidden_bias;
 } nv_mlp_t;

nv_mlp_t *nv_mlp_alloc(int input, int hidden, int k);
void nv_mlp_free(nv_mlp_t **mlp);

float nv_mlp_sigmoid(float a);
void nv_mlp_init(nv_mlp_t *mlp);
void nv_mlp_gaussian_init(nv_mlp_t *mlp, float var, int height, int width, int zdim);

float nv_mlp_train_ex(nv_mlp_t *mlp,
					 const nv_matrix_t *data, const nv_matrix_t *label,
					 float ir, float hr,
					 int start_epoch, int end_epoch, int max_epoch);
float nv_mlp_train_lex(nv_mlp_t *mlp,
					 const nv_matrix_t *data,
					 const nv_matrix_t *label,
					 const nv_matrix_t *t,
					 float ir, float hr,
					 int start_epoch, int end_epoch, int max_epoch);


float nv_mlp_train(nv_mlp_t *mlp, const nv_matrix_t *data, const nv_matrix_t *label);
int nv_mlp_predict_label(const nv_mlp_t *mlp, const nv_matrix_t *x, int xm);
float nv_mlp_predict(const nv_mlp_t *mlp, const nv_matrix_t *x, int xm, int cls);
float nv_mlp_bagging_predict(const nv_mlp_t **mlp, int nmlp, 
							 const nv_matrix_t *x, int xm, int cls);
double nv_mlp_predict_d(const nv_mlp_t *mlp,
					   const nv_matrix_t *x, int xm, int cls);
double nv_mlp_bagging_predict_d(const nv_mlp_t **mlp, int nmlp, 
							   const nv_matrix_t *x, int xm, int cls);

void nv_mlp_train_regression(
					 nv_mlp_t *mlp,
					 const nv_matrix_t *data,
					 const nv_matrix_t *t,
					 float ir, float hr,
					 int start_epoch, int max_epoch);
void nv_mlp_regression(const nv_mlp_t *mlp, const nv_matrix_t *x, int xm, nv_matrix_t *out, int om);
void nv_mlp_dump_c(FILE *out, const nv_mlp_t *mlp, const char *name, int static_variable);
#ifdef __cplusplus
}
#endif

#endif
