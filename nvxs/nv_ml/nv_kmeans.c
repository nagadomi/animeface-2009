#include "nv_core.h"
#include "nv_ml.h"
#include "nv_num.h"

// k-means++

// 最小距離クラス選択
static int 
nv_kmeans_bmc(const nv_matrix_t *mat, int mk,
			  const nv_matrix_t *vec, int vm)
{
	int k;
	int min_k = -1;
	float min_dist = FLT_MAX;

	for (k = 0; k < mk; ++k) {
		float dist = nv_euclidean2(mat, k, vec, vm);
		if (dist < min_dist) {
			min_dist = dist;
			min_k = k;
		}
	}

	return min_k;
}


// K-Means++初期値選択
static void 
nv_kmeans_init(nv_matrix_t *means, int k,
			   const nv_matrix_t *data)
{
	int m, c;
	nv_matrix_t *min_dists = nv_matrix_alloc(1, data->m);
	int rnd = (int)((data->m - 1) * nv_rand());
	float potential = 0.0f;
	int local_tries = 2 + (int)log(data->m);

	// 1つ目
	nv_vector_copy(means, 0, data, rnd);
	for (m = 0; m < data->m; ++m) {
		float dist = nv_euclidean2(means, 0, data, m);
		NV_MAT_V(min_dists, m, 0) = dist;
		potential += dist;
	}
	for (c = 1; c < k; ++c) {
		float min_potential = FLT_MAX;
		int best_index = -1;
		int i;
		for (i = 0; i < local_tries; ++i) {
			float new_potential;
			float rand_val = potential * nv_rand();

			for (m = 0; m < data->m - 1; ++m) {
				if (rand_val <= NV_MAT_V(min_dists, m, 0)) {
					break;
				} else {
					rand_val -= NV_MAT_V(min_dists, m, 0);
				}
			}
			new_potential = 0.0f;
			for (i = 0; i < data->m; ++i) {
				new_potential += min(
					nv_euclidean2(data, m, data, i),
					NV_MAT_V(min_dists, i, 0)
				);
			}
			if (new_potential < min_potential) {
				min_potential = new_potential;
				best_index = m;
			}
		}
		nv_vector_copy(means, c, data, best_index);
		for (m = 0; m < data->m; ++m) {
			NV_MAT_V(min_dists, m, 0) = min(
				nv_euclidean2(data, m, data, best_index),
				NV_MAT_V(min_dists, m, 0)
			);
		}
		potential = min_potential;
	}

	nv_matrix_free(&min_dists);
}

// 最長距離での初期値選択
static void 
nv_kmeans_init_max(nv_matrix_t *means, int k, const nv_matrix_t *data)
{
	int m, c;
	nv_matrix_t *min_dists = nv_matrix_alloc(1, data->m);
	int rnd = (int)((data->m - 1) * nv_rand());

	// 1つ目
	nv_vector_copy(means, 0, data, rnd);
	for (m = 0; m < data->m; ++m) {
		NV_MAT_V(min_dists, m, 0) = nv_euclidean2(means, 0, data, m);
	}
	// 最短距離クラスから一番遠い要素を選択
	for (c = 1; c < k; ++c) {
		float max_dist = -FLT_MAX;
		int max_index = -1;

		for (m = 0; m < data->m; ++m) {
			float dist = min(nv_euclidean2(means, c - 1, data, m), NV_MAT_V(min_dists, m, 0));
			if (dist > max_dist) {
				max_dist = dist;
				max_index = m;
			}
			NV_MAT_V(min_dists, m, 0) = dist;
		}
		nv_vector_copy(means, c, data, max_index);
	}

	nv_matrix_free(&min_dists);
}

// ランダム初期値選択 カブル
static void 
nv_kmeans_init_rand(nv_matrix_t *means, int k, const nv_matrix_t *data)
{
	int c;

	for (c = 0; c < k; ++c) {
		int rnd = (int)((data->m - 1) * nv_rand());
		nv_vector_copy(means, c, data, rnd);
	}
}

int 
nv_kmeans(nv_matrix_t *means,  // k
		  nv_matrix_t *count,  // k
		  nv_matrix_t *labels, // data->m
		  const nv_matrix_t *data,
		  const int k,
		  const int max_epoch)
{
	int m, n, c;
	int processing = 1;
	int converge, epoch;

	nv_matrix_t *old_labels = nv_matrix_alloc(1, data->m);
	nv_matrix_t *sum = nv_matrix_alloc(data->n, k);

	assert(means->n == data->n);
	assert(means->m >= k);
	assert(labels->m >= data->m);

	// 初期値選択
	nv_kmeans_init(means, k, data);

	for (m = 0; m < old_labels->m; ++m) {
		NV_MAT_V(old_labels, m, 0) = -1.0f;
	}

	epoch = 0;
	do {
		nv_matrix_zero(count);
		nv_matrix_zero(sum);
		
		for (m = 0; m < data->m; ++m) {
			int label = nv_kmeans_bmc(means, k, data, m);
			// ラベル決定
			NV_MAT_V(labels, m, 0) = (float)label;
			// カウント
			NV_MAT_V(count, label, 0) += 1.0f;
			// ベクトル合計
			for (n = 0; n < means->n; ++n) {
				NV_MAT_V(sum, label, n) += NV_MAT_V(data, m, n);
			}
		}
		++epoch;

		// 終了判定
		converge = 1;
		for (m = 0; m < data->m; ++m) {
			if (NV_MAT_V(labels, m, 0) != NV_MAT_V(old_labels, m, 0)) {
				converge = 0;
				break;
			}
		}

		if (converge) {
			// 終了
			processing = 0;
		} else {
			// ラベル更新
			nv_matrix_copy(old_labels, 0, labels, 0, old_labels->m);

			// 中央値計算
			for (c = 0; c < k; ++c) {
				if (NV_MAT_V(count, c, 0) != 0.0f) {
					float factor = 1.0f / NV_MAT_V(count, c, 0);
					for (n = 0; n < means->n; ++n) {
						NV_MAT_V(means, c, n) = NV_MAT_V(sum, c, n) * factor;
					}
				}
			}

			// 最大試行回数判定
			if (max_epoch != 0
				&& epoch >= max_epoch)
			{
				// 終了
				processing = 0;
			}
		}
	} while (processing);

	nv_matrix_free(&old_labels);
	nv_matrix_free(&sum);

	return k;
}

