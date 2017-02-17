#include "nv_core.h"
#include "nv_face_feature.h"

// haar的なもの
// 現在1レベルの解像度しか使用してない

static float nv_face_haarlike_diagonal_filter(int type,
											  const nv_matrix_t *sum,
											  int px, int py,
											  float xscale, float yscale)
{
	int i = 0;
	float p = 0.0f, p1 = 0.0f, p2 = 0.0f;
	int area1 = 0, area2 = 0;
	int ystep = NV_ROUND_INT(yscale);
	int xstep = NV_ROUND_INT(xscale);

	if (type == 1) {
		// |＼|
		for (i = 0; i < 7; ++i) {
			int ppx = px + NV_ROUND_INT((1.0f + i) * xscale);
			int ppy = py + NV_ROUND_INT(i * yscale);
			int eex = px + NV_ROUND_INT(8.0f * xscale);
			int eey = py + NV_ROUND_INT((i + 1) * yscale);

			//printf("p1: %d, %d, %d, %d\n", 1+i,8,i,i+1);

			p1 += NV_INTEGRAL_V(sum, ppx, ppy, eex, eey);
			area1 += (eex - ppx) * (eey - ppy);	
		}
		for (i = 1; i < 8; ++i) {
			int ppx = px;
			int ppy = py + NV_ROUND_INT(i * yscale);
			int eex = px + NV_ROUND_INT(i * xscale);
			int eey = py + NV_ROUND_INT((i + 1) * yscale);

			//printf("p2: %d, %d, %d, %d\n", 0,i,i,i+1);
			p2 += NV_INTEGRAL_V(sum, ppx, ppy, eex, eey);
			area2 += (eex - ppx) * (eey - ppy);	
		}
		p = p1 / (area1 * 255.0f) - p2 / (area2 * 255.0f);
	} else {
		// |/|
		for (i = 0; i < 7; ++i) {
			int ppx = px;
			int ppy = py + NV_ROUND_INT(i * yscale);
			int eex = px + NV_ROUND_INT((7.0f - i) * xscale);
			int eey = py + NV_ROUND_INT((i + 1) * yscale);

			//printf("p1: %d, %d, %d, %d\n", 0, 7-i, i, i+1);
			
			p1 += NV_INTEGRAL_V(sum, ppx, ppy, eex, eey);
			area1 += (eex - ppx) * (eey - ppy);	
		}
		for (i = 1; i < 8; ++i) {
			int ppx = px + NV_ROUND_INT((8.0f - i) * xscale);
			int ppy = py + NV_ROUND_INT(i * yscale);
			int eex = px + NV_ROUND_INT(8.0f * xscale);
			int eey = py + NV_ROUND_INT((i + 1) * yscale);

			//printf("p2: %d, %d, %d, %d\n", 8-i, 8, i, i+1);

			p2 += NV_INTEGRAL_V(sum, ppx, ppy, eex, eey);
			area2 += (eex - ppx) * (eey - ppy);	
		}
		p = p1 / (area1 * 255.0f) - p2 / (area2 * 255.0f);
	}

	return p;
}


void nv_face_haarlike(nv_face_haarlike_normalize_e normalize_type,
					 nv_matrix_t *feature, 
					 int feature_m,
					 const nv_matrix_t *sum,
					 int x, int y, int width, int height)
{
	int ix, iy, n;
	float v, vmax, vmin;
	float xscale = width / 32.0f;
	float yscale = height / 32.0f;
	float ystep = yscale;
	float xstep = xscale;
	int hystep = (32 - 8) / 2 * 8;
	int sy = NV_ROUND_INT(4.0f * ystep);
	int sx = NV_ROUND_INT(4.0f * xstep);
	int hy, hx;

	nv_vector_zero(feature, feature_m);

	// level1
#ifdef _OPENMP
//#pragma omp parallel for private(ix)
#endif
	for (iy = 0, hy = 0; iy < 32-8; iy += 2, ++hy) {
		int py = y + NV_ROUND_INT(ystep * iy);
		int ey = py + NV_ROUND_INT(8.0f * ystep);
		const float pty = (ey - py) * 255.0f;
		for (ix = 0, hx = 0; ix < 32-8; ix += 2, ++hx) {
			int px = x + NV_ROUND_INT(xstep * ix);
			int ex = px + NV_ROUND_INT(8.0f * xstep);
			float p1, p2, area, ptx;

			// 全エリア
			area = NV_MAT3D_V(sum, ey, ex, 0)
				- NV_MAT3D_V(sum, ey, px, 0)
				- (NV_MAT3D_V(sum, py, ex, 0) - NV_MAT3D_V(sum, py, px, 0));
			
			// 1
			// [+]
			// [-]
			p1 = NV_MAT3D_V(sum, py + sy, ex, 0)
				   - NV_MAT3D_V(sum, py + sy, px, 0)
				   - (NV_MAT3D_V(sum, py, ex, 0) - NV_MAT3D_V(sum, py, px, 0));
			p2 = area - p1;
			ptx = (ex - px) * 255.0f;
			p1 /= ((py + sy) - py) * ptx;
			p2 /= (ey - (py + sy)) * ptx;
			if (p1 > p2) {
				NV_MAT_V(feature, feature_m, hy * hystep + hx * 8 + 0) = p1 - p2;
			} else {
				NV_MAT_V(feature, feature_m, hy * hystep + hx * 8 + 1) = p2 - p1;
			}

			// 2
			// [+][-]
			p1 = NV_MAT3D_V(sum, ey, px + sx, 0)
				- NV_MAT3D_V(sum, ey, px, 0)
				- (NV_MAT3D_V(sum, py, px + sx, 0) - NV_MAT3D_V(sum, py, px, 0));
			p2 = area - p1;
			p1 /= ((px + sx) - px) * pty;
			p2 /= (ex - (px + sx)) * pty;
			if (p1 > p2) {
				NV_MAT_V(feature, feature_m, hy * hystep + hx * 8 + 2) = p1 - p2;
			} else {
				NV_MAT_V(feature, feature_m, hy * hystep + hx * 8 + 3) = p2 - p1;
			}

			// 3
			p1 = nv_face_haarlike_diagonal_filter(1, sum, px, py, xscale, yscale);
			if (p1 > 0.0f) {
				NV_MAT_V(feature, feature_m, hy * hystep + hx * 8 + 4) = p1;
			} else {
				NV_MAT_V(feature, feature_m, hy * hystep + hx * 8 + 5) = -p1;
			}

			// 4
			p1 = nv_face_haarlike_diagonal_filter(2, sum, px, py, xscale, yscale);
			if (p1 > 0.0f) {
				NV_MAT_V(feature, feature_m, hy * hystep + hx * 8 + 6) = p1;
			} else {
				NV_MAT_V(feature, feature_m, hy * hystep + hx * 8 + 7) = -p1;
			}
		}
	}

	// 正規化
	switch (normalize_type) {
	case NV_NORMALIZE_MAX:
		// 最大値=1.0
		vmax = 0.0f;
		vmin = FLT_MAX;
		for (n = 0; n < feature->n; ++n) {
			if (NV_MAT_V(feature, feature_m, n) > vmax) {
				vmax = NV_MAT_V(feature, feature_m, n);
			}
			if (NV_MAT_V(feature, feature_m, n) != 0.0f
				&& NV_MAT_V(feature, feature_m, n) < vmin) 
			{
				vmin = NV_MAT_V(feature, feature_m, n);
			}
		}
		if (vmax != 0.0f && vmax > vmin) {
			v = 1.0f / (vmax - vmin);
			for (n = 0; n < feature->n; ++n) {
				if (NV_MAT_V(feature, feature_m, n) != 0.0f) {
					NV_MAT_V(feature, feature_m, n) = (NV_MAT_V(feature, feature_m, n) - vmin) * v;
				}
			}
		}
		break;
	case NV_NORMALIZE_NORM:
		// ベクトル NORM=1.0
		v = 0.0f;
		for (n = 0; n < feature->n; ++n) {
			v += NV_MAT_V(feature, feature_m, n) * NV_MAT_V(feature, feature_m, n);
		}
		if (v != 0.0) {
			v = 1.0f / sqrtf(v);
			for (n = 0; n < feature->n; ++n) {
				NV_MAT_V(feature, feature_m, n) *= v;
			}
		}
		break;
	case NV_NORMALIZE_NONE:
	default:
		break;
	}
}
