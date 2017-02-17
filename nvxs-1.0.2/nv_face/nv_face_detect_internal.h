#ifndef __NV_FACE_DETECTION_INTERNAL_H
#define __NV_FACE_DETECTION_INTERNAL_H

typedef struct candidate {
	float x1, y1, x2, y2;
	double z;
	int flag;
	nv_matrix_t *parts;
} nv_candidate;

static int nv_candidate_cmp(const void *p1, const void *p2)
{
	nv_candidate *e1 = (nv_candidate *)p1;
	nv_candidate *e2 = (nv_candidate *)p2;
	float a1 = (e1->x2 - e1->x1) * (e1->y2 - e1->y1);
	float a2 = (e2->x2 - e2->x1) * (e2->y2 - e2->y1);

	if (a1 > a2) {
		return 1;
	} else if (a1 < a2) {
		return -1;
	}
	return 0;
}


static int nv_is_face_edge(int window, float scale, float area) 
{
	float v = area / (255.0f * window * window) * scale * scale * 0.5f;
	if (window < 84.0f) {
		if (0.1f < v && v < window * 0.012f - 0.2f) {
			return 1;
		}
	} else {
		if ((window - 64.0f) * 0.005f < v && v < window * 0.012f - 0.2f) {
			return 1;
		}
	}
	return 0;
}

// # ñ êœÇ≈É\Å[ÉgÇ≥ÇÍÇΩèáÇ…óàÇÈ
static 
float nv_rect_intersect(float lx1, float ly1, float lx2, float ly2,
						float rx1, float ry1, float rx2, float ry2)
{
	float ax = lx1;
	float ay = ly1;
	float axl = lx2;
	float ayl = ly2;
	float bx = rx1;
	float by = ry1;
	float bxl = rx2;
	float byl = ry2;
	float area_i = (rx2 - rx1) * (ry2 - ry1);
	float area_j = 0;

	if (ax <= bx
		&& bxl <= axl
		&& ay <= by
		&& byl <= axl)
	{
		area_j = (lx2 - lx1) * (ly2 - ly1);
		return 1.0f;
	}

	if (ax <= bx && axl >= bx) {
		// xé≤Ç©Ç‘ÇË a b
		if (ay <= by && ayl >= by) {
			float sx = bx;
			float sy = by;
			float ex = min(axl, bxl);
			float ey = min(ayl, byl);
			area_j = (ex - sx) * (ey - sy);
		} else if (ay <= byl && ayl >= byl) {
			float sx = bx;
			float sy = max(ay, by);
			float ex = min(axl, bxl);
			float ey = min(ayl, byl);
			area_j = (ex - sx) * (ey - sy);
		}
	} else if (ax >= bx && ax <= bxl) {
		// xé≤Ç©Ç‘ÇË b a
		if (by >= ay && by <= ayl) {
			float sx = ax;
			float sy = by;
			float ex = min(bxl, axl);
			float ey = min(byl, ayl);
			area_j = (ex - sx) * (ey - sy);
		} else if (byl >= ay && byl <= ayl) {
			float sx = ax;
			float sy = ay;
			float ex = min(bxl, axl);
			float ey = min(byl, ayl);
			area_j = (ex - sx) * (ey - sy);
		}
	}
	return area_j / area_i;
}

#endif
