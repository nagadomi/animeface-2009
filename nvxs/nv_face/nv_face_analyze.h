#ifndef __NV_FACE_ANALYZE_H
#define __NV_FACE_ANALYZE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	// 眼の縦横比 
	// x/y が大きいほう(斜めの場合片目が縦長になるため)
	float eye_ratio;
	// (眼の上～アゴ)/(顔縦幅-顎から下の長さ)
	float face_ratio;
	// 眼から口の間の平均色
	nv_color_t skin_bgr;
	nv_color_t skin_ec;
	// 眼の上から眼の3倍まで幅を色でクラスタリングしたとき
	// 肌と離れている最大要素のクラスの平均色
	nv_color_t hair_bgr;
	nv_color_t hair_ec;
	// 目の色をクラスタリングした肌を除く上位4色
	nv_color_t left_eye_bgr[4];
	nv_color_t right_eye_bgr[4];
	nv_color_t eye_bgr[4];
	// 目の色をクラスタリングした肌を除く上位4色(euclidean_color)
	nv_color_t left_eye_ec[4];
	nv_color_t right_eye_ec[4];
	nv_color_t eye_ec[4];
} nv_face_feature_t;

void 
nv_face_analyze(nv_face_feature_t *feature,
				const nv_face_position_t *face,
				const nv_matrix_t *img);


#ifdef __cplusplus
}
#endif
#endif