// 
#include <ruby.h>
#include <nv_core.h>
#include <nv_ip.h>
#include <nv_ml.h>
#include <nv_face.h>

static VALUE cMagickPixel;
static VALUE cMagick;

#define NV_MAX_FACE 4096
#define NV_BIT_SCALE(bits)  (1.0f / (powf(2.0f, bits) - 1.0f))
#define NV_RMAGICK_DATA_TO_8BIT_FLOAT(x, depth) ((float)(NUM2DBL(x)) * NV_BIT_SCALE(depth) * 255.0f)

#define NV_ANIMEFACE_WINDOW_SIZE  42.592f
#define NV_ANIMEFACE_STEP         4.0f
#define NV_ANIMEFACE_SCALE_FACTOR 1.095f

static void
nv_conv_imager2nv(nv_matrix_t *bgr, nv_matrix_t *gray,
				  VALUE im)
{
	int y, x,  xsize, ysize, depth;
	VALUE c;
	ID rb_id_green, rb_id_blue, rb_id_red, rb_id_pixel_color;
	
	xsize = NUM2INT(rb_funcall(im, rb_intern("columns"), 0));
	ysize = NUM2INT(rb_funcall(im, rb_intern("rows"), 0));
	depth = NUM2INT(rb_const_get(cMagick, rb_intern("QuantumDepth")));

	assert(bgr->rows == gray->rows && bgr->rows == ysize);
	assert(bgr->cols == gray->cols && bgr->cols == xsize);
	
	rb_id_green = rb_intern("green");
	rb_id_blue = rb_intern("blue");
	rb_id_red = rb_intern("red");
	rb_id_pixel_color = rb_intern("pixel_color");

	if (rb_funcall(im, rb_intern("gray?"), 0) == Qfalse) {
		float b, g, r;
		// rgb
		for (y = 0; y < ysize; ++y) {
			for (x = 0; x < xsize; ++x) {
				c = rb_funcall(im, rb_id_pixel_color, 2, INT2FIX(x), INT2FIX(y));
				g = NV_RMAGICK_DATA_TO_8BIT_FLOAT(rb_funcall(c, rb_id_green, 0), depth);
				b = NV_RMAGICK_DATA_TO_8BIT_FLOAT(rb_funcall(c, rb_id_blue, 0), depth);
				r = NV_RMAGICK_DATA_TO_8BIT_FLOAT(rb_funcall(c, rb_id_red, 0), depth);

				NV_MAT3D_V(bgr, y, x, NV_CH_B) = b;
				NV_MAT3D_V(bgr, y, x, NV_CH_G) = g;
				NV_MAT3D_V(bgr, y, x, NV_CH_R) = r;
			}
		}
		nv_gray(gray, bgr);
	} else {
		// gray
		float g;
		for (y = 0; y < ysize; ++y) {
			for (x = 0; x < xsize; ++x) {
				c = rb_funcall(im, rb_id_pixel_color, 2, INT2FIX(x), INT2FIX(y));
				g = NV_RMAGICK_DATA_TO_8BIT_FLOAT(rb_funcall(c, rb_id_green, 0), depth);
				NV_MAT3D_V(gray, y, x, 0) = g;
				NV_MAT3D_V(bgr, y, x, NV_CH_B) = g;
				NV_MAT3D_V(bgr, y, x, NV_CH_G) = g;
				NV_MAT3D_V(bgr, y, x, NV_CH_R) = g;
			}
		}
	}
}

VALUE detect(VALUE im,
			 float min_window_size,
			 float step, float scale_factor)
{
	static const nv_mlp_t *detector_mlp = &nv_face_mlp_face_00;
	static const nv_mlp_t *face_mlp[] = {
		&nv_face_mlp_face_01,
		&nv_face_mlp_face_02,
		NULL
	};
	static const nv_mlp_t *dir_mlp = &nv_face_mlp_dir;
	static const nv_mlp_t *parts_mlp = &nv_face_mlp_parts;
	
	VALUE results = rb_ary_new();
	VALUE facehash, tmphash1, tmphash2, colors, pixel;
	
	int i, j, nface;
	nv_rect_t image_size;
	nv_face_position_t face_pos[NV_MAX_FACE];
	int xsize = NUM2INT(rb_funcall(im, rb_intern("columns"), 0));
	int ysize = NUM2INT(rb_funcall(im, rb_intern("rows"), 0));
	
	nv_matrix_t *bgr = nv_matrix3d_alloc(3, ysize, xsize);
	nv_matrix_t *gray = nv_matrix3d_alloc(1, ysize, xsize);
	nv_matrix_t *smooth = nv_matrix3d_alloc(1, ysize, xsize);
	nv_matrix_t *edge = nv_matrix3d_alloc(1, ysize, xsize);
	nv_matrix_t *gray_integral = nv_matrix3d_alloc(1, ysize + 1, xsize + 1);
	nv_matrix_t *edge_integral = nv_matrix3d_alloc(1, ysize + 1, xsize + 1);

	// initialize
	
	nv_matrix_zero(bgr);
	nv_matrix_zero(gray);
	nv_matrix_zero(edge);
	nv_matrix_zero(gray_integral);
	nv_matrix_zero(edge_integral);
	
	image_size.x = image_size.y = 0;
	image_size.width = gray->cols;
	image_size.height = gray->rows;
	
	// convert format
	nv_conv_imager2nv(bgr, gray, im);
	// edge
	nv_gaussian5x5(smooth, 0, gray, 0);
	nv_laplacian1(edge, smooth, 4.0f);
	// integral
	nv_integral(gray_integral, gray, 0);
	nv_integral(edge_integral, edge, 0);
	
	// detect face
	nface = nv_face_detect(face_pos, NV_MAX_FACE,
						   gray_integral, edge_integral, &image_size,
						   dir_mlp,
						   detector_mlp, face_mlp, 2,
						   parts_mlp,
						   step, scale_factor, min_window_size
		);
	// analyze face 
	for (i = 0; i < nface; ++i) {
		nv_face_feature_t face_feature = {0};
		nv_face_analyze(&face_feature, &face_pos[i], bgr);
		
		/*
		 * likelihood => ,
		 * face => {x=>,y=>,width=>,height=>, skin_color =>, hair_color=>, },
		 * eyes => {left} => { x=>,y=>,width=>,height=> , color => (,,,,)},
		 * eyes => {right} => { x=>,y=>,width=>,height=> , color => (,,,,)}, 
		 * nose => {x=>, y=>, width=>1, height=>1},
		 * mouse => { x=>, y=>, width=>, height=>},
		 * chin => { x=>, y =>, width=>1, height=>1}
		 */
		
		facehash = rb_hash_new();
		
		// likelihood
		rb_hash_aset(facehash, rb_str_new2("likelihood"), 
					 rb_float_new(face_pos[i].likelihood));
		
		// face.x,y,w,h,skin_color(b,g,r),hair_color(b,g,r)
		tmphash1 = rb_hash_new();
		rb_hash_aset(tmphash1, rb_str_new2("x"), INT2FIX(face_pos[i].face.x));
		rb_hash_aset(tmphash1, rb_str_new2("y"), INT2FIX(face_pos[i].face.y));
		rb_hash_aset(tmphash1, rb_str_new2("width"), INT2FIX(face_pos[i].face.width));
		rb_hash_aset(tmphash1, rb_str_new2("height"), INT2FIX(face_pos[i].face.height));
		rb_hash_aset(facehash, rb_str_new2("face"), tmphash1);
		
		pixel = rb_funcall(cMagickPixel, rb_intern("new"), 3,
						   INT2FIX(face_feature.skin_bgr.v[2]),
						   INT2FIX(face_feature.skin_bgr.v[1]),
						   INT2FIX(face_feature.skin_bgr.v[0]));
		rb_hash_aset(facehash, rb_str_new2("skin_color"), pixel);
		
		pixel = rb_funcall(cMagickPixel, rb_intern("new"), 3,
						   INT2FIX(face_feature.hair_bgr.v[2]),
						   INT2FIX(face_feature.hair_bgr.v[1]),
						   INT2FIX(face_feature.hair_bgr.v[0]));
		rb_hash_aset(facehash, rb_str_new2("hair_color"), pixel);
		
		// left_eye x,y,w,h,color[4](b,g,r)
		tmphash1 = rb_hash_new();
		tmphash2 = rb_hash_new();	/* left */
		rb_hash_aset(tmphash2, rb_str_new2("x"), INT2FIX(face_pos[i].left_eye.x));
		rb_hash_aset(tmphash2, rb_str_new2("y"), INT2FIX(face_pos[i].left_eye.y));
		rb_hash_aset(tmphash2, rb_str_new2("width"), INT2FIX(face_pos[i].left_eye.width));
		rb_hash_aset(tmphash2, rb_str_new2("height"), INT2FIX(face_pos[i].left_eye.height));
		
		colors = rb_ary_new();
		for (j = 0; j < 4; ++j) {
			pixel = rb_funcall(cMagickPixel, rb_intern("new"), 3,
							   INT2FIX(face_feature.left_eye_bgr[j].v[2]),
							   INT2FIX(face_feature.left_eye_bgr[j].v[1]),
							   INT2FIX(face_feature.left_eye_bgr[j].v[0]));
			rb_ary_push(colors, pixel);
		}
		rb_hash_aset(tmphash2, rb_str_new2("colors"), colors);
		rb_hash_aset(tmphash1, rb_str_new2("left"), tmphash2);
		
		// right_eye x,y,w,h,color[4](b,g,r)
		tmphash2 = rb_hash_new();
		rb_hash_aset(tmphash2, rb_str_new2("x"), INT2FIX(face_pos[i].right_eye.x));
		rb_hash_aset(tmphash2, rb_str_new2("y"), INT2FIX(face_pos[i].right_eye.y));
		rb_hash_aset(tmphash2, rb_str_new2("width"), INT2FIX(face_pos[i].right_eye.width));
		rb_hash_aset(tmphash2, rb_str_new2("height"), INT2FIX(face_pos[i].right_eye.height));
		
		colors = rb_ary_new();
		for (j = 0; j < 4; ++j) {
			pixel = rb_funcall(cMagickPixel, rb_intern("new"), 3,
							   INT2FIX(face_feature.right_eye_bgr[j].v[2]),
							   INT2FIX(face_feature.right_eye_bgr[j].v[1]),
							   INT2FIX(face_feature.right_eye_bgr[j].v[0]));
			rb_ary_push(colors, pixel);
		}
		rb_hash_aset(tmphash2, rb_str_new2("colors"), colors);
		rb_hash_aset(tmphash1, rb_str_new2("right"), tmphash2);
		
		rb_hash_aset(facehash, rb_str_new2("eyes"), tmphash1);
		
		// nose x,y
		tmphash1 = rb_hash_new();
		rb_hash_aset(tmphash1, rb_str_new2("x"), INT2FIX(face_pos[i].nose.x));
		rb_hash_aset(tmphash1, rb_str_new2("y"), INT2FIX(face_pos[i].nose.y));
		rb_hash_aset(facehash, rb_str_new2("nose"), tmphash1);
		
		// mouth x,y,w,h
		tmphash1 = rb_hash_new();
		rb_hash_aset(tmphash1, rb_str_new2("x"), INT2FIX(face_pos[i].mouth.x));
		rb_hash_aset(tmphash1, rb_str_new2("y"), INT2FIX(face_pos[i].mouth.y));
		rb_hash_aset(tmphash1, rb_str_new2("width"), INT2FIX(face_pos[i].mouth.width));
		rb_hash_aset(tmphash1, rb_str_new2("height"), INT2FIX(face_pos[i].mouth.height));
		
		// chin x, y
		tmphash1 = rb_hash_new();
		rb_hash_aset(tmphash1, rb_str_new2("x"), INT2FIX(face_pos[i].chin.x));
		rb_hash_aset(tmphash1, rb_str_new2("y"), INT2FIX(face_pos[i].chin.y));
		rb_hash_aset(facehash, rb_str_new2("chin"), tmphash1);
		
		rb_ary_push(results, facehash);
	}
	nv_matrix_free(&bgr);
	nv_matrix_free(&gray);
	nv_matrix_free(&smooth);
	nv_matrix_free(&edge);
	nv_matrix_free(&gray_integral);
	nv_matrix_free(&edge_integral);
	
	return results;
}

static
float get_option(VALUE hash, const char *key, float defvalue)
{
	VALUE v = rb_hash_aref(hash, rb_str_new2(key));
	if (v == Qnil) {
		v = rb_hash_aref(hash, ID2SYM(rb_intern(key)));
	}
	if (v == Qnil) {
		return defvalue;
	}

	return (float)NUM2DBL(v);
}

static VALUE
wrap_detect(int argc, VALUE *argv)
{
	float min_window_size, step, scale_factor;
	VALUE im, options;
	
	if (rb_scan_args(argc, argv, "11", &im, &options) == 2) {
		Check_Type(options, T_HASH);
		
		min_window_size = get_option(options, "min_window_size", NV_ANIMEFACE_WINDOW_SIZE);
		step = get_option(options, "step", NV_ANIMEFACE_STEP);
		scale_factor = get_option(options, "scale_factor", NV_ANIMEFACE_SCALE_FACTOR);
	} else {
		min_window_size = NV_ANIMEFACE_WINDOW_SIZE;
		step = NV_ANIMEFACE_STEP;
		scale_factor = NV_ANIMEFACE_SCALE_FACTOR;
	}

	return detect(im, min_window_size, step, scale_factor);
}

void Init_AnimeFace()
{
  VALUE module;
  
  cMagick = rb_const_get(rb_cObject, rb_intern("Magick"));
  cMagickPixel = rb_const_get(cMagick, rb_intern("Pixel"));
  
  module = rb_define_module("AnimeFace");
  rb_define_module_function(module, "detect", wrap_detect, -1);
}
