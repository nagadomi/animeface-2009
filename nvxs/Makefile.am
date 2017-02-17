EXTRA_DIST 	= autogen.sh CLAPACK
SUBDIRS = CLAPACK
lib_LTLIBRARIES = libnvxs.la
libnvxs_la_SOURCES = \
nv_face/nv_face_analyze.c nv_face/nv_face_haarlike.c nv_face/nv_face_detect.c nv_face/nv_face_mlp_static.c \
nv_core/nv_core_matrix.c nv_core/nv_core_util.c \
nv_ml/nv_gaussian.c nv_ml/nv_kmeans.c nv_ml/nv_mlp.c \
nv_num/nv_cov.c nv_num/nv_eigen.c nv_num/nv_vector.c nv_num/nv_distance.c nv_num/nv_matrix.c \
nv_ip/nv_euclidean_color.c nv_ip/nv_gaussian5x5.c nv_ip/nv_gray.c nv_ip/nv_gray_cpu.c nv_ip/nv_integral.c nv_ip/nv_laplacian.c

libnvxs_la_LIBADD = CLAPACK/liblapack_nvxs.a CLAPACK/libblas_nvxs.a CLAPACK/libtmglib_nvxs.a CLAPACK/F2CLIBS/libf2c.a
libnvxs_la_CFLAGS = -O3 -ffast-math -msse -msse2 -Inv_core -Inv_num -Inv_ip -Inv_ml -ICLAPACK/INCLUDE
libnvxs_la_LDFLAGS = -no-undefined -lm

include_HEADERS = \
nv_face/nv_face.h  nv_face/nv_face_feature.h nv_face/nv_face_mlp_static.h nv_face/nv_face_analyze.h nv_face/nv_face_detect_internal.h nv_face/nv_face_detect.h \
nv_core/nv_config.h nv_core/nv_core.h nv_core/nv_core_matrix.h nv_core/nv_core_util.h nv_core/nv_portable.h \
nv_ip/nv_ip.h nv_ip/nv_ip_euclidean_color.h nv_ip/nv_ip_gaussian.h nv_ip/nv_ip_gray.h nv_ip/nv_ip_integral.h nv_ip/nv_ip_laplacian.h \
nv_num/nv_num.h nv_num/nv_num_cov.h nv_num/nv_num_distance.h nv_num/nv_num_eigen.h nv_num/nv_num_lapack.h nv_num/nv_num_matrix.h nv_num/nv_num_vector.h \
nv_ml/nv_ml.h nv_ml/nv_ml_gaussian.h nv_ml/nv_ml_kmeans.h nv_ml/nv_ml_mlp.h \
CLAPACK/INCLUDE/blaswrap.h CLAPACK/INCLUDE/clapack.h CLAPACK/INCLUDE/f2c.h





