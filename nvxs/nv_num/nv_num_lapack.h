#ifndef __NV_NUM_LAPACK_H
#define __NV_NUM_LAPACK_H

#ifdef __cplusplus
extern "C" {
#endif

#include "blaswrap.h"
#include "f2c.h"
#include "clapack.h"


// ilaenv
integer ilaenv_(integer *, char *, char *, integer *, integer *, 
				integer *, integer *);

#ifdef __cplusplus
}
#endif
#endif
