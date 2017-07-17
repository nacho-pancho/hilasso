/**
 * Simplistic linear algebra library, dangerously and unwisely implemented
 * from scratch.  the only advantages are: 100% standalone, implements
 * sparse vector/matrix operatins with a Matlab-compatible structure.  I may
 * later switch implementation of dense functions to CBLAS, which is easily
 * available everywhere, but for now this thing works.
 */
#include "mdls_linalg_types.h"

#ifdef USE_CBLAS
#include "mdls_linalg_fun_blas.h"
#else
#include "mdls_linalg_fun.h"
#endif
