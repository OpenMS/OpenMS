// OpenMS addition: TRON kernels for the SSE2 baseline (see tron_kernels.h). Only built into the
// x86-64 baseline (no AVX, no FMA), the build in which ssl.cpp uses the kernels.
#if (defined(__x86_64__) || defined(_M_X64)) && !defined(__AVX__) && !defined(__FMA__)
#define TRON_KERNELS_NS sse2
#include "tron_kernels_impl.h"
#endif
