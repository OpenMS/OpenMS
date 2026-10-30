// OpenMS addition: TRON kernels with AVX2 and without FMA (see tron_kernels.h); compiled with
// -mavx2 -mno-fma -ffp-contract=off, chosen by ssl.cpp at run time on CPUs with AVX2.
#if defined(__x86_64__) || defined(_M_X64)
#if defined(__AVX2__) && !defined(__FMA__)
#define TRON_KERNELS_NS avx2
#include "tron_kernels_impl.h"
#else
// Not compiled for AVX2 without FMA (MSVC, where CMake adds no /arch:AVX2; a macOS universal build;
// FMA enabled for this file): an empty table, which ssl.cpp never selects. FMA kernels would round
// differently from the reference.
#include "tron_kernels.h"
namespace OpenMS { namespace Internal { namespace Percolator { namespace TronKernels {
const Ops avx2::ops = {"avx2 (not built)", nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
}}}}
#endif
#endif
