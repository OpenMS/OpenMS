/*
 * OpenMS addition (not part of upstream Percolator): run-time check for AVX2 (CPU support and
 * YMM state enabled by the OS), for the TRON kernels (tron_kernels.h). GCC/Clang and MSVC on x86.
 */
#ifndef PERCOLATOR_CPU_AVX2_H
#define PERCOLATOR_CPU_AVX2_H

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#include <intrin.h>
#include <immintrin.h>
#elif (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#include <cpuid.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
inline bool cpuHasAvx2()
{
  int r[4];
  __cpuid(r, 0);
  if (r[0] < 7) return false;
  __cpuid(r, 1);
  if (!(r[2] & (1 << 27)) || !(r[2] & (1 << 28))) return false;  // OSXSAVE, AVX
  if ((_xgetbv(0) & 6) != 6) return false;                        // XMM and YMM state enabled by the OS
  __cpuidex(r, 7, 0);
  return (r[1] & (1 << 5)) != 0;                                  // AVX2
}
#elif (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
inline bool cpuHasAvx2()
{
  unsigned a, b, c, d;
  if (!__get_cpuid(0, &a, &b, &c, &d) || a < 7) return false;
  __get_cpuid(1, &a, &b, &c, &d);
  if (!(c & (1u << 27)) || !(c & (1u << 28))) return false;  // OSXSAVE, AVX
  unsigned lo, hi;
  __asm__ volatile("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
  if ((lo & 6) != 6) return false;  // XMM and YMM state enabled by the OS
  __cpuid_count(7, 0, a, b, c, d);
  return (b & (1u << 5)) != 0;      // AVX2
}
#else
inline bool cpuHasAvx2() { return false; }
#endif

}}}  // namespace OpenMS::Internal::Percolator

#endif
