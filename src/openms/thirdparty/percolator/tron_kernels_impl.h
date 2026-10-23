/*
 * OpenMS addition (not part of upstream Percolator): implementation of the TRON kernels
 * (tron_kernels.h). Included once by tron_kernels_sse2.cpp and once by tron_kernels_avx2.cpp;
 * TRON_KERNELS_NS names the namespace of the table. Uses AVX2 code when __AVX2__ is defined.
 * No Eigen here: the two TUs are compiled with different instruction sets, and inline functions
 * of a shared header would be merged across them (ODR).
 * Must be compiled without FMA and, for the AVX2 TU, with -ffp-contract=off: every product is
 * rounded before it is added, as in the reference.
 */
#include "tron_kernels.h"

#include <bit>
#include <cstddef>
#include <immintrin.h>

#ifdef __FMA__
#error "the TRON kernels must not be compiled with FMA"
#endif

namespace OpenMS { namespace Internal { namespace Percolator { namespace TronKernels { namespace TRON_KERNELS_NS {

#ifdef __AVX2__
constexpr int NB = 4;  // blocks of 4 rows per dot iteration (independent chains)
#else
constexpr int NB = 8;
#endif
constexpr int kPrefetch = 2;  // software prefetch distance in superblocks of the Hv scatter (0 = off)

#ifndef __AVX2__
// ---------------------------------------------------------------- SSE2

// dots of NB blocks (4*NB rows) starting at row q
static inline void dots_nb(const double* A, int n, int np, const double* s, const double* c2, int q, double* sc)
{
  const double* B = A + (std::ptrdiff_t)(q / 4) * np * 4;
  __m128d lo[NB], hi[NB];
  for (int b = 0; b < NB; ++b) { lo[b] = _mm_setzero_pd(); hi[b] = _mm_setzero_pd(); }
  for (int j = 0; j < n; ++j)
  {
    const __m128d sj = _mm_set1_pd(s[j]);
    for (int b = 0; b < NB; ++b)
    {
      const double* p = B + (std::ptrdiff_t)b * np * 4 + 4 * j;
      lo[b] = _mm_add_pd(lo[b], _mm_mul_pd(_mm_load_pd(p), sj));
      hi[b] = _mm_add_pd(hi[b], _mm_mul_pd(_mm_load_pd(p + 2), sj));
    }
  }
  for (int b = 0; b < NB; ++b)
  {
    _mm_storeu_pd(sc + q + 4 * b, _mm_mul_pd(_mm_loadu_pd(c2 + q + 4 * b), lo[b]));
    _mm_storeu_pd(sc + q + 4 * b + 2, _mm_mul_pd(_mm_loadu_pd(c2 + q + 4 * b + 2), hi[b]));
  }
}

// h(j = 4g .. 4g+3) += the rows of one block (in row order) scaled by sc (4 values)
static inline void scatter_blk(const double* B, int g, const double* scp, __m128d& h0, __m128d& h1)
{
  const __m128d s01 = _mm_loadu_pd(scp), s23 = _mm_loadu_pd(scp + 2);
  const double* p = B + 16 * g;
  // columns 4g, 4g+1 -> h0; 4g+2, 4g+3 -> h1
  __m128d a0 = _mm_mul_pd(_mm_load_pd(p + 0), s01), a1 = _mm_mul_pd(_mm_load_pd(p + 4), s01);
  __m128d a2 = _mm_mul_pd(_mm_load_pd(p + 8), s01), a3 = _mm_mul_pd(_mm_load_pd(p + 12), s01);
  __m128d b0 = _mm_mul_pd(_mm_load_pd(p + 2), s23), b1 = _mm_mul_pd(_mm_load_pd(p + 6), s23);
  __m128d b2 = _mm_mul_pd(_mm_load_pd(p + 10), s23), b3 = _mm_mul_pd(_mm_load_pd(p + 14), s23);
  h0 = _mm_add_pd(h0, _mm_unpacklo_pd(a0, a1)); h1 = _mm_add_pd(h1, _mm_unpacklo_pd(a2, a3));  // row 0
  h0 = _mm_add_pd(h0, _mm_unpackhi_pd(a0, a1)); h1 = _mm_add_pd(h1, _mm_unpackhi_pd(a2, a3));  // row 1
  h0 = _mm_add_pd(h0, _mm_unpacklo_pd(b0, b1)); h1 = _mm_add_pd(h1, _mm_unpacklo_pd(b2, b3));  // row 2
  h0 = _mm_add_pd(h0, _mm_unpackhi_pd(b0, b1)); h1 = _mm_add_pd(h1, _mm_unpackhi_pd(b2, b3));  // row 3
}
static inline void hload(const double* H, int g, __m128d& h0, __m128d& h1) { h0 = _mm_loadu_pd(H + 4 * g); h1 = _mm_loadu_pd(H + 4 * g + 2); }
static inline void hstore(double* H, int g, __m128d h0, __m128d h1) { _mm_storeu_pd(H + 4 * g, h0); _mm_storeu_pd(H + 4 * g + 2, h1); }
typedef __m128d HV[2];

// margins of rows [b, e) of BT4-stored X, b and e multiples of 4; per row Eigen's order
static void margins_bt(const double* XB, int n, int np, const double* w, const double* t, int b, int e, double* margin)
{
  const int n2 = n & ~1;
  int i = b;
  for (; i + 8 <= e; i += 8)
  {
    const double* B0 = XB + (std::ptrdiff_t)(i / 4) * np * 4;
    const double* B1 = B0 + np * 4;
    __m128d e0 = _mm_setzero_pd(), o0 = e0, e1 = e0, o1 = e0, e2 = e0, o2 = e0, e3 = e0, o3 = e0;
    for (int j = 0; j < n2; j += 2)
    {
      const __m128d we = _mm_set1_pd(w[j]), wo = _mm_set1_pd(w[j + 1]);
      const double* p0 = B0 + 4 * j;
      const double* p1 = B1 + 4 * j;
      e0 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(p0), we), e0); e1 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(p0 + 2), we), e1);
      e2 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(p1), we), e2); e3 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(p1 + 2), we), e3);
      o0 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(p0 + 4), wo), o0); o1 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(p0 + 6), wo), o1);
      o2 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(p1 + 4), wo), o2); o3 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(p1 + 6), wo), o3);
    }
    __m128d c0 = _mm_add_pd(e0, o0), c1 = _mm_add_pd(e1, o1), c2 = _mm_add_pd(e2, o2), c3 = _mm_add_pd(e3, o3);
    if (n & 1)
    {
      const __m128d wl = _mm_set1_pd(w[n - 1]);
      const int l = 4 * (n - 1);
      c0 = _mm_add_pd(c0, _mm_mul_pd(_mm_load_pd(B0 + l), wl)); c1 = _mm_add_pd(c1, _mm_mul_pd(_mm_load_pd(B0 + l + 2), wl));
      c2 = _mm_add_pd(c2, _mm_mul_pd(_mm_load_pd(B1 + l), wl)); c3 = _mm_add_pd(c3, _mm_mul_pd(_mm_load_pd(B1 + l + 2), wl));
    }
    const __m128d z = _mm_setzero_pd();
    _mm_storeu_pd(margin + i, _mm_mul_pd(_mm_add_pd(z, c0), _mm_loadu_pd(t + i)));
    _mm_storeu_pd(margin + i + 2, _mm_mul_pd(_mm_add_pd(z, c1), _mm_loadu_pd(t + i + 2)));
    _mm_storeu_pd(margin + i + 4, _mm_mul_pd(_mm_add_pd(z, c2), _mm_loadu_pd(t + i + 4)));
    _mm_storeu_pd(margin + i + 6, _mm_mul_pd(_mm_add_pd(z, c3), _mm_loadu_pd(t + i + 6)));
  }
  for (; i < e; ++i)
  {
    const double* x = XB + (std::ptrdiff_t)(i / 4) * np * 4 + (i & 3);
    double ev = 0.0, od = 0.0;
    for (int j = 0; j < n2; j += 2) { ev = x[4 * j] * w[j] + ev; od = x[4 * j + 4] * w[j + 1] + od; }
    double cc = ev + od;
    if (n & 1) cc += x[4 * (n - 1)] * w[n - 1];
    margin[i] = (0.0 + cc) * t[i];
  }
}
#else
// ---------------------------------------------------------------- AVX2 (no FMA)

static inline void dots_nb(const double* A, int n, int np, const double* s, const double* c2, int q, double* sc)
{
  const double* B = A + (std::ptrdiff_t)(q / 4) * np * 4;
  __m256d acc[NB];
  for (int b = 0; b < NB; ++b) acc[b] = _mm256_setzero_pd();
  for (int j = 0; j < n; ++j)
  {
    const __m256d sj = _mm256_broadcast_sd(s + j);
    for (int b = 0; b < NB; ++b) acc[b] = _mm256_add_pd(acc[b], _mm256_mul_pd(_mm256_load_pd(B + (std::ptrdiff_t)b * np * 4 + 4 * j), sj));
  }
  for (int b = 0; b < NB; ++b) _mm256_storeu_pd(sc + q + 4 * b, _mm256_mul_pd(_mm256_loadu_pd(c2 + q + 4 * b), acc[b]));
}

static inline void scatter_blk(const double* B, int g, const double* scp, __m256d& h, __m256d&)
{
  const __m256d sv = _mm256_loadu_pd(scp);
  const double* p = B + 16 * g;
  const __m256d a0 = _mm256_mul_pd(_mm256_load_pd(p + 0), sv), a1 = _mm256_mul_pd(_mm256_load_pd(p + 4), sv);
  const __m256d a2 = _mm256_mul_pd(_mm256_load_pd(p + 8), sv), a3 = _mm256_mul_pd(_mm256_load_pd(p + 12), sv);
  const __m256d t0 = _mm256_unpacklo_pd(a0, a1), t1 = _mm256_unpackhi_pd(a0, a1);  // (r0j0 r0j1 r2j0 r2j1) (r1.. r3..)
  const __m256d t2 = _mm256_unpacklo_pd(a2, a3), t3 = _mm256_unpackhi_pd(a2, a3);  // (r0j2 r0j3 r2j2 r2j3) (r1.. r3..)
  h = _mm256_add_pd(h, _mm256_permute2f128_pd(t0, t2, 0x20));  // row 0
  h = _mm256_add_pd(h, _mm256_permute2f128_pd(t1, t3, 0x20));  // row 1
  h = _mm256_add_pd(h, _mm256_permute2f128_pd(t0, t2, 0x31));  // row 2
  h = _mm256_add_pd(h, _mm256_permute2f128_pd(t1, t3, 0x31));  // row 3
}
static inline void hload(const double* H, int g, __m256d& h0, __m256d&) { h0 = _mm256_loadu_pd(H + 4 * g); }
static inline void hstore(double* H, int g, __m256d h0, __m256d) { _mm256_storeu_pd(H + 4 * g, h0); }
typedef __m256d HV[2];

// margins of rows [b, e) of BT4-stored X, b and e multiples of 4; lane = row, even/odd accumulators per row
static void margins_bt(const double* XB, int n, int np, const double* w, const double* t, int b, int e, double* margin)
{
  const int n2 = n & ~1;
  int i = b;
  for (; i + 16 <= e; i += 16)
  {
    const double* B = XB + (std::ptrdiff_t)(i / 4) * np * 4;
    const std::ptrdiff_t bs = np * 4;
    __m256d e0 = _mm256_setzero_pd(), o0 = e0, e1 = e0, o1 = e0, e2 = e0, o2 = e0, e3 = e0, o3 = e0;
    for (int j = 0; j < n2; j += 2)
    {
      const __m256d we = _mm256_broadcast_sd(w + j), wo = _mm256_broadcast_sd(w + j + 1);
      const double* p = B + 4 * j;
      e0 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p), we), e0); o0 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p + 4), wo), o0);
      e1 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p + bs), we), e1); o1 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p + bs + 4), wo), o1);
      e2 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p + 2 * bs), we), e2); o2 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p + 2 * bs + 4), wo), o2);
      e3 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p + 3 * bs), we), e3); o3 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p + 3 * bs + 4), wo), o3);
    }
    __m256d c0 = _mm256_add_pd(e0, o0), c1 = _mm256_add_pd(e1, o1), c2 = _mm256_add_pd(e2, o2), c3 = _mm256_add_pd(e3, o3);
    if (n & 1)
    {
      const __m256d wl = _mm256_broadcast_sd(w + n - 1);
      const double* p = B + 4 * (n - 1);
      c0 = _mm256_add_pd(c0, _mm256_mul_pd(_mm256_load_pd(p), wl)); c1 = _mm256_add_pd(c1, _mm256_mul_pd(_mm256_load_pd(p + bs), wl));
      c2 = _mm256_add_pd(c2, _mm256_mul_pd(_mm256_load_pd(p + 2 * bs), wl)); c3 = _mm256_add_pd(c3, _mm256_mul_pd(_mm256_load_pd(p + 3 * bs), wl));
    }
    const __m256d z = _mm256_setzero_pd();
    _mm256_storeu_pd(margin + i, _mm256_mul_pd(_mm256_add_pd(z, c0), _mm256_loadu_pd(t + i)));
    _mm256_storeu_pd(margin + i + 4, _mm256_mul_pd(_mm256_add_pd(z, c1), _mm256_loadu_pd(t + i + 4)));
    _mm256_storeu_pd(margin + i + 8, _mm256_mul_pd(_mm256_add_pd(z, c2), _mm256_loadu_pd(t + i + 8)));
    _mm256_storeu_pd(margin + i + 12, _mm256_mul_pd(_mm256_add_pd(z, c3), _mm256_loadu_pd(t + i + 12)));
  }
  for (; i + 4 <= e; i += 4)
  {
    const double* p0 = XB + (std::ptrdiff_t)(i / 4) * np * 4;
    __m256d e0 = _mm256_setzero_pd(), o0 = e0;
    for (int j = 0; j < n2; j += 2)
    {
      e0 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p0 + 4 * j), _mm256_broadcast_sd(w + j)), e0);
      o0 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(p0 + 4 * j + 4), _mm256_broadcast_sd(w + j + 1)), o0);
    }
    __m256d c0 = _mm256_add_pd(e0, o0);
    if (n & 1) c0 = _mm256_add_pd(c0, _mm256_mul_pd(_mm256_load_pd(p0 + 4 * (n - 1)), _mm256_broadcast_sd(w + n - 1)));
    _mm256_storeu_pd(margin + i, _mm256_mul_pd(_mm256_add_pd(_mm256_setzero_pd(), c0), _mm256_loadu_pd(t + i)));
  }
}
#endif

// fun on BT4-stored X: margins (Eigen's order), active rows in order (idx), f += C*d*d serially in row order,
// and (A != nullptr) the active rows packed into A while the chunk of X is still in L1. Returns #active.
static int fun_bt(const double* XB, int n, int np, const double* w, const double* t, const double* C, int m,
                  double* margin, int* idx, double* A, int capA, double* f, int* packed)
{
  const int mp = (m + 31) & ~31;
  const std::ptrdiff_t bs = (std::ptrdiff_t)np * 4;
  double fs = *f;
  int c = 0;
  for (int i = 0; i < mp; i += 32)
  {
    margins_bt(XB, n, np, w, t, i, i + 32, margin);
    unsigned mask = 0;
    const int lim = m - i < 32 ? m - i : 32;
    for (int r = 0; r < lim; ++r) mask |= (unsigned)(margin[i + r] < 1.0) << r;
    while (mask)
    {
      const int r = std::countr_zero(mask);
      mask &= mask - 1;
      const int row = i + r;
      const double d = 1.0 - margin[row];
      fs += C[row] * d * d;
      if (A && c >= capA) A = nullptr;  // buffer too small: grad packs instead
      if (A)
      {
        const double* src = XB + (std::ptrdiff_t)(row >> 2) * bs + (row & 3);
        double* dst = A + (std::ptrdiff_t)(c >> 2) * bs + (c & 3);
        for (int j = 0; j < np; ++j) dst[4 * j] = src[4 * j];
      }
      idx[c++] = row;
    }
  }
  if (A && ((c + 31) & ~31) > capA) A = nullptr;
  *packed = A != nullptr;
  if (A)
  {
    const int nq = (c + 31) & ~31;
    for (int q = c; q < nq; ++q)
    {
      double* dst = A + (std::ptrdiff_t)(q >> 2) * bs + (q & 3);
      for (int j = 0; j < np; ++j) dst[4 * j] = 0.0;
    }
  }
  *f = fs;
  return c;
}

static void build_bt(const double* const* rows, int m, int nf, int np, double* XB)
{
  const int mp = (m + 31) & ~31;
  const std::ptrdiff_t bs = (std::ptrdiff_t)np * 4;
  for (int b = 0; b < mp / 4; ++b)
  {
    double* B = XB + (std::ptrdiff_t)b * bs;
    if (4 * b + 3 < m)
    {
      const double* r0 = rows[4 * b];
      const double* r1 = rows[4 * b + 1];
      const double* r2 = rows[4 * b + 2];
      const double* r3 = rows[4 * b + 3];
      int j = 0;
#ifdef __AVX2__
      for (; j + 4 <= nf; j += 4)
      {
        const __m256d a0 = _mm256_loadu_pd(r0 + j), a1 = _mm256_loadu_pd(r1 + j), a2 = _mm256_loadu_pd(r2 + j), a3 = _mm256_loadu_pd(r3 + j);
        const __m256d t0 = _mm256_unpacklo_pd(a0, a1), t1 = _mm256_unpackhi_pd(a0, a1), t2 = _mm256_unpacklo_pd(a2, a3), t3 = _mm256_unpackhi_pd(a2, a3);
        _mm256_store_pd(B + 4 * j, _mm256_permute2f128_pd(t0, t2, 0x20)); _mm256_store_pd(B + 4 * j + 4, _mm256_permute2f128_pd(t1, t3, 0x20));
        _mm256_store_pd(B + 4 * j + 8, _mm256_permute2f128_pd(t0, t2, 0x31)); _mm256_store_pd(B + 4 * j + 12, _mm256_permute2f128_pd(t1, t3, 0x31));
      }
#else
      for (; j + 2 <= nf; j += 2)
      {
        const __m128d a0 = _mm_loadu_pd(r0 + j), a1 = _mm_loadu_pd(r1 + j), a2 = _mm_loadu_pd(r2 + j), a3 = _mm_loadu_pd(r3 + j);
        _mm_store_pd(B + 4 * j, _mm_unpacklo_pd(a0, a1)); _mm_store_pd(B + 4 * j + 2, _mm_unpacklo_pd(a2, a3));
        _mm_store_pd(B + 4 * j + 4, _mm_unpackhi_pd(a0, a1)); _mm_store_pd(B + 4 * j + 6, _mm_unpackhi_pd(a2, a3));
      }
#endif
      for (; j < nf; ++j) { B[4 * j] = r0[j]; B[4 * j + 1] = r1[j]; B[4 * j + 2] = r2[j]; B[4 * j + 3] = r3[j]; }
      for (int r = 0; r < 4; ++r) B[4 * nf + r] = 1.0;
      for (std::ptrdiff_t k = 4 * (nf + 1); k < 4 * (std::ptrdiff_t)np; ++k) B[k] = 0.0;
    }
    else
    {
      for (int r = 0; r < 4; ++r)
      {
        const int i = 4 * b + r;
        if (i < m)
        {
          for (int j = 0; j < nf; ++j) B[4 * j + r] = rows[i][j];
          B[4 * nf + r] = 1.0;
          for (int j = nf + 1; j < np; ++j) B[4 * j + r] = 0.0;
        }
        else
          for (int j = 0; j < np; ++j) B[4 * j + r] = 0.0;
      }
    }
  }
}

static void pack_bt(const double* XB, int np, const int* idx, int na, int qb, int qe, double* A)
{
  for (int q = qb; q < qe; ++q)
  {
    double* d = A + (std::ptrdiff_t)(q / 4) * np * 4 + (q & 3);
    if (q < na)
    {
      const int i = idx[q];
      const double* x = XB + (std::ptrdiff_t)(i / 4) * np * 4 + (i & 3);
      for (int j = 0; j < np; ++j) d[4 * j] = x[4 * j];
    }
    else
      for (int j = 0; j < np; ++j) d[4 * j] = 0.0;
  }
}

// superblocks of NB blocks outer, column groups inner (independent chains per group)
static void scatter_generic(const double* A, int np, const double* sc, int nq, int gb, int ge, double* out)
{
  for (int j = 4 * gb; j < 4 * ge; ++j) out[j] = 0.0;
  for (int q = 0; q < nq; q += 4 * NB)
  {
    const double* B = A + (std::ptrdiff_t)(q / 4) * np * 4;
    for (int g = gb; g < ge; ++g)
    {
      HV h{};  // AVX2: h[1] unused, but passed by value
      hload(out, g, h[0], h[1]);
      for (int b = 0; b < NB; ++b) scatter_blk(B + (std::ptrdiff_t)b * np * 4, g, sc + q + 4 * b, h[0], h[1]);
      hstore(out, g, h[0], h[1]);
    }
  }
}

static void hv_generic(const double* A, int n, int np, const double* s, const double* c2, int nq, double* Hs, double* sc)
{
  for (int j = 0; j < np; ++j) Hs[j] = 0.0;
  const int G = np / 4;
  for (int q = 0; q < nq; q += 4 * NB)
  {
    dots_nb(A, n, np, s, c2, q, sc);
    const double* B = A + (std::ptrdiff_t)(q / 4) * np * 4;
    for (int g = 0; g < G; ++g)
    {
      HV h{};  // AVX2: h[1] unused, but passed by value
      hload(Hs, g, h[0], h[1]);
      if (kPrefetch && q + (kPrefetch + 1) * 4 * NB <= nq)
        for (int b = 0; b < NB; ++b)
        {
          const char* pf = reinterpret_cast<const char*>(B + (std::ptrdiff_t)(kPrefetch * NB + b) * np * 4 + 16 * g);
          _mm_prefetch(pf, _MM_HINT_T0);
          _mm_prefetch(pf + 64, _MM_HINT_T0);
        }
      for (int b = 0; b < NB; ++b) scatter_blk(B + (std::ptrdiff_t)b * np * 4, g, sc + q + 4 * b, h[0], h[1]);
      hstore(Hs, g, h[0], h[1]);
    }
  }
}

#ifdef __AVX2__
// GN column groups (4 columns each) starting at gb accumulated in registers, rows in order
template <int GN> static void scatter_reg(const double* A, int np, const double* sc, int nq, int gb, double* out)
{
  __m256d h[GN];
  for (int g = 0; g < GN; ++g) h[g] = _mm256_setzero_pd();
  __m256d dummy;
  for (int q = 0; q < nq; q += 4)
  {
    const double* B = A + (std::ptrdiff_t)(q / 4) * np * 4;
    for (int g = 0; g < GN; ++g) scatter_blk(B, gb + g, sc + q, h[g], dummy);
  }
  for (int g = 0; g < GN; ++g) _mm256_storeu_pd(out + 4 * (gb + g), h[g]);
}
static void scatter(const double* A, int np, const double* sc, int nq, int gb, int ge, double* out)
{
  switch (ge - gb)
  {
    case 1: return scatter_reg<1>(A, np, sc, nq, gb, out);
    case 2: return scatter_reg<2>(A, np, sc, nq, gb, out);
    case 3: return scatter_reg<3>(A, np, sc, nq, gb, out);
    case 4: return scatter_reg<4>(A, np, sc, nq, gb, out);
    case 5: return scatter_reg<5>(A, np, sc, nq, gb, out);
    case 6: return scatter_reg<6>(A, np, sc, nq, gb, out);
    case 7: return scatter_reg<7>(A, np, sc, nq, gb, out);
    case 8: return scatter_reg<8>(A, np, sc, nq, gb, out);
    default: return scatter_generic(A, np, sc, nq, gb, ge, out);
  }
}
// fused Hv, G = np/4 column groups in registers; the dots of the next superblock are interleaved with the
// scatter of this one
template <int G> static void hv_reg(const double* A, int n, const double* s, const double* c2, int nq, double* Hs, double* sc)
{
  constexpr int np = 4 * G;
  constexpr std::ptrdiff_t bs = np * 4;
  __m256d h[G];
  for (int g = 0; g < G; ++g) h[g] = _mm256_setzero_pd();
  __m256d dummy;
  if (nq > 0) dots_nb(A, n, np, s, c2, 0, sc);  // no active rows: Hs = 0
  for (int q = 0; q < nq; q += 4 * NB)
  {
    const bool more = q + 4 * NB < nq;
    const bool pf = q + (kPrefetch + 1) * 4 * NB <= nq;
    const double* B = A + (std::ptrdiff_t)(q / 4) * bs;
    const double* Bn = B + NB * bs;
    __m256d acc[NB];
    for (int b = 0; b < NB; ++b) acc[b] = _mm256_setzero_pd();
    int j = 0;
    for (int b = 0; b < NB; ++b)
    {
      for (int g = 0; g < G; ++g)
      {
        if (kPrefetch && pf)
        {
          const char* p = reinterpret_cast<const char*>(B + kPrefetch * NB * bs + b * bs + 16 * g);
          _mm_prefetch(p, _MM_HINT_T0);
          _mm_prefetch(p + 64, _MM_HINT_T0);
        }
        scatter_blk(B + b * bs, g, sc + q + 4 * b, h[g], dummy);
        if (more && j < n)
        {
          const __m256d sj = _mm256_broadcast_sd(s + j);
          for (int bb = 0; bb < NB; ++bb) acc[bb] = _mm256_add_pd(acc[bb], _mm256_mul_pd(_mm256_load_pd(Bn + bb * bs + 4 * j), sj));
          ++j;
        }
      }
    }
    if (more)
    {
      for (; j < n; ++j)
      {
        const __m256d sj = _mm256_broadcast_sd(s + j);
        for (int bb = 0; bb < NB; ++bb) acc[bb] = _mm256_add_pd(acc[bb], _mm256_mul_pd(_mm256_load_pd(Bn + bb * bs + 4 * j), sj));
      }
      for (int b = 0; b < NB; ++b) _mm256_storeu_pd(sc + q + 4 * NB + 4 * b, _mm256_mul_pd(_mm256_loadu_pd(c2 + q + 4 * NB + 4 * b), acc[b]));
    }
  }
  for (int g = 0; g < G; ++g) _mm256_storeu_pd(Hs + 4 * g, h[g]);
}
static void hv(const double* A, int n, int np, const double* s, const double* c2, int nq, double* Hs, double* sc)
{
  switch (np / 4)
  {
    case 1: return hv_reg<1>(A, n, s, c2, nq, Hs, sc);
    case 2: return hv_reg<2>(A, n, s, c2, nq, Hs, sc);
    case 3: return hv_reg<3>(A, n, s, c2, nq, Hs, sc);
    case 4: return hv_reg<4>(A, n, s, c2, nq, Hs, sc);
    case 5: return hv_reg<5>(A, n, s, c2, nq, Hs, sc);
    case 6: return hv_reg<6>(A, n, s, c2, nq, Hs, sc);
    case 7: return hv_reg<7>(A, n, s, c2, nq, Hs, sc);
    case 8: return hv_reg<8>(A, n, s, c2, nq, Hs, sc);
    case 9: return hv_reg<9>(A, n, s, c2, nq, Hs, sc);
    case 10: return hv_reg<10>(A, n, s, c2, nq, Hs, sc);
    default: return hv_generic(A, n, np, s, c2, nq, Hs, sc);
  }
}
#else
static void scatter(const double* A, int np, const double* sc, int nq, int gb, int ge, double* out) { scatter_generic(A, np, sc, nq, gb, ge, out); }
static void hv(const double* A, int n, int np, const double* s, const double* c2, int nq, double* Hs, double* sc) { hv_generic(A, n, np, s, c2, nq, Hs, sc); }
#endif

}  // namespace TRON_KERNELS_NS

const Ops TRON_KERNELS_NS::ops = {
#ifdef __AVX2__
  "avx2",
#else
  "sse2",
#endif
  TRON_KERNELS_NS::margins_bt, TRON_KERNELS_NS::pack_bt, TRON_KERNELS_NS::fun_bt, TRON_KERNELS_NS::build_bt,
  TRON_KERNELS_NS::scatter, TRON_KERNELS_NS::hv};

}}}}  // namespace OpenMS::Internal::Percolator::TronKernels
