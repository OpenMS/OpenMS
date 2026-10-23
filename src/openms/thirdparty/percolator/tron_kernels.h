/*
 * OpenMS addition (not part of upstream Percolator or percolator-tron): kernels of the
 * L2-loss SVC function object of the TRON solver (ssl.cpp), bit-identical to the reference
 * TronFunction (Eigen row-major gemv in the SSE2 baseline, serial loops otherwise).
 *
 * Compiled twice, without Eigen: tron_kernels_sse2.cpp (baseline flags) and
 * tron_kernels_avx2.cpp (-mavx2, no FMA, no contraction); ssl.cpp picks one table at run time.
 *
 * Layout "BT4": rows in blocks of 4; block b holds np x 4 doubles, element (q, j) at
 *   A[(q/4)*np*4 + j*4 + q%4]; np = n rounded up to a multiple of 4; padded columns and rows are 0.
 * Every value is computed in the order of the reference:
 *   margin_i = (0.0 + ((e_i + o_i) [+ x_i,n-1 * w_n-1 if n odd])) * t_i, e/o = serial sums over the
 *              even/odd columns < 2*(n/2) (Eigen's Packet2d accumulator, then predux and the scalar tail)
 *   dot_q    = serial sum_{j<n} (x_qj * s_j) from 0.0;  sc_q = c2_q * dot_q  (c2 = 2.0*C)
 *   out_j    = serial sum over q in order of (sc_q * x_qj) from 0.0
 * Adding +0 or -0 (padded rows and columns) never changes these sums: an accumulator that starts
 * at +0.0 can never become -0 (round to nearest, no flush to zero). Padded rows do change Hv when
 * s holds an Inf or NaN (0 * Inf = NaN): FastTronFunction::Hv then loops over the active rows itself.
 *
 * The contract: the expressions of TronFunction in ssl.cpp, as compiled for the x86-64 baseline
 * (SSE2, no FMA; Eigen 5.0.1-dev), which these kernels reproduce:
 *   fun:  margin_ = (X_ * w).cwiseProduct(t_);           // Eigen row-major gemv, then * t_i
 *         f = 0.5 * w.squaredNorm();
 *         for i < m: if (margin_(i) < 1.0) { d = 1.0 - margin_(i); f += C_(i) * d * d; }
 *   grad: g = 0; for active i in order: z = 2.0 * C_(i) * t_(i) * (margin_(i) - 1.0);
 *                                       for j < n: g[j] += z * row_i[j];
 *         for j < n: g[j] += w[j];
 *   Hv:   Hs = 0; for active i in order: xTs = 0.0; for j < n: xTs += row_i[j] * s[j];
 *                                       scale = 2.0 * C_(i) * xTs; for j < n: Hs[j] += scale * row_i[j];
 *         for j < n: Hs[j] += s[j];
 *   tron(): Outputs = X * w (Eigen gemv). A change to any of these must be made in the kernels too;
 *   PercolatorTron_test and the self-check of tron() (at the first use of a kernel table) compare
 *   the two.
 */
#ifndef PERCOLATOR_TRON_KERNELS_H
#define PERCOLATOR_TRON_KERNELS_H

namespace OpenMS { namespace Internal { namespace Percolator { namespace TronKernels {

struct Ops
{
  const char* name;
  // X stored in BT4 (rows padded with zeros to a multiple of 32): margins of rows [b, e), multiples of 4
  void (*margins_bt)(const double* XB, int n, int np, const double* w, const double* t, int b, int e, double* margin);
  // rows idx[q], q in [qb, qe) (multiples of 4), of BT4-stored X into BT4 A; q >= na are zero rows
  void (*pack_bt)(const double* XB, int np, const int* idx, int na, int qb, int qe, double* A);
  // margins of all m rows + active rows in order (idx) + f += C_i * d * d serially in row order (*f in/out)
  // + (A != nullptr, capacity capA rows) the active rows packed into A, zero rows up to a multiple of 32;
  // *packed = whether A was filled. Returns the number of active rows.
  int (*fun_bt)(const double* XB, int n, int np, const double* w, const double* t, const double* C, int m,
                double* margin, int* idx, double* A, int capA, double* f, int* packed);
  // BT4 X from m row pointers (nf features each), bias column nf = 1.0, zero columns up to np, zero rows up
  // to a multiple of 32
  void (*build_bt)(const double* const* rows, int m, int nf, int np, double* XB);
  // out[4g .. 4g+3] for groups g in [gb, ge): sum over rows q < nq (multiple of 32) in order of sc[q] * A(q, j)
  void (*scatter)(const double* A, int np, const double* sc, int nq, int gb, int ge, double* out);
  // Hs[0..np) = scatter(sc), sc[q] = c2[q] * dot(A row q, s); sc is scratch of nq (multiple of 32) values
  void (*hv)(const double* A, int n, int np, const double* s, const double* c2, int nq, double* Hs, double* sc);
};

namespace sse2 { extern const Ops ops; }
namespace avx2 { extern const Ops ops; }

}}}}  // namespace OpenMS::Internal::Percolator::TronKernels

#endif
