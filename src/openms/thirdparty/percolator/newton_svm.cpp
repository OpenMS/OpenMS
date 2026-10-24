/* See newton_svm.h. */
#include "newton_svm.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <vector>

#include <Eigen/Dense>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

namespace {

// rows per block: the unit of the partial sums (added in block order) and of the tasks
constexpr int kBlockRows = 8192;
// problems of at least this many blocks run their blocks as tasks within a parallel region
constexpr int kMinTaskBlocks = 8;

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::VectorXd;

// f(b) for the blocks b = 0 .. numBlocks - 1; as tasks of the team within a parallel region
template <class F>
void forEachBlock(int numBlocks, const F& f) {
#ifdef _OPENMP
  if (numBlocks >= kMinTaskBlocks && omp_in_parallel() && omp_get_num_threads() > 1) {
    // a few tasks per thread (libgomp runs a taskloop undeferred if the team has too many tasks)
    const int numTasks = std::min(numBlocks, 2 * omp_get_num_threads());
#pragma omp taskloop num_tasks(numTasks) default(shared)
    for (int b = 0; b < numBlocks; ++b) f(b);
    return;
  }
#endif
  for (int b = 0; b < numBlocks; ++b) f(b);
}

struct Problem {
  const double* const* rows;  // m rows of n - 1 features
  const double* Y;            // labels (+-1)
  int m;
  int n;          // columns: features and the bias
  double c2pos;   // 2 C_i of the rows with Y > 0: cpos / lambda
  double c2neg;   // and of the others
  int numBlocks() const { return (m + kBlockRows - 1) / kBlockRows; }
  int blockBegin(int b) const { return b * kBlockRows; }
  int blockEnd(int b) const { return std::min(m, (b + 1) * kBlockRows); }
  double c2(int i) const { return Y[i] > 0.0 ? c2pos : c2neg; }
};

// x = (row, 1)
inline void loadRow(const Problem& p, int i, double* x) {
  std::memcpy(x, p.rows[i], sizeof(double) * static_cast<std::size_t>(p.n - 1));
  x[p.n - 1] = 1.0;
}

// H[j][k] += sum_r a_r[j] x_r[k] for k >= j (upper triangle of an n x n row-major matrix), r = 0..3
inline void rank4(double* H, int n, const double* x, const double* a) {
  const double* x0 = x;
  const double* x1 = x + n;
  const double* x2 = x + 2 * n;
  const double* x3 = x + 3 * n;
  for (int j = 0; j < n; ++j) {
    const double a0 = a[j], a1 = a[n + j], a2 = a[2 * n + j], a3 = a[3 * n + j];
    double* h = H + static_cast<std::size_t>(j) * n;
    for (int k = j; k < n; ++k) h[k] += a0 * x0[k] + a1 * x1[k] + a2 * x2[k] + a3 * x3[k];
  }
}

// collects rows (x_r, a_r = c_r x_r) in groups of four for rank4
struct Rank4Batch {
  double x[4 * kNewtonMaxColumns];
  double a[4 * kNewtonMaxColumns];
  int r = 0;
  void add(const Problem& p, int i, double c, double* H) {
    double* xr = x + r * p.n;
    double* ar = a + r * p.n;
    loadRow(p, i, xr);
    for (int j = 0; j < p.n; ++j) ar[j] = c * xr[j];
    if (++r == 4) {
      rank4(H, p.n, x, a);
      r = 0;
    }
  }
  void flush(const Problem& p, double* H) {
    if (r == 0) return;
    for (; r < 4; ++r) {
      std::fill(x + r * p.n, x + (r + 1) * p.n, 0.0);
      std::fill(a + r * p.n, a + (r + 1) * p.n, 0.0);
    }
    rank4(H, p.n, x, a);
    r = 0;
  }
};

// out[i] = x_i . w for the rows i0 .. i1 - 1, four rows at a time (each sum in column order)
inline void dots(const Problem& p, int i0, int i1, const double* w, double* out) {
  const int nf = p.n - 1;
  const double bias = w[nf];
  int i = i0;
  for (; i + 4 <= i1; i += 4) {
    const double* r0 = p.rows[i];
    const double* r1 = p.rows[i + 1];
    const double* r2 = p.rows[i + 2];
    const double* r3 = p.rows[i + 3];
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    for (int j = 0; j < nf; ++j) {
      s0 += r0[j] * w[j];
      s1 += r1[j] * w[j];
      s2 += r2[j] * w[j];
      s3 += r3[j] * w[j];
    }
    out[i] = s0 + bias;
    out[i + 1] = s1 + bias;
    out[i + 2] = s2 + bias;
    out[i + 3] = s3 + bias;
  }
  for (; i < i1; ++i) {
    const double* r0 = p.rows[i];
    double s0 = 0.0;
    for (int j = 0; j < nf; ++j) s0 += r0[j] * w[j];
    out[i] = s0 + bias;
  }
}

// a breakpoint of the line search: at d, phi'(delta) = L + delta R changes by dL + delta dR
struct Breakpoint {
  double d;
  double dL;
  double dR;
};

// the breakpoints are counted in buckets of quarter octaves of delta (bucket 0: [0, 2^-16), the last
// one: [2^16, inf)), so that the line search needs to select only among those of one bucket
constexpr int kBucketOctaves = 16;
constexpr int kBuckets = 2 + 2 * kBucketOctaves * 4;
constexpr double kQuarter[4] = {0.5, 0.59460355750136053, 0.70710678118654752, 0.84089641525371454};

inline int bucketOf(double d) {
  if (!(d >= 0x1p-16)) return 0;
  if (d >= 0x1p16) return kBuckets - 1;
  int e;
  const double f = std::frexp(d, &e);  // d = f 2^e, f in [0.5, 1), e - 1 in [-16, 15]
  const int q = (f >= kQuarter[1]) + (f >= kQuarter[2]) + (f >= kQuarter[3]);
  return 1 + (e - 1 + kBucketOctaves) * 4 + q;
}

// the lower end of bucket k >= 1 (bucketOf(d) >= k exactly if d >= bucketBegin(k))
inline double bucketBegin(int k) {
  if (k == kBuckets - 1) return 0x1p16;
  const int e = (k - 1) / 4 - kBucketOctaves + 1;
  return std::ldexp(kQuarter[(k - 1) % 4], e);
}

struct BucketSum {
  double dL = 0.0;
  double dR = 0.0;
  int count = 0;
};

// the smallest delta >= 0 with L(delta) + delta R(delta) >= 0, where L and R take the changes of
// the breakpoints below delta; crossed = the number of breakpoints below it. Selection by
// partitioning (expected linear time), the last few sorted.
double lineSearch(std::vector<Breakpoint>& bp, double L, double R, std::size_t& crossed) {
  crossed = 0;
  auto first = bp.begin();
  auto last = bp.end();
  while (last - first > 32) {
    const double a = first->d, b = first[(last - first) / 2].d, c = (last - 1)->d;
    const double pivot = std::max(std::min(a, b), std::min(std::max(a, b), c));  // median of three
    auto mid = std::partition(first, last, [pivot](const Breakpoint& q) { return q.d <= pivot; });
    double L2 = L, R2 = R;
    for (auto q = first; q != mid; ++q) {
      L2 += q->dL;
      R2 += q->dR;
    }
    if (L2 + pivot * R2 >= 0.0) {
      // the root is at most pivot: among the breakpoints below it (phi' is continuous, those at
      // the pivot do not change it there)
      last = std::partition(first, mid, [pivot](const Breakpoint& q) { return q.d < pivot; });
    } else {
      crossed += static_cast<std::size_t>(mid - first);
      L = L2;
      R = R2;
      first = mid;
    }
  }
  std::sort(first, last, [](const Breakpoint& x, const Breakpoint& y) { return x.d < y.d; });
  for (auto q = first; q != last; ++q) {
    if (L + q->d * R >= 0.0) break;
    L += q->dL;
    R += q->dR;
    ++crossed;
  }
  return R > 0.0 ? std::max(0.0, -L / R) : 1.0;
}

}  // namespace

int newtonSvm(const AlgIn& set, const options& opts, double cpos, double cneg, double eps,
              double* w, double* outputs, NewtonStats* stats) {
  const int n = set.n;
  const int m = set.m;
  if (n < 1 || n > kNewtonMaxColumns || m < 1 || !(opts.lambda > 0.0)) return -1;
  Problem p{set.vals, set.Y, m, n, cpos / opts.lambda, cneg / opts.lambda};
  const int nb = p.numBlocks();
  const std::size_t nn = static_cast<std::size_t>(n) * n;

  std::vector<double> o(m), obar(m);              // X w and X wbar
  std::vector<double> partH(nb * nn), partB(static_cast<std::size_t>(nb) * n);
  std::vector<double> partL(nb), partR(nb);
  std::vector<long long> partActive(nb);
  std::vector<int> partCount(nb);
  std::vector<BucketSum> partBuckets(static_cast<std::size_t>(nb) * kBuckets);
  std::vector<std::size_t> partOffset(nb);
  std::vector<Breakpoint> bpBlocks(m), bp;

  Vector wv = Eigen::Map<const Vector>(w, n);
  const Vector wStart = wv;
  NewtonStats st;
  double gnorm0 = 0.0;
  const bool cold = wv.isZero(0.0);
  if (cold) {
    std::fill(o.begin(), o.end(), 0.0);  // all rows active: the first Hessian pass gives g(0) = -b
  } else {
    // X w, and g(0) = -sum_i 2 C_i t_i x_i
    forEachBlock(nb, [&](int b) {
      const int i0 = p.blockBegin(b), i1 = p.blockEnd(b);
      dots(p, i0, i1, w, o.data());
      double* pb = partB.data() + static_cast<std::size_t>(b) * n;
      std::fill(pb, pb + n, 0.0);
      double x[kNewtonMaxColumns];
      for (int i = i0; i < i1; ++i) {
        loadRow(p, i, x);
        const double ct = p.c2(i) * set.Y[i];
        for (int j = 0; j < n; ++j) pb[j] += ct * x[j];
      }
    });
    Vector g0 = Vector::Zero(n);
    for (int b = 0; b < nb; ++b) g0 += Eigen::Map<const Vector>(partB.data() + static_cast<std::size_t>(b) * n, n);
    gnorm0 = g0.norm();
    ++st.passes;
  }

  Matrix H(n, n);
  Vector rhs(n), wbar(n), pdir(n);
  bool haveH = false;
  int ret = 0;
  for (int iter = 0;; ++iter) {
    if (!haveH) {
      // H = I + sum_{t_i o_i < 1} 2 C_i x_i x_i^T, rhs = sum 2 C_i t_i x_i
      forEachBlock(nb, [&](int b) {
        double* ph = partH.data() + b * nn;
        double* pb = partB.data() + static_cast<std::size_t>(b) * n;
        std::fill(ph, ph + nn, 0.0);
        std::fill(pb, pb + n, 0.0);
        Rank4Batch batch;
        long long active = 0;
        for (int i = p.blockBegin(b), i1 = p.blockEnd(b); i < i1; ++i) {
          const double t = set.Y[i];
          if (t * o[i] < 1.0) {
            const double c = p.c2(i);
            batch.add(p, i, c, ph);
            const double* x = batch.x + ((batch.r + 3) & 3) * n;  // the row just added
            const double ct = c * t;
            for (int j = 0; j < n; ++j) pb[j] += ct * x[j];
            ++active;
          }
        }
        batch.flush(p, ph);
        partActive[b] = active;
      });
      H.setZero();
      rhs.setZero();
      for (int b = 0; b < nb; ++b) {
        const double* ph = partH.data() + b * nn;
        for (int j = 0; j < n; ++j) {
          for (int k = j; k < n; ++k) H(j, k) += ph[static_cast<std::size_t>(j) * n + k];
        }
        rhs += Eigen::Map<const Vector>(partB.data() + static_cast<std::size_t>(b) * n, n);
        st.hessianRows += partActive[b];
      }
      for (int j = 0; j < n; ++j) {
        H(j, j) += 1.0;
        for (int k = j + 1; k < n; ++k) H(k, j) = H(j, k);
      }
      haveH = true;
      if (cold && iter == 0) gnorm0 = rhs.norm();
    }
    // g(w) = w + sum_I 2 C_i x_i (x_i . w - t_i) = H w - rhs
    const double gnorm = (H * wv - rhs).norm();
    st.gradRatio = gnorm0 > 0.0 ? gnorm / gnorm0 : 0.0;
    if (gnorm <= eps * gnorm0) {
      ret = 1;
      break;
    }
    if (iter >= opts.mfnitermax) break;
    if (!std::isfinite(gnorm)) {
      ret = -1;
      break;
    }

    // the Newton point: the minimum of the regularized least squares problem on the active rows
    Eigen::LLT<Matrix> llt(H);
    if (llt.info() != Eigen::Success) {
      ret = -1;
      break;
    }
    wbar = llt.solve(rhs);
    pdir = wbar - wv;
    ++st.iterations;

    // X wbar, and the breakpoints of phi(delta) = F(w + delta (wbar - w))
    forEachBlock(nb, [&](int b) {
      const int i0 = p.blockBegin(b), i1 = p.blockEnd(b);
      dots(p, i0, i1, wbar.data(), obar.data());
      double sL = 0.0, sR = 0.0;
      Breakpoint* q = bpBlocks.data() + i0;
      BucketSum* buckets = partBuckets.data() + static_cast<std::size_t>(b) * kBuckets;
      std::fill(buckets, buckets + kBuckets, BucketSum());
      int count = 0;
      for (int i = i0; i < i1; ++i) {
        const double t = set.Y[i];
        const double a = 1.0 - t * o[i];
        const double s = t * (obar[i] - o[i]);
        const double c = p.c2(i);
        Breakpoint br;
        if (a > 0.0) {
          sL += c * a * s;
          sR += c * s * s;
          if (!(s > 0.0)) continue;
          br = Breakpoint{a / s, c * a * s, -c * s * s};  // leaves
        } else if (s < 0.0) {
          br = Breakpoint{a / s, -c * a * s, c * s * s};  // enters
        } else {
          continue;
        }
        q[count++] = br;
        BucketSum& bs = buckets[bucketOf(br.d)];
        bs.dL += br.dL;
        bs.dR += br.dR;
        ++bs.count;
      }
      partL[b] = sL;
      partR[b] = sR;
      partCount[b] = count;
    });
    ++st.passes;
    double L = wv.dot(pdir), R = pdir.squaredNorm();
    for (int b = 0; b < nb; ++b) {
      L -= partL[b];
      R += partR[b];
    }
    // the bucket of the root: phi'(delta) at the lower ends of the buckets, with the changes of the
    // buckets below (the sums of a bucket in block order)
    std::size_t crossed = 0;
    int kRoot = kBuckets - 1;
    for (int k = 0; k < kBuckets - 1; ++k) {
      BucketSum tot;
      for (int b = 0; b < nb; ++b) {
        const BucketSum& bs = partBuckets[static_cast<std::size_t>(b) * kBuckets + k];
        tot.dL += bs.dL;
        tot.dR += bs.dR;
        tot.count += bs.count;
      }
      const double next = bucketBegin(k + 1);
      if (L + tot.dL + next * (R + tot.dR) >= 0.0) {
        kRoot = k;
        break;
      }
      L += tot.dL;
      R += tot.dR;
      crossed += static_cast<std::size_t>(tot.count);
    }
    // its breakpoints, in block order
    std::size_t numRoot = 0;
    for (int b = 0; b < nb; ++b) {
      partOffset[b] = numRoot;
      numRoot += static_cast<std::size_t>(partBuckets[static_cast<std::size_t>(b) * kBuckets + kRoot].count);
    }
    bp.resize(numRoot);
    forEachBlock(nb, [&](int b) {
      const Breakpoint* q = bpBlocks.data() + p.blockBegin(b);
      const Breakpoint* qEnd = q + partCount[b];
      Breakpoint* out = bp.data() + partOffset[b];
      for (; q != qEnd; ++q) {
        if (bucketOf(q->d) == kRoot) *out++ = *q;
      }
    });
    std::size_t crossedInBucket = 0;
    const double delta = lineSearch(bp, L, R, crossedInBucket);
    crossed += crossedInBucket;
    if (!std::isfinite(delta)) {
      ret = -1;
      break;
    }
    wv += delta * pdir;
    forEachBlock(nb, [&](int b) {
      for (int i = p.blockBegin(b), i1 = p.blockEnd(b); i < i1; ++i) o[i] += delta * (obar[i] - o[i]);
    });
    // without a breakpoint below delta, the active rows (and so H and rhs) are the same
    haveH = crossed == 0;
  }
  if (ret < 0) {
    Eigen::Map<Vector>(w, n) = wStart;
    return -1;
  }
  Eigen::Map<Vector>(w, n) = wv;
  if (outputs) std::copy(o.begin(), o.end(), outputs);
  if (stats) *stats = st;
  return ret;
}

bool ldaSvm(const AlgIn& set, const options& opts, double cpos, double cneg, double* w) {
  const int n = set.n;
  const int m = set.m;
  if (n < 2 || n > kNewtonMaxColumns || m < 1 || !(opts.lambda > 0.0)) return false;
  Problem p{set.vals, set.Y, m, n, cpos / opts.lambda, cneg / opts.lambda};
  const int nb = p.numBlocks();
  const std::size_t nn = static_cast<std::size_t>(n) * n;
  // per block and class: sum x x^T of x = (row, 1), i.e. also the sums of the features and the count
  std::vector<double> part(2 * nb * nn);
  forEachBlock(nb, [&](int b) {
    double* hp = part.data() + 2 * b * nn;
    double* hn = hp + nn;
    std::fill(hp, hp + 2 * nn, 0.0);
    Rank4Batch pos, neg;
    for (int i = p.blockBegin(b), i1 = p.blockEnd(b); i < i1; ++i) {
      if (set.Y[i] > 0.0) pos.add(p, i, 1.0, hp);
      else neg.add(p, i, 1.0, hn);
    }
    pos.flush(p, hp);
    neg.flush(p, hn);
  });
  Matrix S[2] = {Matrix::Zero(n, n), Matrix::Zero(n, n)};
  for (int b = 0; b < nb; ++b) {
    for (int c = 0; c < 2; ++c) {
      const double* h = part.data() + (2 * b + c) * nn;
      for (int j = 0; j < n; ++j) {
        for (int k = j; k < n; ++k) S[c](j, k) += h[static_cast<std::size_t>(j) * n + k];
      }
    }
  }
  const int nf = n - 1;
  for (int c = 0; c < 2; ++c) {
    for (int j = 0; j < n; ++j) {
      for (int k = j + 1; k < n; ++k) S[c](k, j) = S[c](j, k);
    }
  }
  const double cnt[2] = {S[0](nf, nf), S[1](nf, nf)};
  if (cnt[0] < 2.0 || cnt[1] < 2.0) return false;
  // pooled within-class scatter and the class means of the features
  Vector mu[2];
  Matrix Sw = Matrix::Zero(nf, nf);
  for (int c = 0; c < 2; ++c) {
    mu[c] = S[c].col(nf).head(nf) / cnt[c];
    Sw += S[c].topLeftCorner(nf, nf) - cnt[c] * mu[c] * mu[c].transpose();
  }
  Sw /= (cnt[0] + cnt[1] - 2.0);
  const Vector diff = mu[0] - mu[1];
  // a small ridge for constant or collinear features
  double ridge = 1e-10 * std::max(Sw.trace() / nf, 1e-300);
  Vector wt;
  for (int attempt = 0;; ++attempt) {
    Eigen::LLT<Matrix> llt(Sw + ridge * Matrix::Identity(nf, nf));
    if (llt.info() == Eigen::Success) {
      wt = llt.solve(diff);
      if (wt.allFinite()) break;
    }
    if (attempt == 6) return false;
    ridge *= 1e3;
  }
  // gamma and b0 minimizing F((gamma wt, b0)) with all rows active (Keerthi & DeCoste, eq. 8/9):
  // [|wt|^2 + sum c s^2, sum c s; sum c s, 1 + sum c] (gamma, b0) = (sum c t s, sum c t), s_i = wt . x_i,
  // c_i = 2 C_i, from the sums of the classes
  const double cc[2] = {p.c2pos, p.c2neg};
  const double tt[2] = {1.0, -1.0};
  double css = 0.0, cs = 0.0, cts = 0.0, csum = 0.0, ctsum = 0.0;
  for (int c = 0; c < 2; ++c) {
    const double sumS = wt.dot(S[c].col(nf).head(nf));
    css += cc[c] * wt.dot(S[c].topLeftCorner(nf, nf) * wt);
    cs += cc[c] * sumS;
    cts += cc[c] * tt[c] * sumS;
    csum += cc[c] * cnt[c];
    ctsum += cc[c] * tt[c] * cnt[c];
  }
  const double p11 = wt.squaredNorm() + css, p12 = cs, p22 = 1.0 + csum;
  const double det = p11 * p22 - p12 * p12;
  if (!(det > 0.0)) return false;
  const double gamma = (p22 * cts - p12 * ctsum) / det;
  const double b0 = (p11 * ctsum - p12 * cts) / det;
  if (!std::isfinite(gamma) || !std::isfinite(b0)) return false;
  for (int j = 0; j < nf; ++j) w[j] = gamma * wt(j);
  w[nf] = b0;
  return true;
}

double svmObjective(const AlgIn& set, double lambda, double cpos, double cneg, const double* w) {
  const int n = set.n;
  double f = 0.0;
  for (int j = 0; j < n; ++j) f += w[j] * w[j];
  f *= 0.5;
  const double inv2l = 1.0 / (2.0 * lambda);
  for (int i = 0; i < set.m; ++i) {
    double s = 0.0;
    for (int j = 0; j < n - 1; ++j) s += set.vals[i][j] * w[j];
    s += w[n - 1];
    const double margin = set.Y[i] * s;
    if (margin < 1.0) {
      const double d = 1.0 - margin;
      f += (set.Y[i] > 0.0 ? cpos : cneg) * inv2l * d * d;
    }
  }
  return f;
}

}}}  // namespace OpenMS::Internal::Percolator
