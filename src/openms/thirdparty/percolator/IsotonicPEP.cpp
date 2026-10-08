#include "Globals.h"
#include <vector>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <numeric>    // for std::iota
#include <cassert>
#include <iostream>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "MonotoneRegressor.h"
#include "IsotonicPEP.h"

namespace OpenMS { namespace Internal { namespace Percolator {

using clock_type = std::chrono::high_resolution_clock;

namespace {
std::vector<double> q_values_to_raw_pep(const std::vector<double>& q_values) {
  const size_t n = q_values.size();
  std::vector<double> raw_pep(n);
  if (n == 0) {
    return raw_pep;
  }

  double prev_qn = 0.0;
  for (size_t i = 0; i < n; ++i) {
    assert((i + 1 == n) || (q_values[i] <= q_values[i + 1]));
    const double qn = q_values[i] * static_cast<double>(i + 1);
    raw_pep[i] = qn - prev_qn;
    prev_qn = qn;
  }
  return raw_pep;
}
}  // namespace

InferPEP::InferPEP(bool use_ispline)
{
    if (use_ispline) {
        regressor_ptr_ = make_monotone_regressor(RegressorType::ISPLINE_TRR);
        if (VERB > 1) {
            std::cerr << "Performing isotonic regression using I-Splines" << std::endl;
        }                
    } else {
        regressor_ptr_ = make_monotone_regressor(RegressorType::PAVA);
        if (VERB > 1) {
            std::cerr << "Performing isotonic regression using PAVA" << std::endl;
        }
    }
}

std::vector<double> InferPEP::q_to_pep(const std::vector<double>& q_values) {
    qs = q_values;
    std::vector<double> raw_pep = q_values_to_raw_pep(q_values);
    // Bayesian pseudocount: spread half a false discovery uniformly across
    // all peptides.  This prevents exact-zero PEPs in the leading tail
    // (where q = 0) and gives the monotone smoother something to fit.
    if (!raw_pep.empty()) {
      const double pseudo = 0.5 / static_cast<double>(raw_pep.size());
      for (auto& v : raw_pep) v += pseudo;
    }
    const double eps = 1e-10;
    return regressor_ptr_->fit_y(raw_pep, /*clip_lo=*/eps);
}

std::vector<double> InferPEP::qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores) {
    if (q_values.size() != scores.size()) {
      throw std::invalid_argument("InferPEP::qns_to_pep: q_values and scores size mismatch");
    }
    qs = q_values;
    std::vector<double> raw_pep = q_values_to_raw_pep(q_values);
    if (!raw_pep.empty()) {
      const double pseudo = 0.5 / static_cast<double>(raw_pep.size());
      for (auto& v : raw_pep) v += pseudo;
    }
    const double eps = 1e-10;
    return regressor_ptr_->fit_xy(scores, raw_pep, /*clip_lo=*/eps);
}

std::vector<double>
InferPEP::tdc_to_pep(const std::vector<double>& is_decoy,
                     const std::vector<double>& scores) {
  if (!scores.empty() && scores.size() != is_decoy.size()) {
    throw std::invalid_argument("InferPEP::tdc_to_pep: scores and is_decoy size mismatch");
  }
  const bool print_timing = (VERB > 2);
  const auto start = print_timing ? clock_type::now() : clock_type::time_point{};
  if (VERB > 2) std::cerr << "[TIMING] entering tdc_to_pep\n";

  const double epsilon = 1e-20;

  // Use scores directly in fit_xy path (no synthetic anchor point).
  const bool will_use_fit_xy = !scores.empty();

  std::vector<double> is_dec;
  std::vector<double> sc;
  if (will_use_fit_xy) {
    assert(scores.size() == is_decoy.size());
    is_dec = is_decoy;
    sc = scores;

    // A zero-decoy anchor just above the best observed score lets the monotone
    // fit decay toward zero in the extreme high-score tail instead of forcing
    // a positive floor from sparse decoys deeper in the list.
    const auto mm = std::minmax_element(sc.begin(), sc.end());
    const double score_span = *mm.second - *mm.first;
    const double delta = std::max(score_span * 1e-6, 1e-12);
    sc.insert(sc.begin(), *mm.second + delta);
    is_dec.insert(is_dec.begin(), 0.0);
  } else {
    is_dec = is_decoy;
  }

  // Fit decoy rate p(decoy | x)
  const std::vector<double> decoy_rate = will_use_fit_xy
      ? regressor_ptr_->fit_xy(sc, is_dec, /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon)
      : regressor_ptr_->fit_y(is_dec,      /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon);

  const size_t offset = will_use_fit_xy ? 1u : 0u;
  std::vector<double> pep_iso(is_decoy.size());
  for (size_t i = 0; i < pep_iso.size(); ++i) {
    double p = decoy_rate[i + offset];
    double pep = p / (1.0 - p);  
    pep = std::max(0.0, std::min(1.0, pep));
    pep_iso[i] = pep;
  }

  if (print_timing) {
    const double duration = std::chrono::duration_cast<std::chrono::duration<double>>(clock_type::now() - start).count();
    std::cerr << "[TIMING] tdc_to_pep duration: " << duration << " seconds\n";
  }

  return pep_iso;
}

}}}  // namespace OpenMS::Internal::Percolator
