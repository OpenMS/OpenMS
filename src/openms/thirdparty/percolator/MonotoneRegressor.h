#pragma once
#include <vector>
#include <memory>


namespace OpenMS { namespace Internal { namespace Percolator {

// Constants for tuning
constexpr double DEFAULT_LAMBDA = 1e-4;
constexpr double DEFAULT_SMOOTH_LAMBDA = 1e-3;
constexpr int DEFAULT_NUM_KNOTS = 50;
constexpr double DEFAULT_SKEW_FACTOR = 0.75;


struct MonotoneParams {
  // Shared
  double clip_lo = 0.0;
  double clip_hi = 1.0;
  // I-spline TRR specific
  double ridge_lambda = DEFAULT_LAMBDA;
  double smooth_lambda = DEFAULT_SMOOTH_LAMBDA;
  bool y_decreasing_in_x = true;
  int ispline_degree = 3;
  bool include_intercept = true;
  int intercept_col = 0; // if include_intercept, usually 0
  std::vector<double> knots; // knot policy fills this
};

class MonotoneRegressor {
public:
  virtual ~MonotoneRegressor() = default;

  // Fit y ~ f(x) with optional weights and output range clamp
  virtual std::vector<double> fit_xy(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     double clip_lo = 0.0,
                                     double clip_hi = 1.0) = 0;

  // Fit y ~ f(index) (no x); useful for the tdc path when no scores
  virtual std::vector<double> fit_y(const std::vector<double>& y,
                                    double clip_lo = 0.0,
                                    double clip_hi = 1.0) = 0;

  virtual void set_params(const MonotoneParams& p) { params_ = p; }
  virtual MonotoneParams params() const { return params_; }

protected:
  MonotoneParams params_;
};

enum class RegressorType {
  PAVA,          // isotonic/PAVA implementation
  ISPLINE_TRR    // trust-region reflective constrained ridge I-spline
};

std::unique_ptr<MonotoneRegressor> make_monotone_regressor(RegressorType type);

}}}  // namespace OpenMS::Internal::Percolator
