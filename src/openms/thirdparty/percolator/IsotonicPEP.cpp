#include "Globals.h"
#include <vector>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <numeric>    // for std::iota
#include <cassert>
#include <iostream>

#ifdef _MSC_VER
  #define NOMINMAX
  #include <float.h>
  #define _finite(x) _finite(x)
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse> 
#include <chrono>
#include "IsotonicPEP.h"


/**
 * Perform standard PAVA (pool-adjacent-violators) to enforce a non-decreasing sequence
 * on the input 'values'.  
 * Here we extended the method to also merge block if average is outside 
 * the range [min_value, max_value] 
 *
 * Each final "block" is repeated block.count times in the output.
 *
 * Typically used in 1D for y-values only.
 */
std::vector<double> PavaRegression::pavaNonDecreasingRanged(
    const std::vector<double>& values,
    const double min_value,
    const double max_value
) const
{
    const int n = static_cast<int>(values.size());
    if (n == 0) {
        return {};
    }

    // We'll store blocks on a stack
    struct MergeBlock {
        double sum;
        int    count;
        double avg;
    };

    std::vector<MergeBlock> stack;
    stack.reserve(n);

    // 1. Left to right
    for (int i = 0; i < n; ++i) {
        MergeBlock newBlock { values[i], 1, values[i] };
        stack.push_back(newBlock);

        // 2. Merge while there's a violation of non-decreasing
        while (stack.size() > 1) {
            auto &top    = stack.back();
            auto &secTop = stack[stack.size() - 2];
            double mergedSum   = secTop.sum + top.sum;
            int    mergedCount = secTop.count + top.count;
            double mergedAvg   = mergedSum / mergedCount;

            // if (( secTop.avg > top.avg) || ( mergedAvg < min_value ) || ( mergedAvg > max_value )) {
            if ( secTop.avg > top.avg) {
                stack.pop_back();
                stack.pop_back();
                stack.push_back({ mergedSum, mergedCount, mergedAvg });
            } else {
                break;
            }
        }
    }

    // 3. Expand final solution
    std::vector<double> result;
    result.reserve(n); // approximate; actual size = sum(counts)

    for (auto &b : stack) {
        double val = std::max( std::min(b.avg, max_value), min_value);
        for (int c = 0; c < b.count; ++c) {
            result.push_back(val);
        }
    }
    return result;
}

double IsplineRegression::cubic_ispline(double x, double left, double right) const {
    if (x < left) return 0.0;
    if (x >= right) return 1.0;
    double u = (x - left) / (right - left);
    return 3 * u * u - 2 * u * u * u;
}

namespace util {
    template <typename T>
    const T& clamp(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

std::vector<double> IsplineRegression::fit_y(const std::vector<double>& y, double min_val, double max_val) const {
    std::vector<double> x(y.size());
    std::iota(x.begin(), x.end(), 0.0);
    return fit_xy(x, y, min_val, max_val);
}

IsplineRegression::BinnedData IsplineRegression::bin_data(
    const std::vector<double>& x, const std::vector<double>& y, int max_bins) const
{
    size_t n = x.size();
    if (n == 0 || max_bins <= 0) return {{}, {}, {}};

    std::vector<double> x_binned, y_binned, weights;

    double target_bin_size = static_cast<double>(n) / max_bins;
    double next_bin_threshold = target_bin_size;

    size_t bin_start = 0;
    for (size_t i = 0; i < n; ++i) {
        double bin_progress = static_cast<double>(i + 1);  // since i is inclusive
        if (bin_progress >= next_bin_threshold || i == n - 1) {
            size_t bin_end = i + 1;  // exclusive
            size_t bin_size = bin_end - bin_start;

            double x_sum = std::accumulate(x.begin() + bin_start, x.begin() + bin_end, 0.0);
            double y_sum = std::accumulate(y.begin() + bin_start, y.begin() + bin_end, 0.0);

            double x_avg = x_sum / bin_size;
            double y_avg = y_sum / bin_size;

            if (VERB > 4) {
                std::cerr << "Creating bin from " << bin_start << " to " << i
                          << ", x_avg: " << x_avg << ", y_avg: " << y_avg << "\n";
            }

            x_binned.push_back(x_avg);
            y_binned.push_back(y_avg);
            weights.push_back(static_cast<double>(bin_size));

            bin_start = bin_end;
            next_bin_threshold += target_bin_size;
        }
    }

    return {x_binned, y_binned, weights};
}

std::vector<double> IsplineRegression::compute_adaptive_knots(const std::vector<double>& x, const std::vector<double>& /* y */, int num_knots) const {
    
    std::vector<double> knots;
    knots.push_back(x.front());
    for (int i = 1; i < num_knots; ++i) {
        double q = 1.0 - std::pow(1.0 - static_cast<double>(i) / num_knots, skew_factor_);
        size_t idx = static_cast<size_t>(q * (x.size() - 1));
        knots.push_back(x[idx]);
    }
    knots.push_back(x.back());
    if (VERB > 3) {
        std::cerr << "Knots: ";
        for (const auto& k : knots) std::cerr << k << " ";
        std::cerr << "\n";
    }
    return knots;
}

Eigen::VectorXd IsplineRegression::fit_spline(const BinnedData& data, const std::vector<double>& knots, double lambda) const {
    int n = data.x.size(), k = knots.size();
    Eigen::MatrixXd X(n, k);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < k - 1; ++j)
            X(i, j) = cubic_ispline(data.x[i], knots[j], knots[j + 1]);
    X.col(k - 1).setOnes();

    Eigen::VectorXd y = Eigen::Map<const Eigen::VectorXd>(data.y.data(), n);
    Eigen::VectorXd w = Eigen::Map<const Eigen::VectorXd>(data.weights.data(), n);
    Eigen::VectorXd coeffs = (X.transpose() * w.asDiagonal() * X + lambda * Eigen::MatrixXd::Identity(k, k)).ldlt().solve(X.transpose() * w.asDiagonal() * y);

    if (VERB > 3) std::cerr << "Spline coefficients: " << coeffs.transpose() << "\n";

    return coeffs;
}

std::vector<double> IsplineRegression::fit_xy(const std::vector<double>& x, const std::vector<double>& y, double min_val, double max_val) const {
    if (x.empty() || y.empty()) {
        if (VERB > 0) std::cerr << "[WARNING] fit_xy called with empty input.\n";
        return {};
    }
    assert(x.size() == y.size());
    auto data = bin_data(x, y, num_bins_);
    // auto data = adaptive_bin(x, y, num_bins);
    auto knots = compute_adaptive_knots(data.x, data.y, std::min(50, (int)std::sqrt(data.x.size())));

    Eigen::VectorXd coeffs = fit_spline(data, knots, lambda_);

    std::vector<double> result(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        double pred = coeffs.tail(1)(0);
        for (size_t j = 0; j < knots.size() - 1; ++j)
            pred += coeffs(j) * cubic_ispline(x[i], knots[j], knots[j + 1]);
        result[i] = util::clamp(pred, min_val, max_val);
    }

    return PavaRegression().fit_xy(x, result, min_val, max_val);
}

std::vector<double> InferPEP::q_to_pep(const std::vector<double>& q_values) {
    qs = q_values;

    std::vector<double> qn(q_values.size());
    std::vector<int> indices(q_values.size());
    std::iota(indices.begin(), indices.end(), 1);

    for (size_t i = 0; i < q_values.size(); ++i) {
        qn[i] = q_values[i] * indices[i];
        assert((i == q_values.size() - 1) || (q_values[i] <= q_values[i + 1]));
    }
    std::vector<double> raw_pep(qn.size());
    raw_pep[0] = qn[0];
    for (size_t i = 1; i < qn.size(); ++i) {
        raw_pep[i] = qn[i] - qn[i - 1];
    }   
    return regressor_ptr_->fit_y(raw_pep, epsilon_, 1. - epsilon_);
}

std::vector<double> InferPEP::qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores) {
    qs = q_values;

    std::vector<double> qn(q_values.size());
    std::vector<int> indices(q_values.size());
    std::iota(indices.begin(), indices.end(), 1);

    for (size_t i = 0; i < q_values.size(); ++i) {
        qn[i] = q_values[i] * indices[i];
        assert((i == q_values.size() - 1) || (q_values[i] <= q_values[i + 1]));
    }

    std::vector<double> raw_pep(qn.size());
    raw_pep[0] = qn[0];
    for (size_t i = 1; i < qn.size(); ++i) {
        raw_pep[i] = qn[i] - qn[i - 1];
    }

    return regressor_ptr_->fit_xy(scores, raw_pep, epsilon_, 1. - epsilon_);
}

std::vector<double> InferPEP::tdc_to_pep(const std::vector<double>& is_decoy, const std::vector<double>& scores) {

    using namespace std::chrono;
    auto start = std::chrono::high_resolution_clock::now();
    if (VERB > 2)
        std::cerr << "[TIMING] entering tdc_to_pep\n";

    auto is_dec = is_decoy;
    is_dec.insert(is_dec.begin(), 0.5);

    std::vector<double> decoy_rate;
    if (!scores.empty()) {
        if (VERB > 2)
            std::cerr << "[TIMING] choosing fit_xy\n";
        auto sc = scores;
        sc.insert(sc.begin(), sc[0]);
        decoy_rate = regressor_ptr_->fit_xy(sc, is_dec, epsilon_, 1. - epsilon_);
    } else {
        if (VERB > 2)
            std::cerr << "[TIMING] choosing fit_y\n";
        decoy_rate = regressor_ptr_->fit_y(is_dec,  epsilon_, 1. - epsilon_);
    }
    decoy_rate.erase(decoy_rate.begin());

    std::vector<double> pep_iso;
    // auto i_d = is_decoy.begin(); double tar = 0.; double dec = 0.; double errors = 0.;
    for (auto& dp : decoy_rate) {
        if (dp > 1. - epsilon_)
            dp = 1. - epsilon_;
        double pep = dp / (1 - dp);
        if (pep > 1.)
            pep = 1.;
        pep_iso.push_back(pep);
        // if (*i_d>0.5) {dec+=1.;} else {tar+=1.; errors += pep;} i_d++;
        // cerr << dp << "\t" << pep << "\t" << tar << "\t" << dec << "\t" << errors << "\t" << (dec+0.5)/(tar) << "\t" << errors/(tar) << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    if (VERB > 2)
        std::cerr << "[TIMING] tdc_to_pep duration: " << duration << " seconds\n";

    return pep_iso;
}
