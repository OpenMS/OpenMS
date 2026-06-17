#ifndef ISOTONIC_REGRESSION_H_
#define ISOTONIC_REGRESSION_H_

#include <chrono>
#include <iostream> // for std::cerr, std::endl
#include <memory>   // for std::unique_ptr, std::make_unique
#include <algorithm>
#include <cmath>
#include "MonotoneRegressor.h"

namespace OpenMS { namespace Internal { namespace Percolator {

class MonotoneRegressor;

class InferPEP {
    public:
        InferPEP(bool use_ispline);
    
        std::vector<double> q_to_pep(const std::vector<double>& q_values);
        // Calibrate PEPs from q-values (qns) or from target/decoy labels (TDC)
        std::vector<double> qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores);
        std::vector<double> tdc_to_pep(const std::vector<double>& is_decoy, const std::vector<double>& scores = {});

        double interpolate(const double q_value, const double q1, const double q2, const double pep1, const double pep2) const {
            if (std::abs(q2 - q1) <= 1e-12) {
                return std::max(0.0, std::min(1.0, pep2));
            }
            double interp_pep = pep1 + (q_value - q1) * (pep2 - pep1) / (q2 - q1);
            return std::max(0.0, std::min(1.0, interp_pep));
        }
        
    private:
        std::unique_ptr<MonotoneRegressor> regressor_ptr_;
        std::vector<double> qs;
        // std::vector<double> pep_iso;
    };

}}}  // namespace OpenMS::Internal::Percolator

    #endif /* ISOTONICPEP_H_ */
