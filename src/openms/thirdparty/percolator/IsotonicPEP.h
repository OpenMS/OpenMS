#ifndef ISOTONIC_REGRESSION_H_
#define ISOTONIC_REGRESSION_H_

#include "Globals.h"
#include <chrono>
#include <iostream> // for std::cerr, std::endl
#include <memory>   // for std::unique_ptr, std::make_unique
#include <algorithm> // for std::min, std::max

#ifdef _MSC_VER
  #define NOMINMAX
  #include <float.h>
  #define _finite(x) _finite(x)
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse> 

// Constants for tuning (Point 7)
constexpr int DEFAULT_NUM_BINS = 10000;
constexpr double DEFAULT_LAMBDA = 1e-6;
constexpr int DEFAULT_NUM_KNOTS = 50;
constexpr double DEFAULT_SKEW_FACTOR = 0.75;
constexpr double DEFAULT_EPSILON = 1e-20;

template <typename T>
const T& clamp(const T& v, const T& lo, const T& hi) {
    return (v < lo) ? lo : (hi < v) ? hi : v;
}

class IsotonicRegression {
public:
    IsotonicRegression(double epsilon = DEFAULT_EPSILON) : epsilon_(epsilon) {};
    virtual ~IsotonicRegression() = default;

    virtual std::vector<double> fit_y(const std::vector<double>& y, 
                                      double min_val = std::numeric_limits<double>::lowest(),
                                      double max_val = std::numeric_limits<double>::max()) const = 0;

    virtual std::vector<double> fit_xy(const std::vector<double>& x, const std::vector<double>& y,
                                       double min_val = std::numeric_limits<double>::lowest(),
                                       double max_val = std::numeric_limits<double>::max()) const = 0;
protected:
    double epsilon_;   
};

class PavaRegression : public IsotonicRegression {
protected:
    //--------------------------------------------------------
    // A simple block structure for merging
    //--------------------------------------------------------
    struct Block {
        double sum   = 0.0;
        int    count = 0;
        double avg   = 0.0;
    };

    virtual std::vector<double> pavaNonDecreasingRanged(
        const std::vector<double>& values,
        const double min_value = std::numeric_limits<double>::min(),
        const double max_value = std::numeric_limits<double>::max()
    ) const;

public:
    std::vector<double> fit_y(const std::vector<double>& y, double min_val, double max_val) const override
    { return pavaNonDecreasingRanged(y, min_val, max_val); }
    std::vector<double> fit_xy(const std::vector<double>&/* x */, const std::vector<double>& y,
                               double min_val, double max_val) const override { return fit_y(y, min_val, max_val); };

};


class IsplineRegression : public IsotonicRegression {
    public:
        IsplineRegression(int num_bins = DEFAULT_NUM_BINS,
            double lambda = DEFAULT_LAMBDA,
            int max_knots = DEFAULT_NUM_KNOTS,
            double skew_factor = DEFAULT_SKEW_FACTOR,
            double epsilon = DEFAULT_EPSILON) : 
                IsotonicRegression(epsilon), 
                    num_bins_(num_bins),
                    lambda_(lambda),
                    max_knots_(max_knots),
                    skew_factor_(skew_factor) {}

        std::vector<double> fit_xy(const std::vector<double>& x, const std::vector<double>& y,
            double min_val = std::numeric_limits<double>::lowest(),
            double max_val = std::numeric_limits<double>::max()) const override;
        std::vector<double> fit_y(const std::vector<double>& y, 
            double min_val = std::numeric_limits<double>::lowest(),
            double max_val = std::numeric_limits<double>::max()) const override;

        struct BinnedData {
            std::vector<double> x, y, weights;
        };
        double cubic_ispline(double x, double left, double right) const;

    protected:
    
        BinnedData bin_data(const std::vector<double>& x, const std::vector<double>& y, int max_bins) const;
        std::vector<double> compute_adaptive_knots(const std::vector<double>& x, const std::vector<double>& y, int num_knots) const;
        Eigen::VectorXd fit_spline(const BinnedData& data, const std::vector<double>& knots, double lambda) const;

        int num_bins_;
        double lambda_;
        int max_knots_;
        double skew_factor_; 
    };
    

class InferPEP {
    public:
        InferPEP(bool use_ispline, double epsilon = DEFAULT_EPSILON) 
            : regressor_ptr_(use_ispline ? std::unique_ptr<IsotonicRegression>(new IsplineRegression()) : 
                                           std::unique_ptr<IsotonicRegression>(new PavaRegression())), 
              epsilon_(epsilon) {
            if (VERB > 1) {
                std::cerr << "Performing isotonic regression using " 
                          << (use_ispline ? "I-Splines" : "PAVA") << std::endl;
            }
        }
        virtual ~InferPEP() = default;        
        std::vector<double> q_to_pep(const std::vector<double>& q_values);
        std::vector<double> qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores);
        std::vector<double> tdc_to_pep(const std::vector<double>& is_decoy, const std::vector<double>& scores = {});

        double interpolate(const double q_value, const double q1, const double q2, const double pep1, const double pep2) const {
            double interp_pep = pep1 + (q_value - q1) * (pep2 - pep1) / (q2 - q1);
            return interp_pep;
        }
        
    private:
        std::unique_ptr<IsotonicRegression> regressor_ptr_;
        std::vector<double> qs;
        std::vector<double> pep_iso;
        double epsilon_;
    };

    #endif /* ISOTONICPEP_H_ */