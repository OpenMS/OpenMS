// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EGHTRACEFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EGHTRACEFITTER_H

#include <numeric> // for "accumulate"
#include <sstream>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
#include "OpenMS/MATH/GSL_WRAPPER/gsl_wrapper.h"

namespace OpenMS
{

  /**
   * @brief A RT Profile fitter using an Exponential Gaussian Hybrid background model
   *
   * Lan K, Jorgenson JW.
   * <b>A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks.</b>
   * <em>Journal of Chromatography A.</em> 915 (1-2)p. 1-13.
   * Available at: http://linkinghub.elsevier.com/retrieve/pii/S0021967301005945
   *
   * @htmlinclude OpenMS_EGHTraceFitter.parameters
   *
   * @experimental Needs further testing on real data. Note that the tests are currently also focused on testing the EGH as replacement for the gaussian.
   */
  template <class PeakType>
  class EGHTraceFitter :
    public TraceFitter<PeakType>
  {
public:
    /** Functor for LM Optimization */
    class EGHTraceFunctor : public TraceFitter<PeakType>::GenericFunctor
    {
    public:
      EGHTraceFunctor(int dimensions,
          const typename TraceFitter<PeakType>::ModelData* data)
      : TraceFitter<PeakType>::GenericFunctor(dimensions, data->traces_ptr->getPeakCount()), m_data(data) {}

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec)
      {
        double H  = x(0);
        double tR = x(1);
        double sigma_square = x(2);
        double tau = x(3);

        double t_diff, t_diff2, denominator = 0.0;

        double fegh = 0.0;

        UInt count = 0;
        for (Size t = 0; t < m_data->traces_ptr->size(); ++t)
        {
          const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> & trace = m_data->traces_ptr->at(t);
          DoubleReal weight = m_data->weighted ? trace.theoretical_int : 1.0;
          for (Size i = 0; i < trace.peaks.size(); ++i)
          {
            DoubleReal rt = trace.peaks[i].first;

            t_diff = rt - tR;
            t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

            denominator = 2 * sigma_square + tau * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

            if (denominator > 0.0)
            {
              fegh =  m_data->traces_ptr->baseline + trace.theoretical_int * H * exp(-t_diff2 / denominator);
            }
            else
            {
              fegh = 0.0;
            }

            fvec(count) = (fegh - trace.peaks[i].second->getIntensity()) * weight;
            ++count;
          }
        }
        return 0;
      }
      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J)
      {
        double H  = x(0);
        double tR = x(1);
        double sigma_square = x(2);
        double tau = x(3);

        double derivative_H, derivative_tR, derivative_sigma_square, derivative_tau = 0.0;
        double t_diff, t_diff2, exp1, denominator = 0.0;

        UInt count = 0;
        for (Size t = 0; t < m_data->traces_ptr->size(); ++t)
        {
          const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> & trace = m_data->traces_ptr->at(t);
          DoubleReal weight = m_data->weighted ? trace.theoretical_int : 1.0;
          for (Size i = 0; i < trace.peaks.size(); ++i)
          {
            DoubleReal rt = trace.peaks[i].first;

            t_diff = rt - tR;
            t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

            denominator = 2 * sigma_square + tau * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

            if (denominator > 0)
            {
              exp1 = exp(-t_diff2 / denominator);

              // \partial H f_{egh}(t) = \exp\left( \frac{-\left(t-t_R \right)}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right)
              derivative_H = trace.theoretical_int * exp1;

              // \partial t_R f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{\left( 4 \sigma_{g}^{2} + \tau \left(t-t_R \right) \right) \left(t-t_R \right)}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
              derivative_tR = trace.theoretical_int * H * exp1 * (((4 * sigma_square + tau * t_diff) * t_diff) / (denominator * denominator));

              // \partial \sigma_{g}^{2} f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)^2}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{ 2 \left(t - t_R\right)^2}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
              derivative_sigma_square = trace.theoretical_int * H * exp1 * ((2 * t_diff2) / (denominator * denominator));

              // \partial \tau f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)^2}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{ \left(t - t_R\right)^3}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
              derivative_tau = trace.theoretical_int * H * exp1 * ((t_diff * t_diff2) / (denominator * denominator));
            }
            else
            {
              derivative_H = 0.0;
              derivative_tR = 0.0;
              derivative_sigma_square = 0.0;
              derivative_tau = 0.0;
            }

            // set the jacobian matrix
            J(count, 0) = derivative_H * weight;
            J(count, 1) = derivative_tR * weight;
            J(count, 2) = derivative_sigma_square * weight;
            J(count, 3) = derivative_tau * weight;
            ++count;
          }
        }
        return 0;
      }
    protected:
      const typename TraceFitter<PeakType>::ModelData* m_data;
    };

    EGHTraceFitter()
    {
      //setName("EGHTraceFitter");
    }

    EGHTraceFitter(const EGHTraceFitter& other) :
      TraceFitter<PeakType>(other)
    {
      this->height_ = other.height_;
      this->apex_rt_ = other.apex_rt_;
      this->sigma_ = other.sigma_;
      this->tau_ = other.tau_;

      this->sigma_5_bound_ = other.sigma_5_bound_;

      updateMembers_();
    }

    EGHTraceFitter& operator=(const EGHTraceFitter& source)
    {
      TraceFitter<PeakType>::operator=(source);

      this->height_ = source.height_;
      this->apex_rt_ = source.apex_rt_;
      this->sigma_ = source.sigma_;
      this->tau_ = source.tau_;

      this->sigma_5_bound_ = source.sigma_5_bound_;

      updateMembers_();

      return *this;
    }

    virtual ~EGHTraceFitter()
    {
    }

    // override important methods
    void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>& traces)
    {
      setInitialParameters_(traces);

      Eigen::VectorXd x_init(NUM_PARAMS_);
      x_init(0) = height_;
      x_init(1) = apex_rt_;
      x_init(2) = sigma_square_;
      x_init(3) = tau_;

      typename TraceFitter<PeakType>::ModelData data;
      data.traces_ptr = &traces;
      data.weighted = this->weighted_;
      EGHTraceFunctor functor (NUM_PARAMS_, &data);

      TraceFitter<PeakType>::optimize_(x_init, functor);
    }

    DoubleReal getLowerRTBound() const
    {
      return sigma_5_bound_.first;
    }

    DoubleReal getTau() const
    {
      return tau_;
    }

    DoubleReal getUpperRTBound() const
    {
      return sigma_5_bound_.second;
    }

    DoubleReal getHeight() const
    {
      return height_;
    }

    DoubleReal getSigma() const
    {
      return sigma_;
    }

    DoubleReal getCenter() const
    {
      return apex_rt_;
    }

    bool checkMaximalRTSpan(const DoubleReal max_rt_span)
    {
      return (sigma_5_bound_.second - sigma_5_bound_.first) > max_rt_span * region_rt_span_;
    }

    virtual bool checkMinimalRTSpan(const std::pair<DoubleReal, DoubleReal>& rt_bounds, const DoubleReal min_rt_span)
    {
      return (rt_bounds.second - rt_bounds.first) < min_rt_span * (sigma_5_bound_.second - sigma_5_bound_.first);
    }

    DoubleReal computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace, Size k)
    {
      double rt = trace.peaks[k].first;
      double t_diff, t_diff2, denominator = 0.0;
      double fegh = 0.0;

      t_diff = rt - apex_rt_;
      t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

      denominator = 2 * sigma_ * sigma_ + tau_ * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

      if (denominator > 0.0)
      {
        fegh =  trace.theoretical_int * height_ * exp(-t_diff2 / denominator);
      }

      return fegh;
    }

    virtual DoubleReal getArea()
    {
      // equation 21 from Lan & Jorgenson paper:
      DoubleReal abs_tau = fabs(tau_);
      DoubleReal phi = atan(abs_tau / sigma_);
      DoubleReal epsilon = EPSILON_COEFS[0];
      DoubleReal phi_pow = phi;
      for (Size i = 1; i < 7; ++i) {
        epsilon += phi_pow * EPSILON_COEFS[i];
        phi_pow *= phi;
      }
      // 0.62... is approx. sqrt(pi / 8):
      return height_ * (sigma_ * 0.6266571 + abs_tau) * epsilon;
    }

    DoubleReal getFWHM() const
    {
      std::pair<DoubleReal, DoubleReal> bounds = getAlphaBoundaries_(0.5);
      return bounds.second - bounds.first;
    }

    virtual String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace, const char function_name, const DoubleReal baseline, const DoubleReal rt_shift)
    {
      std::stringstream s;
      s << String(function_name)  << "(x)= " << baseline << " + ";
      s << "("; // the overall bracket
      s << "((" << 2 * sigma_ * sigma_ << " + " << tau_ << " * (x - " << (rt_shift + apex_rt_) << " )) > 0) ? "; // condition
      s <<  (trace.theoretical_int *  height_) << " * exp(-1 * (x - " << (rt_shift + apex_rt_) << ")**2 " <<
      "/" <<
      " ( " << 2 * sigma_ * sigma_ << " + " << tau_ << " * (x - " << (rt_shift + apex_rt_) << " )))";
      s << " : 0)";
      return String(s.str());
    }

protected:
    DoubleReal apex_rt_;
    DoubleReal height_;

    DoubleReal sigma_;
    DoubleReal tau_;

    std::pair<DoubleReal, DoubleReal> sigma_5_bound_;

    DoubleReal region_rt_span_;

    /// Coefficients to calculate the proportionality factor for the peak area
    static const DoubleReal EPSILON_COEFS[];

    static const Size NUM_PARAMS_ = 4;

    /**
     * @brief Return an ordered pair of the positions where the EGH reaches a height of alpha * height of the EGH
     *
     * @param alpha The alpha at which the boundaries should be computed
     */
    std::pair<DoubleReal, DoubleReal> getAlphaBoundaries_(const DoubleReal alpha) const
    {
      std::pair<DoubleReal, DoubleReal> bounds;
      // solved equations A.2 and A.3 from the Lan & Jorgenson paper (Appendix
      // A) for the boundaries A_alpha and B_alpha:
      DoubleReal L = log(alpha);
      DoubleReal s = sqrt(((L * tau_) * (L * tau_) / 4) - 2 * L * sigma_ * sigma_);

      DoubleReal s1, s2;
      s1 = (-1 * (L * tau_) / 2) + s;
      s2 = (-1 * (L * tau_) / 2) - s;

      // the smaller one (should be < 0) = lower bound
      bounds.first = apex_rt_ + std::min(s1, s2);
      // bigger one (should be > 0) = upper bound
      bounds.second = apex_rt_ + std::max(s1, s2);

      return bounds;
    }

    void getOptimizedParameters_(const Eigen::VectorXd& x_init)
    {
      height_ =  x_init(0);
      apex_rt_ =  x_init(1);
      sigma_square_ =  x_init(2);
      tau_ =  x_init(3);

      // we set alpha to 0.04 which is conceptually equal to
      // 2.5 sigma for lower and upper bound
      sigma_5_bound_ = getAlphaBoundaries_(0.043937);
    }

    void setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>& traces)
    {
      LOG_DEBUG << "EGHTraceFitter->setInitialParameters(...)" << std::endl;
      LOG_DEBUG << "Number of traces: " << traces.size() << std::endl;

      // aggregate data; some peaks (where intensity is zero) can be missing!
      // mapping: RT -> total intensity over all mass traces
      std::map<DoubleReal, DoubleReal> total_intensities;
      for (typename FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>::iterator t_it = traces.begin(); t_it != traces.end(); ++t_it)
      {
        for (typename std::vector<std::pair<DoubleReal, const PeakType*> >::iterator p_it = t_it->peaks.begin(); p_it != t_it->peaks.end(); ++p_it)
        {
          total_intensities[p_it->first] += p_it->second->getIntensity();
        }
      }

      // compute moving average for smoothing:
      const Size N = total_intensities.size();
      const Size LEN = 2; // window size: 2 * LEN + 1
      std::vector<DoubleReal> totals(N + 2 * LEN); // pad with zeros at ends
      Int index = LEN;
      // LOG_DEBUG << "Summed intensities:\n";
      for (std::map<DoubleReal, DoubleReal>::iterator it =
             total_intensities.begin(); it != total_intensities.end(); ++it)
      {
        totals[index++] = it->second;
        // LOG_DEBUG << it->second << std::endl;
      }
      std::vector<DoubleReal> smoothed(N);
      Size max_index = 0; // index of max. smoothed intensity
      // LOG_DEBUG << "Smoothed intensities:\n";
      DoubleReal sum = std::accumulate(&totals[LEN], &totals[2 * LEN], 0.0);
      for (Size i = 0; i < N; ++i)
      {
        sum += totals[i + 2 * LEN];
        smoothed[i] = sum / (2 * LEN + 1);
        sum -= totals[i];
        if (smoothed[i] > smoothed[max_index]) max_index = i;
        // LOG_DEBUG << smoothed[i] << std::endl;
      }
      LOG_DEBUG << "Maximum at index " << max_index << std::endl;
      height_ = smoothed[max_index] - traces.baseline;
      LOG_DEBUG << "height: " << height_ << std::endl;
      std::map<DoubleReal, DoubleReal>::iterator it = total_intensities.begin();
      std::advance(it, max_index);
      apex_rt_ = it->first;
      LOG_DEBUG << "apex_rt: " << apex_rt_ << std::endl;
      region_rt_span_ = (total_intensities.rbegin()->first - 
                         total_intensities.begin()->first);
      LOG_DEBUG << "region_rt_span: " << region_rt_span_ << std::endl;

      // find RT values where intensity is at half-maximum:
      index = max_index;
      while ((index > 0) && (smoothed[index] > height_ * 0.5)) --index;
      DoubleReal left_height = smoothed[index];
      it = total_intensities.begin();
      std::advance(it, index);
      DoubleReal left_rt = it->first;
      LOG_DEBUG << "Left half-maximum at index " << index << ", RT " << left_rt
                << std::endl;
      index = max_index;
      while ((index < Int(N - 1)) && (smoothed[index] > height_ * 0.5)) ++index;
      DoubleReal right_height = smoothed[index];
      it = total_intensities.end();
      std::advance(it, index - Int(N));
      DoubleReal right_rt = it->first;
      LOG_DEBUG << "Right half-maximum at index " << index << ", RT "
                << right_rt << std::endl;

      DoubleReal A = apex_rt_ - left_rt;
      DoubleReal B = right_rt - apex_rt_;
      //LOG_DEBUG << "A: " << A << std::endl;
      //LOG_DEBUG << "B: " << B << std::endl;

      // compute estimates for tau / sigma based on A and B:
      DoubleReal alpha = (left_height + right_height) * 0.5 / height_; // ~0.5
      DoubleReal log_alpha = log(alpha);

      tau_ = -1 / log_alpha * (B - A);
      //EGH function fails when tau==0
      if(tau_ == 0)
        tau_ = std::numeric_limits<double>::epsilon();
      LOG_DEBUG << "tau: " << tau_ << std::endl;
      sigma_ = sqrt(-0.5 / log_alpha * B * A);
      LOG_DEBUG << "sigma: " << sigma_ << std::endl;
    }

    virtual void updateMembers_()
    {
      TraceFitter<PeakType>::updateMembers_();
    }

  };

  // from table 1 in the Lan & Jorgenson paper:
  template <class PeakType>
  const DoubleReal EGHTraceFitter<PeakType>::EPSILON_COEFS[] = 
  {4.0, -6.293724, 9.232834, -11.342910, 9.123978, -4.173753, 0.827797};

} // namespace OpenMS

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDTRACEFITTERGAUSS_H
