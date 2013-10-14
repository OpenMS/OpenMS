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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Marc Sturm$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSTRACEFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSTRACEFITTER_H

#include <sstream>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
#include "OpenMS/MATH/GSL_WRAPPER/gsl_wrapper.h"

namespace OpenMS
{

  /**
   * @brief Fitter for RT profiles using a Gaussian background model
   *
   * @htmlinclude OpenMS_GaussTraceFitter.parameters
   *
   * @todo More docu
   */
  template <typename PeakType>
  class GaussTraceFitter :
    public TraceFitter<PeakType>
  {
public:
    GaussTraceFitter()
    {
      //setName("GaussTraceFitter");
    }

    GaussTraceFitter(const GaussTraceFitter& other) :
      TraceFitter<PeakType>(other)
    {
      this->height_ = other.height_;
      this->x0_ = other.x0_;
      this->sigma_ = other.sigma_;

      updateMembers_();
    }

    GaussTraceFitter& operator=(const GaussTraceFitter& source)
    {
      TraceFitter<PeakType>::operator=(source);

      this->height_ = source.height_;
      this->x0_ = source.x0_;
      this->sigma_ = source.sigma_;

      updateMembers_();

      return *this;
    }

    virtual ~GaussTraceFitter()
    {
    }

    // override important methods
    void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>& traces)
    {
      LOG_DEBUG << "Traces length: " << traces.size() << "\n";
      setInitialParameters_(traces);

      Eigen::VectorXd x_init(NUM_PARAMS_);
      x_init(0) = height_;
      x_init(1) = x0_;
      x_init(2) = sigma_;

      typename TraceFitter<PeakType>::ModelData data;
      data.traces_ptr = &traces;
      data.weighted = this->weighted_;
      GaussTraceFunctor functor (NUM_PARAMS_, &data);

      TraceFitter<PeakType>::optimize_(x_init, functor);
    }

    DoubleReal getLowerRTBound() const
    {
      return x0_ - 2.5 * sigma_;
    }

    DoubleReal getUpperRTBound() const
    {
      return x0_ + 2.5 * sigma_;
    }

    DoubleReal getHeight() const
    {
      return height_;
    }

    DoubleReal getCenter() const
    {
      return x0_;
    }

    DoubleReal getFWHM() const
    {
      return 2.0 * sigma_;
    }

    /**
     * @brief Returns the sigma of the fitted gaussian model
     */
    DoubleReal getSigma() const
    {
      return sigma_;
    }

    bool checkMaximalRTSpan(const DoubleReal max_rt_span)
    {
      return 5.0 * sigma_ > max_rt_span * region_rt_span_;
    }

    bool checkMinimalRTSpan(const std::pair<DoubleReal, DoubleReal>& rt_bounds, const DoubleReal min_rt_span)
    {
      return (rt_bounds.second - rt_bounds.first) < (min_rt_span * 5.0 * sigma_);
    }

    DoubleReal computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace, Size k)
    {
      return trace.theoretical_int *  height_ * exp(-0.5 * pow(trace.peaks[k].first - x0_, 2) / pow(sigma_, 2));
    }

    DoubleReal getArea()
    {
      // area under the curve, 2.5 is approx. sqrt(2 * pi):
      return 2.5 * height_ * sigma_;
    }

    String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace, const char function_name, const DoubleReal baseline, const DoubleReal rt_shift)
    {
      std::stringstream s;
      s << String(function_name)  << "(x)= " << baseline << " + ";
      s << (trace.theoretical_int *  height_) << " * exp(-0.5*(x-" << (rt_shift + x0_) << ")**2/(" << sigma_ << ")**2)";
      return String(s.str());
    }

protected:
    DoubleReal sigma_;
    DoubleReal x0_;
    DoubleReal height_;
    DoubleReal region_rt_span_;

    static const Size NUM_PARAMS_ = 3;

    void getOptimizedParameters_(const Eigen::VectorXd& x_init)
    {
      height_ = x_init(0);
      x0_ = x_init(1);
      sigma_ = std::fabs(x_init(2));
    }

    class GaussTraceFunctor : public TraceFitter<PeakType>::GenericFunctor
    {
    public:
      GaussTraceFunctor(int dimensions,
          const typename TraceFitter<PeakType>::ModelData* data)
      : TraceFitter<PeakType>::GenericFunctor(dimensions, data->traces_ptr->getPeakCount()), m_data(data) {}

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec)
      {
        double height = x(0);
        double x0 = x(1);
        double sig = x(2);
        double c_fac = -0.5 / pow(sig, 2);

        Size count = 0;
        for (Size t = 0; t < m_data->traces_ptr->size(); ++t)
        {
          const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> & trace = (*m_data->traces_ptr)[t];
          double weight = m_data->weighted ? trace.theoretical_int : 1.0;
          for (Size i = 0; i < trace.peaks.size(); ++i)
          {
            fvec(count) = (m_data->traces_ptr->baseline + trace.theoretical_int * height
                * exp(c_fac * pow(trace.peaks[i].first - x0, 2)) - trace.peaks[i].second->getIntensity()) * weight;
            ++count;
          }
        }

        return 0;
      }
      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J)
      {
        double height = x(0);
        double x0 = x(1);
        double sig = x(2);
        double sig_sq = pow(sig, 2);
        double sig_3 = pow(sig, 3);
        double c_fac = -0.5 / sig_sq;

        Size count = 0;
        for (Size t = 0; t < m_data->traces_ptr->size(); ++t)
        {
          const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> & trace = (*m_data->traces_ptr)[t];
          DoubleReal weight = m_data->weighted ? trace.theoretical_int : 1.0;
          for (Size i = 0; i < trace.peaks.size(); ++i)
          {
            DoubleReal rt = trace.peaks[i].first;
            DoubleReal e = exp(c_fac * pow(rt - x0, 2));
            J(count, 0) = trace.theoretical_int * e * weight;
            J(count, 1) = trace.theoretical_int * height * e * (rt - x0) / sig_sq * weight;
            J(count, 2) = 0.125 * trace.theoretical_int * height * e * pow(rt - x0, 2) / sig_3 *weight;
            ++count;
          }
        }
        return 0;
      }
    protected:
      const typename TraceFitter<PeakType>::ModelData* m_data;
    };

    void setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>& traces)
    {
      LOG_DEBUG << "GaussTraceFitter->setInitialParameters(..)" << std::endl;
      LOG_DEBUG << "Traces length: " << traces.size() << std::endl;
      LOG_DEBUG << "Max trace: " << traces.max_trace << std::endl;

      // initial values for externals
      height_ = traces[traces.max_trace].max_peak->getIntensity() - traces.baseline;
      LOG_DEBUG << "height: " << height_ << std::endl;
      x0_ = traces[traces.max_trace].max_rt;
      LOG_DEBUG << "x0: " << x0_ << std::endl;
      region_rt_span_ = traces[traces.max_trace].peaks.back().first - traces[traces.max_trace].peaks[0].first;
      LOG_DEBUG << "region_rt_span_: " << region_rt_span_ << std::endl;
      sigma_ = region_rt_span_ / 20.0;
      LOG_DEBUG << "sigma_: " << sigma_ << std::endl;
    }

    virtual void updateMembers_()
    {
      TraceFitter<PeakType>::updateMembers_();
    }

  };

} // namespace OpenMS

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDTRACEFITTERGAUSS_H
