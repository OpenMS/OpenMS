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
#include <numeric> // for "accumulate"

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

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
      LOG_DEBUG << "Traces length: " << traces.size() << std::endl;
      setInitialParameters_(traces);

      double x_init[NUM_PARAMS_] = {height_, x0_, sigma_};

      Size num_params = NUM_PARAMS_;

      TraceFitter<PeakType>::optimize_(traces, num_params, x_init,
                                       &(GaussTraceFitter<PeakType>::residual_),
                                       &(GaussTraceFitter<PeakType>::jacobian_),
                                       &(GaussTraceFitter<PeakType>::evaluate_));
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
      return 2.35482 * sigma_; // 2 * sqrt(2 * log(2)) * sigma
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

    DoubleReal getValue(DoubleReal rt) const
    {
      return height_ * exp(-0.5 * pow(rt - x0_, 2) / pow(sigma_, 2));
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

    void getOptimizedParameters_(gsl_multifit_fdfsolver* s)
    {
      height_ = gsl_vector_get(s->x, 0);
      x0_ = gsl_vector_get(s->x, 1);
      sigma_ = std::fabs(gsl_vector_get(s->x, 2));
    }

    static Int residual_(const gsl_vector* param, void* data, gsl_vector* f)
    {
      typename TraceFitter<PeakType>::ModelData* model_data = static_cast<typename TraceFitter<PeakType>::ModelData*>(data);
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>* traces = model_data->traces_ptr;

      double height = gsl_vector_get(param, 0);
      double x0 = gsl_vector_get(param, 1);
      double sig = gsl_vector_get(param, 2);
      double c_fac = -0.5 / pow(sig, 2);

      Size count = 0;
      for (Size t = 0; t < traces->size(); ++t)
      {
        const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace = (*traces)[t];
        DoubleReal weight = model_data->weighted ? trace.theoretical_int : 1.0;
        for (Size i = 0; i < trace.peaks.size(); ++i)
        {
          DoubleReal fgauss = traces->baseline + trace.theoretical_int * height * exp(c_fac * pow(trace.peaks[i].first - x0, 2));
          gsl_vector_set(f, count, (fgauss - trace.peaks[i].second->getIntensity()) * weight);
          ++count;
        }
      }
      return GSL_SUCCESS;
    }

    static Int jacobian_(const gsl_vector* param, void* data, gsl_matrix* J)
    {
      typename TraceFitter<PeakType>::ModelData* model_data = static_cast<typename TraceFitter<PeakType>::ModelData*>(data);
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>* traces = model_data->traces_ptr;

      double height = gsl_vector_get(param, 0);
      double x0 = gsl_vector_get(param, 1);
      double sig = gsl_vector_get(param, 2);
      double sig_sq = pow(sig, 2);
      double sig_3 = pow(sig, 3);
      double c_fac = -0.5 / sig_sq;

      Size count = 0;
      for (Size t = 0; t < traces->size(); ++t)
      {
        const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace = (*traces)[t];
        DoubleReal weight = model_data->weighted ? trace.theoretical_int : 1.0;
        for (Size i = 0; i < trace.peaks.size(); ++i)
        {
          DoubleReal rt = trace.peaks[i].first;
          DoubleReal e = exp(c_fac * pow(rt - x0, 2));
          gsl_matrix_set(J, count, 0, trace.theoretical_int * e * weight);
          gsl_matrix_set(J, count, 1, trace.theoretical_int * height * e * (rt - x0) / sig_sq * weight);
          gsl_matrix_set(J, count, 2, 0.125 * trace.theoretical_int * height * e * pow(rt - x0, 2) / sig_3 * weight);
          ++count;
        }
      }
      return GSL_SUCCESS;
    }

    static Int evaluate_(const gsl_vector* param, void* data, gsl_vector* f, gsl_matrix* J)
    {
      residual_(param, data, f);
      jacobian_(param, data, J);
      return GSL_SUCCESS;
    }

    void setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>& traces)
    {
      LOG_DEBUG << "GaussTraceFitter->setInitialParameters(...)" << std::endl;
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
      if (N <= LEN + 1) // not enough distinct x values for smoothing
      {
        // throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-MovingAverage", "Too few time points for smoothing with window size " + String(2 * LEN + 1));
        for (Size i = 0; i < N; ++i)
        {
          smoothed[i] = totals[i + LEN];
          if (smoothed[i] > smoothed[max_index]) max_index = i;
        }
      }
      else // compute moving average for smoothing
      {       
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
      }
      LOG_DEBUG << "Maximum at index " << max_index << std::endl;
      height_ = smoothed[max_index] - traces.baseline;
      LOG_DEBUG << "height: " << height_ << std::endl;
      std::map<DoubleReal, DoubleReal>::iterator it = total_intensities.begin();
      std::advance(it, max_index);
      x0_ = it->first;
      LOG_DEBUG << "x0: " << x0_ << std::endl;
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
      
      DoubleReal delta_x = right_rt - left_rt;
      DoubleReal alpha = (left_height + right_height) * 0.5 / height_; // ~0.5
      if (alpha >= 1) sigma_ = 1.0; // degenerate case, all values are the same
      else sigma_ = delta_x * 0.5 / sqrt(-2.0 * log(alpha));
      LOG_DEBUG << "sigma: " << sigma_ << std::endl;
    }

    virtual void updateMembers_()
    {
      TraceFitter<PeakType>::updateMembers_();
    }

    void printState_(SignedSize iter, gsl_multifit_fdfsolver* s)
    {
      LOG_DEBUG << "iter: " << iter << " " 
                << "height: " << gsl_vector_get(s->x, 0) << " "
                << "x0: " << gsl_vector_get(s->x, 1) << " "
                << "sigma: " << std::fabs(gsl_vector_get(s->x, 2)) << " "
                << "|f(x)| = " << gsl_blas_dnrm2(s->f) << std::endl;
    }

  };

} // namespace OpenMS

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDTRACEFITTERGAUSS_H
