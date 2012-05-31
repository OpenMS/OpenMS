// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EGHTRACEFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EGHTRACEFITTER_H

#include <sstream>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

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
  template<class PeakType>
  class EGHTraceFitter
    : public TraceFitter<PeakType>
  {
  public:
    EGHTraceFitter()
    {
      //setName("EGHTraceFitter");
    }

    EGHTraceFitter(const EGHTraceFitter& other)
      : TraceFitter<PeakType>(other)
    {
      this->height_ = other.height_;
      this->apex_rt_ = other.apex_rt_;
      this->sigma_square_ = other.sigma_square_;
      this->tau_ = other.tau_;

      this->sigma_5_bound_ = other.sigma_5_bound_;
      this->fwhm_bound_ = other.fwhm_bound_;

      updateMembers_();
    }

    EGHTraceFitter& operator = (const EGHTraceFitter& source)
    {
      TraceFitter<PeakType>::operator = (source);

      this->height_ = source.height_;
      this->apex_rt_ = source.apex_rt_;
      this->sigma_square_ = source.sigma_square_;
      this->tau_ = source.tau_;

      this->sigma_5_bound_ = source.sigma_5_bound_;
      this->fwhm_bound_ = source.fwhm_bound_;

      updateMembers_();

      return *this;
    }

    virtual ~EGHTraceFitter()
    {
    }


    // override important methods
    void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType> & traces)
    {
      setInitialParameters_(traces);

      double x_init[NUM_PARAMS_] = {height_, apex_rt_, sigma_square_, tau_};

      Size num_params = NUM_PARAMS_;

      TraceFitter<PeakType>::optimize_(traces, num_params , x_init ,
          &( EGHTraceFitter<PeakType>::residual_ ) ,
          &( EGHTraceFitter<PeakType>::jacobian_ ) ,
          &( EGHTraceFitter<PeakType>::evaluate_ ));
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

    DoubleReal getSigmaSquare() const
    {
      return sigma_square_;
    }

    DoubleReal getCenter() const
    {
      return apex_rt_;
    }

    bool checkMaximalRTSpan(const DoubleReal max_rt_span)
    {
      return ( (sigma_5_bound_.second - sigma_5_bound_.first) > max_rt_span*region_rt_span_);
    }

    virtual bool checkMinimalRTSpan(const std::pair<DoubleReal,DoubleReal> & rt_bounds, const DoubleReal min_rt_span)
    {
      return (rt_bounds.second-rt_bounds.first) < min_rt_span * (sigma_5_bound_.second - sigma_5_bound_.first);
    }

    DoubleReal computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> & trace, Size k)
    {
      double rt = trace.peaks[k].first;
      double t_diff,t_diff2,denominator = 0.0;
      double fegh = 0.0;

      t_diff = rt - apex_rt_;
      t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

      denominator = 2 * sigma_square_ + tau_ * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

      if(denominator > 0.0)
      {
        fegh =  trace.theoretical_int * height_ * exp(- t_diff2 / denominator);
      }

      return fegh;
    }

    virtual DoubleReal getFeatureIntensityContribution()
    {
      return height_ * (fwhm_bound_.second - fwhm_bound_.first);
    }

    DoubleReal getFWHM() const
    {

      std::pair<DoubleReal, DoubleReal> bounds = getAlphaBoundaries_(0.5);
      return bounds.second - bounds.first;
    }

    virtual String getGnuplotFormula(FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> const & trace, const char function_name, const DoubleReal baseline, const DoubleReal rt_shift)
    {
      std::stringstream s;
      s << String(function_name)  << "(x)= " << baseline << " + ";
      s << "("; // the overall bracket
      s << "((" << 2 * sigma_square_ << " + " << tau_ << " * (x - " << (rt_shift + apex_rt_) << " )) > 0) ? "; // condition
      s <<  (trace.theoretical_int *  height_) << " * exp(-1 * (x - " << (rt_shift + apex_rt_) << ")**2 " <<
          "/" <<
          " ( " << 2 * sigma_square_ << " + " << tau_ << " * (x - " << (rt_shift + apex_rt_) << " )))";
      s << " : 0)";
      return String(s.str());
    }

  protected:
    DoubleReal apex_rt_;
    DoubleReal height_;

    DoubleReal sigma_square_;
    DoubleReal tau_;

    std::pair<DoubleReal, DoubleReal> sigma_5_bound_;
    std::pair<DoubleReal, DoubleReal> fwhm_bound_;

    DoubleReal region_rt_span_;

    static const Size NUM_PARAMS_ = 4;

    /**
     * @brief Return an ordered pair of the positions where the EGH reaches a height of alpha * height of the EGH
     *
     * @param alpha The alpha at which the boundaries should be computed
     */
    std::pair<DoubleReal, DoubleReal> getAlphaBoundaries_(const DoubleReal alpha) const
    {
      std::pair<DoubleReal, DoubleReal> bounds;
      DoubleReal L = log(alpha);
      DoubleReal s = sqrt(
          ((L*tau_)*(L*tau_) / 4) - 2*L*sigma_square_
      );

      DoubleReal s1,s2;
      s1 = (-1 * (L * tau_) / 2) + s;
      s2 = (-1 * (L * tau_) / 2) - s;

      // the smaller one (should be < 0) = lower bound
      bounds.first = apex_rt_ + std::min(s1,s2);
      // bigger one (should be > 0) = upper bound
      bounds.second = apex_rt_ + std::max(s1,s2);

      return bounds;
    }

    void getOptimizedParameters_(gsl_multifit_fdfsolver * fdfsolver)
    {
      height_ =  gsl_vector_get(fdfsolver->x, 0);
      apex_rt_ =  gsl_vector_get(fdfsolver->x, 1);
      sigma_square_ =  gsl_vector_get(fdfsolver->x, 2);
      tau_ =  gsl_vector_get(fdfsolver->x, 3);

      // we set alpha to 0.04 which is conceptually equal to
      // 2.5 sigma for lower and upper bound
      sigma_5_bound_ = getAlphaBoundaries_(0.043937);
      // this is needed for the intensity contribution -> this is the 1.25 sigma region
      fwhm_bound_ = getAlphaBoundaries_(0.45783);
    }

    static Int residual_(const gsl_vector* param, void* data, gsl_vector* f)
    {
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>* traces = static_cast<FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>*>(data);

      double H  = gsl_vector_get( param, 0 );
      double tR = gsl_vector_get( param, 1 );
      double sigma_square = gsl_vector_get( param, 2 );
      double tau = gsl_vector_get( param, 3 );

      double t_diff,t_diff2,denominator = 0.0;

      double fegh = 0.0;

      UInt count = 0;
      for (Size t = 0 ; t < traces->size() ; ++t)
      {
        FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace = traces->at(t);
        for (Size i = 0 ; i < trace.peaks.size() ; ++i)
        {
          DoubleReal rt = trace.peaks[i].first;

          t_diff = rt - tR;
          t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

          denominator = 2 * sigma_square + tau * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

          if(denominator > 0.0)
          {
            fegh =  traces->baseline + trace.theoretical_int * H * exp(- t_diff2 / denominator);
          }
          else
          {
            fegh = 0.0;
          }

          gsl_vector_set( f, count, ( fegh - trace.peaks[i].second->getIntensity() ) );
          ++count;
        }
      }
      return GSL_SUCCESS;
    }

    static Int jacobian_(const gsl_vector* param, void* data, gsl_matrix* J)
    {
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>* traces = static_cast<FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>*>(data);

      double H  = gsl_vector_get( param, 0 );
      double tR = gsl_vector_get( param, 1 );
      double sigma_square = gsl_vector_get( param, 2 );
      double tau = gsl_vector_get( param, 3 );

      double derivative_H, derivative_tR, derivative_sigma_square, derivative_tau = 0.0;
      double t_diff,t_diff2,exp1,denominator = 0.0;

      UInt count = 0;
      for (Size t = 0; t < traces->size() ; ++t)
      {
        FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace = traces->at(t);
        for (Size i = 0; i < trace.peaks.size() ; ++i)
        {
          DoubleReal rt = trace.peaks[i].first;

          t_diff = rt - tR;
          t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

          denominator = 2 * sigma_square + tau * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

         if( denominator > 0)
         {
           exp1 = exp(- t_diff2 / denominator);

           // \partial H f_{egh}(t) = \exp\left( \frac{-\left(t-t_R \right)}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right)
           derivative_H = trace.theoretical_int * exp1;

           // \partial t_R f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{\left( 4 \sigma_{g}^{2} + \tau \left(t-t_R \right) \right) \left(t-t_R \right)}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
           derivative_tR = trace.theoretical_int * H * exp1 * ( ((4 * sigma_square + tau * t_diff) * t_diff) / (denominator * denominator));

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
         gsl_matrix_set( J, count, 0, derivative_H );
         gsl_matrix_set( J, count, 1, derivative_tR );
         gsl_matrix_set( J, count, 2, derivative_sigma_square );
         gsl_matrix_set( J, count, 3, derivative_tau );

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

    void setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType> & traces)
    {
      LOG_DEBUG << "EGHTraceFitter->setInitialParameters(..)" << std::endl;
      LOG_DEBUG << "Traces length: " << traces.size() << std::endl;
      LOG_DEBUG << "Max trace: " << traces.max_trace << std::endl;

      // initial values for externals
      height_ = traces[traces.max_trace].max_peak->getIntensity() - traces.baseline;
      LOG_DEBUG << "height: " << height_ << std::endl;
      apex_rt_ = traces[traces.max_trace].max_rt;
      LOG_DEBUG << "apex_rt: " << apex_rt_ << std::endl;
      region_rt_span_ = traces[traces.max_trace].peaks.back().first-traces[traces.max_trace].peaks[0].first;
      LOG_DEBUG << "region_rt_span_: " << region_rt_span_ << std::endl;

      const PeakType * max_peak = traces[traces.max_trace].peaks.begin()->second;
      Size max_pos = 0;

      for (Size i = 1; i < traces[traces.max_trace].peaks.size() ; ++i)
      {
        if (traces[traces.max_trace].peaks[i].second->getIntensity()>max_peak->getIntensity())
        {
          max_peak = traces[traces.max_trace].peaks[i].second;
          max_pos = i;
        }
      }

      Size i = max_pos;
      LOG_DEBUG << "max_pos: " << max_pos << std::endl;
      if (traces[traces.max_trace].peaks.size() < 3)
      {
        // TODO: abort the whole thing here??
        //       because below we REQUIRE at least three peaks!!!
      }

      Size filter_max_pos = traces[traces.max_trace].peaks.size() - 2;

      // compute a smoothed value for the maxima
      // if the maximum is close to the borders, we need to think of something...
      DoubleReal smoothed_height;
      if ((max_pos < 2) || (max_pos+2 >= traces[traces.max_trace].peaks.size()))
      {
        // ... too close to border... no smoothing
        smoothed_height = traces[traces.max_trace].peaks[max_pos].second->getIntensity();
        // TODO: does this trace even make sense?! why wasn't it extended it further? or should we have skipped it beforehand?
      }
      else
      {
        smoothed_height = (traces[traces.max_trace].peaks[max_pos - 2].second->getIntensity()
                            + traces[traces.max_trace].peaks[max_pos - 1].second->getIntensity()
                            + traces[traces.max_trace].peaks[max_pos].second->getIntensity()
                            + traces[traces.max_trace].peaks[max_pos + 1].second->getIntensity()
                            + traces[traces.max_trace].peaks[max_pos + 2].second->getIntensity() ) / 5.0;
      }

      // use  moving average filter to avoid bad initial values
      // moving average of size 5
      // TODO: optimize windows size
      while(i > 2 && i < filter_max_pos)
      {
        // compute smoothed
        DoubleReal smoothed = (traces[traces.max_trace].peaks[i - 2].second->getIntensity()
                               + traces[traces.max_trace].peaks[i - 1].second->getIntensity()
                               + traces[traces.max_trace].peaks[i].second->getIntensity()
                               + traces[traces.max_trace].peaks[i + 1].second->getIntensity()
                               + traces[traces.max_trace].peaks[i + 2].second->getIntensity() ) / 5.0;

        if(smoothed / smoothed_height < 0.5) break;
        else --i;
      }
      LOG_DEBUG << "Left alpha at " << i << " with " << traces[traces.max_trace].peaks[i].first << std::endl;
      double A = apex_rt_ - traces[traces.max_trace].peaks[i].first;

      i = max_pos;
      while(i < filter_max_pos && i > 2)
      {
        DoubleReal smoothed = (traces[traces.max_trace].peaks[i - 2].second->getIntensity()
                               + traces[traces.max_trace].peaks[i - 1].second->getIntensity()
                               + traces[traces.max_trace].peaks[i].second->getIntensity()
                               + traces[traces.max_trace].peaks[i + 1].second->getIntensity()
                               + traces[traces.max_trace].peaks[i + 2].second->getIntensity() ) / 5.0;

        if(smoothed / smoothed_height < 0.5) break;
        else ++i;
      }
      LOG_DEBUG << "Right alpha at " << i << " with " << traces[traces.max_trace].peaks[i].first << std::endl;
      double B = traces[traces.max_trace].peaks[i].first - apex_rt_;

      //LOG_DEBUG << "A: " << A << std::endl;
      //LOG_DEBUG << "B: " << B << std::endl;

      // compute estimates for tau / sigma_square based on A/B
      double log_alpha = log(0.5);

      tau_ = ( -1 / log_alpha ) * (B - A);
      LOG_DEBUG << "tau: " << tau_ << std::endl;
      sigma_square_ = ( -1 / (2 * log_alpha) ) * (B * A);
      LOG_DEBUG << "sigma_square: " << sigma_square_ << std::endl;
    }

    virtual void updateMembers_()
    {
      TraceFitter<PeakType>::updateMembers_();
    }

    void printState_(SignedSize iter, gsl_multifit_fdfsolver * s )
    {
      LOG_DEBUG << "iter: " << iter << " "
          << "height: " << gsl_vector_get( s->x, 0 ) << " "
          << "apex_rt: " << gsl_vector_get( s->x, 1 ) << " "
          << "sigma_square: " << gsl_vector_get( s->x, 2 ) << " "
          << "tau: " << gsl_vector_get( s->x, 3 ) << " "
          << "|f(x)| = " << gsl_blas_dnrm2( s->f ) << std::endl;
    }

  };

} // namespace OpenMS

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDTRACEFITTERGAUSS_H
