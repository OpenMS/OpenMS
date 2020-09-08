// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <numeric> // for "accumulate"

namespace OpenMS
{
  const Size GaussTraceFitter::NUM_PARAMS_ = 3;

  GaussTraceFitter::GaussTraceFitter()
  {
    //setName("GaussTraceFitter");
  }

  GaussTraceFitter::GaussTraceFitter(const GaussTraceFitter& other) :
    TraceFitter(other)
  {
    this->height_ = other.height_;
    this->x0_ = other.x0_;
    this->sigma_ = other.sigma_;

    updateMembers_();
  }

  GaussTraceFitter& GaussTraceFitter::operator=(const GaussTraceFitter& source)
  {
    TraceFitter::operator=(source);

    this->height_ = source.height_;
    this->x0_ = source.x0_;
    this->sigma_ = source.sigma_;

    updateMembers_();

    return *this;
  }

  GaussTraceFitter::~GaussTraceFitter()
  {
  }

  void GaussTraceFitter::fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces)
  {
    OPENMS_LOG_DEBUG << "Traces length: " << traces.size() << "\n";
    setInitialParameters_(traces);

    Eigen::VectorXd x_init(NUM_PARAMS_);
    x_init(0) = height_;
    x_init(1) = x0_;
    x_init(2) = sigma_;

    TraceFitter::ModelData data;
    data.traces_ptr = &traces;
    data.weighted = this->weighted_;
    GaussTraceFunctor functor(NUM_PARAMS_, &data);

    TraceFitter::optimize_(x_init, functor);
  }

  double GaussTraceFitter::getLowerRTBound() const
  {
    return x0_ - 2.5 * sigma_;
  }

  double GaussTraceFitter::getUpperRTBound() const
  {
    return x0_ + 2.5 * sigma_;
  }

  double GaussTraceFitter::getHeight() const
  {
    return height_;
  }

  double GaussTraceFitter::getCenter() const
  {
    return x0_;
  }

  double GaussTraceFitter::getFWHM() const
  {
    return 2.35482 * sigma_; // 2 * sqrt(2 * log(2)) * sigma
  }

  double GaussTraceFitter::getSigma() const
  {
    return sigma_;
  }

  bool GaussTraceFitter::checkMaximalRTSpan(const double max_rt_span)
  {
    return 5.0 * sigma_ > max_rt_span * region_rt_span_;
  }

  bool GaussTraceFitter::checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span)
  {
    return (rt_bounds.second - rt_bounds.first) < (min_rt_span * 5.0 * sigma_);
  }

  double GaussTraceFitter::getValue(double rt) const
  {
    return height_ * exp(-0.5 * pow(rt - x0_, 2) / pow(sigma_, 2));
  }

  double GaussTraceFitter::getArea()
  {
    // area under the curve, 2.5... is approx. sqrt(2 * pi):
    return 2.506628 * height_ * sigma_;
  }

  String GaussTraceFitter::getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift)
  {
    std::stringstream s;
    s << String(function_name)  << "(x)= " << baseline << " + ";
    s << (trace.theoretical_int *  height_) << " * exp(-0.5*(x-" << (rt_shift + x0_) << ")**2/(" << sigma_ << ")**2)";
    return String(s.str());
  }

  void GaussTraceFitter::getOptimizedParameters_(const Eigen::VectorXd& x_init)
  {
    height_ = x_init(0);
    x0_ = x_init(1);
    sigma_ = std::fabs(x_init(2));
  }

  GaussTraceFitter::GaussTraceFunctor::GaussTraceFunctor(int dimensions,
                                                         const TraceFitter::ModelData* data) :
    TraceFitter::GenericFunctor(dimensions,
                                static_cast<int>(data->traces_ptr->getPeakCount())),
    m_data(data)
  {
  }

  int GaussTraceFitter::GaussTraceFunctor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec)
  {
    double height = x(0);
    double x0 = x(1);
    double sig = x(2);
    double c_fac = -0.5 / pow(sig, 2);

    Size count = 0;
    for (Size t = 0; t < m_data->traces_ptr->size(); ++t)
    {
      const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace = (*m_data->traces_ptr)[t];
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
  int GaussTraceFitter::GaussTraceFunctor::df(const Eigen::VectorXd& x, Eigen::MatrixXd& J)
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
      const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace = (*m_data->traces_ptr)[t];
      double weight = m_data->weighted ? trace.theoretical_int : 1.0;
      for (Size i = 0; i < trace.peaks.size(); ++i)
      {
        double rt = trace.peaks[i].first;
        double e = exp(c_fac * pow(rt - x0, 2));
        J(count, 0) = trace.theoretical_int * e * weight;
        J(count, 1) = trace.theoretical_int * height * e * (rt - x0) / sig_sq * weight;
        J(count, 2) = 0.125* trace.theoretical_int* height* e* pow(rt - x0, 2) / sig_3 * weight;
        ++count;
      }
    }
    return 0;
  }

  void GaussTraceFitter::setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces)
  {
    OPENMS_LOG_DEBUG << "GaussTraceFitter->setInitialParameters(...)" << std::endl;
    OPENMS_LOG_DEBUG << "Number of traces: " << traces.size() << std::endl;

    // aggregate data; some peaks (where intensity is zero) can be missing!
    // mapping: RT -> total intensity over all mass traces
    std::list<std::pair<double, double> > total_intensities;
    traces.computeIntensityProfile(total_intensities);

    const Size N = total_intensities.size();
    const Size LEN = 2; // window size: 2 * LEN + 1

    std::vector<double> totals(N + 2 * LEN); // pad with zeros at ends
    Int index = LEN;
    // OPENMS_LOG_DEBUG << "Summed intensities:\n";
    for (std::list<std::pair<double, double> >::iterator it =
           total_intensities.begin(); it != total_intensities.end(); ++it)
    {
      totals[index++] = it->second;
      // OPENMS_LOG_DEBUG << it->second << std::endl;
    }

    std::vector<double> smoothed(N);
    Size max_index = 0; // index of max. smoothed intensity
    if (N <= LEN + 1) // not enough distinct x values for smoothing
    {
      // throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-MovingAverage", "Too few time points for smoothing with window size " + String(2 * LEN + 1));
      for (Size i = 0; i < N; ++i)
      {
        smoothed[i] = totals[i + LEN];
        if (smoothed[i] > smoothed[max_index]) max_index = i;
      }
    }
    else // compute moving average for smoothing
    {
      // OPENMS_LOG_DEBUG << "Smoothed intensities:\n";
      double sum = std::accumulate(&totals[LEN], &totals[2 * LEN], 0.0);
      for (Size i = 0; i < N; ++i)
      {
        sum += totals[i + 2 * LEN];
        smoothed[i] = sum / (2 * LEN + 1);
        sum -= totals[i];
        if (smoothed[i] > smoothed[max_index]) max_index = i;
        // OPENMS_LOG_DEBUG << smoothed[i] << std::endl;
      }
    }
    OPENMS_LOG_DEBUG << "Maximum at index " << max_index << std::endl;
    height_ = smoothed[max_index] - traces.baseline;
    OPENMS_LOG_DEBUG << "height: " << height_ << std::endl;
    std::list<std::pair<double, double> >::iterator it = total_intensities.begin();
    std::advance(it, max_index);
    x0_ = it->first;
    OPENMS_LOG_DEBUG << "x0: " << x0_ << std::endl;
    region_rt_span_ = (total_intensities.rbegin()->first -
                       total_intensities.begin()->first);
    OPENMS_LOG_DEBUG << "region_rt_span: " << region_rt_span_ << std::endl;

    // find RT values where intensity is at half-maximum:
    index = static_cast<Int>(max_index);
    while ((index > 0) && (smoothed[index] > height_ * 0.5))
      --index;
    double left_height = smoothed[index];
    it = total_intensities.begin();
    std::advance(it, index);
    double left_rt = it->first;
    OPENMS_LOG_DEBUG << "Left half-maximum at index " << index << ", RT " << left_rt
              << std::endl;
    index = static_cast<Int>(max_index);
    while ((index < Int(N - 1)) && (smoothed[index] > height_ * 0.5))
      ++index;
    double right_height = smoothed[index];
    it = total_intensities.end();
    std::advance(it, index - Int(N));
    double right_rt = it->first;
    OPENMS_LOG_DEBUG << "Right half-maximum at index " << index << ", RT "
              << right_rt << std::endl;

    double delta_x = right_rt - left_rt;
    double alpha = (left_height + right_height) * 0.5 / height_; // ~0.5
    if (alpha >= 1) sigma_ = 1.0; // degenerate case, all values are the same
    else sigma_ = delta_x * 0.5 / sqrt(-2.0 * log(alpha));
    OPENMS_LOG_DEBUG << "sigma: " << sigma_ << std::endl;
  }

  void GaussTraceFitter::updateMembers_()
  {
    TraceFitter::updateMembers_();
  }

} // namespace OpenMS
