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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <boost/math/special_functions/acosh.hpp>


#ifdef DEBUG_DECONV
#include <iostream>
#include <fstream>
#endif

namespace OpenMS
{


  const double OptimizePeakDeconvolution::dist_ = 1.003;

  //TODO: the operator() and the df function need heavy refactoring!!!
  struct OPDFunctor
  {
    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

    OPDFunctor(unsigned dimensions, unsigned numDataPoints, const OptimizePeakDeconvolution::Data* data) :
      m_inputs(dimensions), m_values(numDataPoints), m_data(data)
    {}

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec)
    {
      //TODO: holding the parameters to be optimized and additional values in the same vector is
      //      most likely not the best idea. should be split in two vectors.
      //
      // x contains the parameters to be optimized.
      // The first two entries are the left and right width, respectively.They are equal
      // for all peaks. Then the height and position of all peaks are stored.
      //
      // m_data might contain any additional parameters. We handle these using class members
      // instead.
      // The vector f is supposed to contain the result when we return from this function.
      const std::vector<double>& signal = m_data->signal;
      const std::vector<double>& positions = m_data->positions;
      const std::vector<PeakShape>& peaks = m_data->peaks;
      const OptimizationFunctions::PenaltyFactorsIntensity& penalties = m_data->penalties;
      Int charge = m_data->charge;
      double leftwidth = x(0);
      double rightwidth = x(1);
      //double posP1 = x(2);

      // iterate over all points of the signal
      for (Size current_point = 0; current_point < positions.size(); current_point++)
      {
        double computed_signal = 0.;
        double current_position = positions[current_point];
        double experimental_signal = signal[current_point];

        //iterate over all peaks
        for (Size current_peak = 0; current_peak < peaks.size(); current_peak++)
        {
          //Store the current parameters for this peak
          double p_height = x(2 + 2 * current_peak);
          double p_position = x(2 + 2 * current_peak + 1);
          double p_width = (current_position <= p_position) ? leftwidth : rightwidth;

          //is it a Lorentz or a Sech - Peak?
          if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
          {
            computed_signal += p_height / (1. + pow(p_width * (current_position - p_position), 2));
          }
          else // It's a Sech - Peak
          {
            computed_signal += p_height / pow(cosh(p_width * (current_position - p_position)), 2);
          }
        }
        fvec(current_point) = computed_signal - experimental_signal;
      }

      // penalties : especially negative heights have to be penalised
      double penalty = 0.;

      double penalty_pos = penalties.pos;
      double penalty_lwidth = penalties.lWidth;
      double penalty_rwidth = penalties.rWidth;
      double penalty_intensity = penalties.height;


      //iterate over all peaks again to compute the penalties
      for (Size current_peak = 0; current_peak < peaks.size(); current_peak++)
      {
        double p_position = x(2 + 2 * current_peak + 1);
        if (current_peak < peaks.size() - 1)
        {

          double next_p_position  = x(2 + 2 * current_peak + 3);
          // if distance between peaks does not match the peptide mass rule
          if (fabs(fabs(p_position - next_p_position) - 1.003 / charge) > 0.05)
          {
            // penalize it
            penalty +=  penalty_pos * 10000
                       * pow(fabs(fabs(p_position - next_p_position) - 1.003 / charge), 2);
          }
        }
        double old_position = peaks[current_peak].mz_position;
        double old_width_l  = peaks[current_peak].left_width;
        double old_width_r = peaks[current_peak].right_width;
        double old_height = peaks[current_peak].height;

        double p_width_l = x(0);
        double p_width_r = x(1);
        double p_height = x(2 + 2 * current_peak);

        if (p_height <  1)
        {
          penalty += 100000* penalty_intensity* pow(fabs(p_height - old_height), 2);

        }
        if (p_width_l < 0)
        {
          penalty += penalty_lwidth * peaks.size() * 10000 * pow(fabs(p_width_l - old_width_l), 2);
        }
        else if (p_width_l < 1.5)
          penalty += 10000 * pow(fabs(p_width_l - old_width_l), 2);
        if (p_width_r < 0)
        {
          penalty += penalty_rwidth * peaks.size() * 10000 * pow(fabs(p_width_r - old_width_r), 2);
        }
        else if (p_width_r < 1.5)
          penalty += 10000 * pow(fabs(p_width_r - old_width_r), 2);
        if (fabs(old_position - p_position) > 0.1)
        {
          penalty += 10000* penalty_pos* pow(fabs(old_position - p_position), 2);
        }
      }
      fvec(fvec.size() - 1) = penalty;

      return 0;
    }

    // compute Jacobian matrix for the different parameters
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J)
    {
      // For the conventions on x and params c.f. the commentary in residual()
      //
      // The matrix J is supposed to contain the result when we return from this function.
      // Note: Jacobian is expected as follows:
      //                    - each row corresponds to one data point
      //                    - each column corresponds to one parameter
      const std::vector<double>& positions = m_data->positions;
      const std::vector<PeakShape>& peaks = m_data->peaks;
      const OptimizationFunctions::PenaltyFactorsIntensity& penalties = m_data->penalties;
      Int charge = m_data->charge;

      double leftwidth = x(0);
      double rightwidth = x(1);

      //TODO: is the reset needed?
      J.setZero();

      // iterate over all points of the signal
      for (Size current_point = 0; current_point < positions.size(); current_point++)
      {
        double current_position = positions[current_point];

        // iterate over all peaks
        for (Size current_peak = 0; current_peak < peaks.size(); current_peak++)
        {
          //Store the current parameters for this peak
          double p_height = x(2 + 2 * current_peak);
          double p_position = x(2 + 2 * current_peak + 1);
          double p_width = (current_position <= p_position) ? leftwidth : rightwidth;

          //is it a Lorentz or a Sech - Peak?
          if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
          {
            double diff = current_position - p_position;
            double denom_inv = 1. / (1. + pow(p_width * diff, 2));

            double ddl_left
              = (current_position <= p_position) ? -2* p_height* pow(diff, 2) * p_width * pow(denom_inv, 2) : 0;

            double ddl_right
              = (current_position  > p_position) ? -2* p_height* pow(diff, 2) * p_width * pow(denom_inv, 2) : 0;

            // left and right width are the same for all peaks,
            // the sums of the derivations over all peaks are stored in the first two columns
            J(current_point, 0) = J(current_point, 0) + ddl_left;
            J(current_point, 1) = J(current_point, 1) + ddl_right;

            double ddx0    = 2* p_height* pow(p_width, 2) * diff * pow(denom_inv, 2);

            // partial derivation with respect to intensity
            J(current_point, 2 + 2 * current_peak) = denom_inv;

            // partial derivation with respect to the mz-position
            J(current_point, 2 + 2 * current_peak + 1) = ddx0;
          }
          else // It's a Sech - Peak
          {
            double diff      = current_position - p_position;
            double denom_inv = 1. / cosh(p_width * diff);

            // The remaining computations are not stable if denom_inv == 0. In that case, we are far away from the peak
            // and can assume that all derivatives vanish
            double sinh_term = (fabs(denom_inv) < 1e-6) ? 0.0 : sinh(p_width * diff);


            double ddl_left  = (current_position <= p_position)
                               ? -2* p_height* sinh_term* diff* pow(denom_inv, 3) :
                               0;
            double ddl_right = (current_position  > p_position)
                               ? -2* p_height* sinh_term* diff* pow(denom_inv, 3) :
                               0;

            J(current_point, 0) = J(current_point, 0) + ddl_left;
            J(current_point, 1) = J(current_point, 1) + ddl_right;

            double ddx0      = 2* p_height* p_width* sinh_term* pow(denom_inv, 3);

            J(current_point, 2 + 2 * current_peak) = pow(denom_inv, 2);
            J(current_point, 2 + 2 * current_peak + 1) = ddx0;
          }
        }
      }

      /** Now iterate over all peaks again to compute the
       *  penalties.
       */
      for (Size current_peak = 0; current_peak < peaks.size(); current_peak++)
      {


        double penalty_p = 0;
        double p_position = x(2 + 2 * current_peak + 1);
        if (current_peak < peaks.size() - 1)
        {

          double next_p_position  = x(2 + 2 * current_peak + 3);
          // if distance between peaks does not match the peptide mass rule
          if (fabs(fabs(p_position - next_p_position) - 1.003 / charge) > 0.05)
          {
            // penalize it
            penalty_p += penalties.pos * 20000
                         * fabs(fabs(p_position - next_p_position) - 1.003 / charge);

          }
        }
        std::cout << "Eigen penalty_p " << penalty_p << std::endl;
        double p_width_left = x(0);
        double p_width_right = x(1);
        double p_height = x(2 + 2 * current_peak);

        double old_position = peaks[current_peak].mz_position;
        double old_width_left  = peaks[current_peak].left_width;
        double old_width_right = peaks[current_peak].right_width;
        double old_height = peaks[current_peak].height;

        double penalty_h = 0., penalty_l = 0., penalty_r = 0.;
        if (p_height < 1)
        {
          penalty_h += 100000 * 2 * penalties.height * (fabs(p_height) - fabs(old_height));
        }

        if (p_width_left < 0)
        {
          penalty_l += peaks.size() * 2 * penalties.lWidth * 10000 * (fabs(p_width_left - old_width_left));
        }
        else if (p_width_left < 1.5)
          penalty_l += 2 * penalties.lWidth * 10000 * pow(fabs(p_width_left - old_width_left), 2);
        if (p_width_right < 0)
        {
          penalty_r += peaks.size() * 2 * penalties.rWidth * 10000 * (fabs(p_width_right - old_width_right));
        }
        else if (p_width_right < 1.5)
          penalty_r += 2 * penalties.rWidth * 10000 * pow(fabs(p_width_right - old_width_right), 2);
        if (fabs(old_position - p_position) > 0.1)
        {
          penalty_p += 10000 * penalties.pos * 2 * fabs(old_position - p_position);
        }

        J(positions.size(), 2 + 2 * current_peak) = 100 * penalty_h;
        J(positions.size(), 0) = 100 * penalty_l;
        J(positions.size(), 1) = 100 * penalty_r;
        J(positions.size(), 2 + 2 * current_peak + 1) = 100 * penalty_p;
      }
      for (int i = 0; i < J.rows(); ++i)
      {
        for (int j = 0; j < J.cols(); ++j)
          std::cout << J(i, j) << " ";
        std::cout << std::endl;
      }
      std::cout << std::endl;
      return 0;
    }

    const int m_inputs, m_values;
    const OptimizePeakDeconvolution::Data* m_data;
  };



  OptimizePeakDeconvolution::OptimizePeakDeconvolution() :
    DefaultParamHandler("OptimizePeakDeconvolution"), charge_(1)
  {

    defaults_.setValue("max_iteration", 10, "maximal number of iterations for the fitting step");
    defaults_.setValue("eps_abs", 1e-04, "if the absolute error gets smaller than this value the fitting is stopped", ListUtils::create<String>("advanced"));
    defaults_.setValue("eps_rel", 1e-04, "if the relative error gets smaller than this value the fitting is stopped", ListUtils::create<String>("advanced"));

    defaults_.setValue("penalties:left_width", 0.0, "penalty term for the fitting of the left width:" \
                                                    "If the left width gets too broad or negative during the fitting it can be penalized.");
    defaults_.setValue("penalties:right_width", 0.0, "penalty term for the fitting of the right width:" \
                                                     "If the right width gets too broad or negative during the fitting it can be penalized.");
    defaults_.setValue("penalties:height", 0.0, "penalty term for the fitting of the intensity:" \
                                                "If it gets negative during the fitting it can be penalized.");
    defaults_.setValue("penalties:position", 0.0, "penalty term for the fitting of the peak position:" \
                                                  "If the position changes more than 0.5Da during the fitting it can be penalized as well as " \
                                                  "discrepancies of the peptide mass rule.");

    defaults_.setValue("fwhm_threshold", 1.0, "If a peaks is broader than fwhm_threshold, it is assumed that it contains another peaks and an additional peak is added.");

    defaultsToParam_();
  }

  void OptimizePeakDeconvolution::updateMembers_()
  {
    penalties_.rWidth = (float)param_.getValue("penalties:right_width");
    penalties_.lWidth = (float)param_.getValue("penalties:left_width");
    penalties_.height = (float)param_.getValue("penalties:height");
    penalties_.pos    = (float)param_.getValue("penalties:position");

  }

  bool OptimizePeakDeconvolution::optimize(std::vector<PeakShape>& peaks, Data& data)
  {

    if (peaks.empty())
      return true;


#ifdef DEBUG_DECONV
    std::cout << "peaksanzahl:" << peaks.size();
    std::cout << "\tpeaks[0].mz_position:" << peaks[0].mz_position << std::endl;

    for (Size j = 0; j < peaks.size(); ++j)
    {
      std::cout << "\tpeaks[j].mz_position:" << peaks[j].mz_position;
      std::cout << "\tpeaks[j].height:" << peaks[j].height << std::endl;
      std::cout << "\tpeaks[j].left_width:" << peaks[j].left_width;
      std::cout << "\tpeaks[j].right_width:" << peaks[j].right_width << std::endl << std::endl;
    }

    for (Size j = 0; j < data.positions.size(); ++j)
    {
      std::cout << "positions[" << j << "]=" << data.positions[j] << std::endl;
    }

#endif

    // the input peaks are stored in a temporary vector
    std::vector<PeakShape> temp_shapes = peaks;

    Size global_peak_number = 0;

    double min(std::numeric_limits<double>::max());
    Int bestCharge = 0;
    Size bestNumPeaks = 0;
    Eigen::VectorXd bestResult(2 + 2 * data.peaks.size());
    bestResult.setZero();

    // try three different charge states : charge-1, charge, charge +1
    // take the best solution
    Int chargeState = (charge_ - 1 > 1) ? charge_ - 1 : charge_;
    Int firstChargeState = chargeState;
#ifdef DEBUG_DECONV
    std::cout << "charge " << chargeState << " max_charge" << charge_ + 1
              << "\tpeaks.size() " << peaks.size() << std::endl;
#endif
    bestCharge = chargeState;
    bestNumPeaks = peaks.size();
    for (; chargeState < charge_ + 2; ++chargeState)
    {

      setNumberOfPeaks_(data, temp_shapes, chargeState);
      Eigen::VectorXd x_init(2 + 2 * data.peaks.size());
      for (Size i = 0; i < data.peaks.size(); i++)
      {
        x_init(2 + 2 * i) = data.peaks[i].height;
        x_init(3 + 2 * i) = data.peaks[i].mz_position;
      }
      // Initialize the parameters for the optimization
      // all peaks shall have the same width
      double wl = data.peaks[0].left_width;
      double wr = data.peaks[0].right_width;
      if (boost::math::isnan(wl))
      {
        for (Size i = 0; i < data.peaks.size(); ++i)
        {
          data.peaks[i].left_width = 1;
        }
        wl = 1.;
      }
      if (boost::math::isnan(wr))
      {
        for (Size i = 0; i < data.peaks.size(); ++i)
        {
          data.peaks[i].right_width = 1;
        }
        wr = 1.;
      }
      x_init(0) = wl;
      x_init(1) = wr;
      data.penalties = penalties_;
      data.charge = chargeState;
      unsigned numDataPoints = std::max(data.positions.size() + 1, 2 + 2 * data.peaks.size());
      OPDFunctor functor(2, numDataPoints, &data);
      Eigen::LevenbergMarquardt<OPDFunctor> lmSolver(functor);
      Eigen::LevenbergMarquardt<OPDFunctor>::Parameters config;
      config.maxfev = (Int)param_.getValue("max_iteration");
      lmSolver.parameters = config;
      Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

      //the states are poorly documented. after checking the source, we believe that
      //all states except NotStarted, Running and ImproperInputParameters are good
      //termination states.
      if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-OptimizePeakDeconvolution", "Could not fit the curve to the data: Error " + String(status));
      }
      double chi = lmSolver.fnorm;
      if ((chargeState == firstChargeState) || (chi < min))
      {

        bestResult = x_init;
        min = chi;
        bestCharge = chargeState;
        bestNumPeaks = data.peaks.size();
      }
    }
    global_peak_number += bestNumPeaks;
    // iterate over all peaks and store the optimized values in peaks
    if (bestNumPeaks > 0)
    {
      peaks.resize(bestNumPeaks);
      for (Size current_peak = 0; current_peak < bestNumPeaks; current_peak++)
      {

        // Store the current parameters for this peak

        peaks[current_peak].left_width  = bestResult(0);
        peaks[current_peak].right_width = bestResult(1);

        peaks[current_peak].height = bestResult(2 + 2 * current_peak);
        peaks[current_peak].mz_position = bestResult(2 + 2 * current_peak + 1);



        // compute the area
        // is it a Lorentz or a Sech - Peak?
        if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
        {
          PeakShape p = peaks[current_peak];
          double x_left_endpoint = p.mz_position + 1 / p.left_width * sqrt(p.height / 1 - 1);
          double x_right_endpoint = p.mz_position + 1 / p.right_width * sqrt(p.height / 1 - 1);
#ifdef DEBUG_DECONV
          std::cout << "x_left_endpoint " << x_left_endpoint << " x_right_endpoint " << x_right_endpoint << std::endl;
          std::cout << "p.height" << p.height << std::endl;
#endif
          double area_left = -p.height / p.left_width * atan(p.left_width * (x_left_endpoint - p.mz_position));
          double area_right = -p.height / p.right_width * atan(p.right_width * (p.mz_position - x_right_endpoint));
          peaks[current_peak].area = area_left + area_right;

        }
        else                   //It's a Sech - Peak
        {
          PeakShape p = peaks[current_peak];
          double x_left_endpoint = p.mz_position + 1 / p.left_width * boost::math::acosh(sqrt(p.height / 0.001));
          double x_right_endpoint = p.mz_position + 1 / p.right_width * boost::math::acosh(sqrt(p.height / 0.001));
#ifdef DEBUG_DECONV
          std::cout << "x_left_endpoint " << x_left_endpoint << " x_right_endpoint " << x_right_endpoint << std::endl;
          std::cout << "p.height" << p.height << std::endl;
#endif
          double area_left = -p.height / p.left_width * (sinh(p.left_width * (p.mz_position - x_left_endpoint))
                                                         / cosh(p.left_width * (p.mz_position - x_left_endpoint)));
          double area_right = -p.height / p.right_width * (sinh(p.right_width * (p.mz_position - x_right_endpoint))
                                                           / cosh(p.right_width * (p.mz_position - x_right_endpoint)));
          peaks[current_peak].area = area_left + area_right;

        }

      }
    }
    charge_ = bestCharge;

    return true;
  }

  Size OptimizePeakDeconvolution::getNumberOfPeaks_(Int charge, std::vector<PeakShape>& temp_shapes, Data& data)
  {
    double dist = dist_ / charge;

    data.peaks.clear();

    Size shape = 0;
#ifdef DEBUG_DECONV
    std::cout << "temp_shapes[0].mz_position " << temp_shapes[0].mz_position
              << "\t dist " << dist << "\tp_index " << shape << std::endl;
#endif
    // while the peak's position is smaller than the last considered position
    // take the peak for optimization
    while ((temp_shapes[0].mz_position + shape * dist <
            data.positions[data.positions.size() - 1]) &&
           (shape < temp_shapes.size()))
    {
      data.peaks.push_back(temp_shapes[shape]);
#ifdef DEBUG_DECONV
      std::cout << "temp_shapes[0].mz_position + p_index*dist = " << temp_shapes[0].mz_position + shape * dist << std::endl;
#endif
      ++shape;
    }

    return shape;

  }

  void OptimizePeakDeconvolution::setNumberOfPeaks_(Data& data, const std::vector<PeakShape>& temp_shapes, Int charge)
  {
    double dist = dist_ / charge;

    data.peaks.clear();
#ifdef DEBUG_DECONV
    std::cout << "temp_shapes[0].mz_position " << temp_shapes[0].mz_position
              << "\t dist " << dist << "\tp_index " << shape << std::endl;
#endif
    // while the peak's position is smaller than the last considered position
    // take the peak for optimization
    Size shape = 0;
    while ((temp_shapes[0].mz_position + shape * dist <
            data.positions[data.positions.size() - 1]) &&
           (shape < temp_shapes.size()))
    {
      data.peaks.push_back(temp_shapes[shape]);
#ifdef DEBUG_DECONV
      std::cout << "temp_shapes[0].mz_position + p_index*dist = " << temp_shapes[0].mz_position + shape * dist << std::endl;
#endif
      shape++;
    }
  }

}
