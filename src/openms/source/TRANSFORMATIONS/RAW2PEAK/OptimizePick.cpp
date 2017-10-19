// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
#include <algorithm>
#include <cmath>

#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>


using std::max;


namespace OpenMS
{

  OptimizePick::OptimizePick(
      const struct OptimizationFunctions::PenaltyFactors & penalties,
      const int max_iteration)
  {

    penalties_ = penalties;

    max_iteration_ = max_iteration;

#ifdef DEBUG_PEAK_PICKING
    std::cout << "max iteration " << max_iteration_
              << "\n penalty factor pos " << penalties.pos
              << "\n penalty factor left width " << penalties.lWidth
              << "\n penalty factor right width " << penalties.rWidth
              << std::endl;
#endif

  }

  OptimizePick::~OptimizePick()
  {
  }

  void OptimizePick::optimize(std::vector<PeakShape> & peaks, Data & data)
  {
    if (peaks.empty())
      return;

    size_t global_peak_number = 0;
    data.peaks.assign(peaks.begin(), peaks.end());

    size_t num_dimensions = 4 * data.peaks.size();
    Eigen::VectorXd x_init (num_dimensions);
    x_init.setZero();
    // We have to initialize the parameters for the optimization
    for (size_t i = 0; i < data.peaks.size(); i++)
    {
      PeakShape current_peak = data.peaks[i];
      double h  = current_peak.height;
      double wl = current_peak.left_width;
      double wr = current_peak.right_width;
      double p  = current_peak.mz_position;
      if (boost::math::isnan(wl))
      {
        data.peaks[i].left_width = 1;
        wl = 1.;
      }
      if (boost::math::isnan(wr))
      {
        data.peaks[i].right_width = 1;
        wr = 1.;
      }
      x_init(4 * i) = h;
      x_init(4 * i + 1) = wl;
      x_init(4 * i + 2) = wr;
      x_init(4 * i + 3) = p;
    }

    data.penalties = penalties_;

    unsigned num_data_points = std::max(data.positions.size() + 1, num_dimensions);
    OptPeakFunctor functor (num_dimensions, num_data_points, &data);
    Eigen::LevenbergMarquardt<OptPeakFunctor> lmSolver (functor);
    lmSolver.parameters.maxfev = max_iteration_;
    Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);
    //the states are poorly documented. after checking the source, we believe that
    //all states except NotStarted, Running and ImproperInputParameters are good
    //termination states.
    if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
    {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-OptimizePeak:", "Could not fit the data: Error " + String(status));
    }

    // iterate over all peaks and store the optimized values in peaks
    for (size_t current_peak = 0; current_peak < data.peaks.size(); current_peak++)
    {
      // Store the current parameters for this peak
      peaks[global_peak_number + current_peak].height =  x_init(4 * current_peak);
      peaks[global_peak_number + current_peak].mz_position = x_init(4 * current_peak + 3);
      peaks[global_peak_number + current_peak].left_width = x_init(4 * current_peak + 1);
      peaks[global_peak_number + current_peak].right_width = x_init(4 * current_peak + 2);

      // compute the area
      // is it a Lorentz or a Sech - Peak?
      if (peaks[global_peak_number + current_peak].type == PeakShape::LORENTZ_PEAK)
      {
        PeakShape p = peaks[global_peak_number + current_peak];
        double x_left_endpoint = p.mz_position - 1 / p.left_width * sqrt(p.height / 1 - 1);
        double x_right_endpoint = p.mz_position + 1 / p.right_width * sqrt(p.height / 1 - 1);
        double area_left = -p.height / p.left_width * atan(p.left_width * (x_left_endpoint - p.mz_position));
        double area_right = -p.height / p.right_width * atan(p.right_width * (p.mz_position - x_right_endpoint));
        peaks[global_peak_number + current_peak].area = area_left + area_right;
#ifdef DEBUG_PEAK_PICKING
        std::cout << "Lorentz " << area_left << " " << area_right
                  << " " << peaks[global_peak_number + current_peak].area << std::endl;
#endif
      }
      else  //It's a Sech - Peak
      {
        PeakShape p = peaks[global_peak_number + current_peak];
        double x_left_endpoint = p.mz_position - 1 / p.left_width * boost::math::acosh(sqrt(p.height / 0.001));
        double x_right_endpoint = p.mz_position + 1 / p.right_width * boost::math::acosh(sqrt(p.height / 0.001));
        double area_left = p.height / p.left_width * (sinh(p.left_width * (p.mz_position - x_left_endpoint)) / cosh(p.left_width * (p.mz_position - x_left_endpoint)));
        double area_right = -p.height / p.right_width * (sinh(p.right_width * (p.mz_position - x_right_endpoint)) / cosh(p.right_width * (p.mz_position - x_right_endpoint)));
        peaks[global_peak_number + current_peak].area = area_left + area_right;
#ifdef DEBUG_PEAK_PICKING
        std::cout << "Sech " << area_left << " " << area_right
                  << " " << peaks[global_peak_number + current_peak].area << std::endl;
        std::cout << p.mz_position << " " << x_left_endpoint << " " << x_right_endpoint << std::endl;
#endif
      }
    }
    //global_peak_number += data.peaks.size();

  }

  int OptimizePick::OptPeakFunctor::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec)
  {
    const std::vector<double> & signal = m_data->signal;
    const std::vector<double> & positions = m_data->positions;
    const std::vector<PeakShape> & peaks = m_data->peaks;
    const OptimizationFunctions::PenaltyFactors & penalties = m_data->penalties;
    // iterate over all points of the signal
    for (size_t current_point = 0; current_point < positions.size(); current_point++)
    {
      double computed_signal = 0.;
      double current_position = positions[current_point];
      double experimental_signal = signal[current_point];

      // iterate over all peaks
      for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
      {
        // Store the current parameters for this peak
        double p_height = x(4 * current_peak);
        double p_position = x(4 * current_peak + 3);
        double p_width
            = (current_position <= p_position) ? x(4 * current_peak + 1)
                                 : x(4 * current_peak + 2);

        // is it a Lorentz or a Sech - Peak?
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

    double penalty = 0.;
    double penalty_pos    = penalties.pos;
    double penalty_lwidth = penalties.lWidth;
    double penalty_rwidth = penalties.rWidth;

    // iterate over all peaks again to compute the penalties
    for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
    {
      double old_position = peaks[current_peak].mz_position;
      double old_width_l = peaks[current_peak].left_width;
      double old_width_r = peaks[current_peak].right_width;
      double p_position = x(4 * current_peak + 3);
      double p_width_l = x(4 * current_peak + 1);
      double p_width_r = x(4 * current_peak + 2);

      //penalty += pow(p_position - old_position, 2) + pow(p_width_l - old_width_l, 2) + pow(p_width_r - old_width_r, 2);
      penalty += penalty_pos * pow(p_position - old_position, 2)
        + penalty_lwidth * pow(p_width_l - old_width_l, 2)
        + penalty_rwidth * pow(p_width_r - old_width_r, 2);
    }

    fvec(positions.size()) = 100 * penalty;

    return 0;
  }
  // compute Jacobian matrix for the different parameters
  int OptimizePick::OptPeakFunctor::df(const Eigen::VectorXd &x, Eigen::MatrixXd &J)
  {
    std::cout << "rows: " << J.rows() << " colums: " << J.cols() << std::endl;//DEBUG
    const std::vector<double> & positions = m_data->positions;
    const std::vector<PeakShape> & peaks = m_data->peaks;
    const OptimizationFunctions::PenaltyFactors & penalties = m_data->penalties;
    // iterate over all points of the signal
    for (size_t current_point = 0; current_point < positions.size(); current_point++)
    {
      double current_position = positions[current_point];

      // iterate over all peaks
      for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
      {
        // Store the current parameters for this peak
        double p_height = x(4 * current_peak);
        double p_position = x(4 * current_peak + 3);
        double p_width = (current_position <= p_position) ? x(4 * current_peak + 1)
                          : x(4 * current_peak + 2);

        // is it a Lorentz or a Sech - Peak?
        if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
        {
          double diff = current_position - p_position;
          double denom_inv = 1. / (1. + pow(p_width * diff, 2));

          double ddl_left = (current_position <= p_position)
              ? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2) : 0;

          double ddl_right = (current_position  > p_position)
              ? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2) : 0;

          double ddx0 = -2 * p_height * pow(p_width, 2) * diff * pow(denom_inv, 2);

          J(current_point, 4 * current_peak) = denom_inv;
          J(current_point, 4 * current_peak + 1) = ddl_left;
          J(current_point, 4 * current_peak + 2) = ddl_right;
          J(current_point, 4 * current_peak + 3) = ddx0;
        }
        else // It's a Sech - Peak
        {
          double diff = current_position - p_position;
          double denom_inv = 1. / cosh(p_width * diff);

          // The remaining computations are not stable if denom_inv == 0. In that case, we are far away from the peak
          // and can assume that all derivatives vanish
          double sinh_term = (fabs(denom_inv) < 1e-6) ? 0.0 : sinh(p_width * diff);
          double ddl_left  = (current_position <= p_position)
              ? -2 * p_height * sinh_term * diff * pow(denom_inv, 3) : 0;
          double ddl_right = (current_position  > p_position)
              ? -2 * p_height * sinh_term * diff * pow(denom_inv, 3) : 0;
          double ddx0      = 2 * p_height * p_width * sinh_term * pow(denom_inv, 3);

          J(current_point, 4 * current_peak) = pow(denom_inv, 2);
          J(current_point, 4 * current_peak + 1) = ddl_left;
          J(current_point, 4 * current_peak + 2) = ddl_right;
          J(current_point, 4 * current_peak + 3) = ddx0;
        }
      }
    }

    // Now iterate over all peaks again to compute the penalties.
    for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
    {
      double p_width_left = x(4 * current_peak + 1);
      double p_width_right = x(4 * current_peak + 2);
      double p_position = x(4 * current_peak + 3);

      double old_width_left = peaks[current_peak].left_width;
      double old_width_right = peaks[current_peak].right_width;
      double old_position = peaks[current_peak].mz_position;


      double penalty_l = 2. * penalties.lWidth * (p_width_left - old_width_left);
      double penalty_r = 2. * penalties.rWidth * (p_width_right - old_width_right);
      double penalty_p = 0;
      if (fabs(p_position - old_position) < 0.2)
      {
        penalty_p = 2. * penalties.pos * (p_position - old_position);
      }

      J(positions.size(), 4 * current_peak) = 0.;
      J(positions.size(), 4 * current_peak + 1) = 100 * penalty_l;
      J(positions.size(), 4 * current_peak + 2) = 100 * penalty_r;
      J(positions.size(), 4 * current_peak + 3) = 100 * penalty_p;
    }
    return 0;
  }
}//namespace
