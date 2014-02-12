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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_OPTIMIZEPICK_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_OPTIMIZEPICK_H

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <fstream>
#include <vector>

namespace OpenMS
{
  /** @brief Namespace for all functions and classes needed for the gsl levenberg-marquard algorithm.

      We have to use function pointers for the gsl and can't put them into
      a class, so we provide an extra namespace.
  */
  namespace OptimizationFunctions
  {
    /// Raw data vector type
    typedef std::vector<Peak1D> RawDataVector;
    /// Raw data iterator type
    typedef RawDataVector::iterator PeakIterator;

    /** @brief Class for the penalty factors used during the optimization.

        A great deviation (squared deviation) of a peak shape's position or its left or right width parameter can be penalised.
        In each iteration the penalty (for each peak shape) is computed by:
                penalty = penalty_pos * pow(p_position - old_position, 2)
                        + penalty_lwidth * pow(p_width_l - old_width_l, 2)
                        + penalty_rwidth * pow(p_width_r - old_width_r, 2);
    */
    struct OPENMS_DLLAPI PenaltyFactors
    {
      PenaltyFactors() :
        pos(0), lWidth(0), rWidth(0) {}
      PenaltyFactors(const PenaltyFactors & p) :
        pos(p.pos), lWidth(p.lWidth), rWidth(p.rWidth) {}
      inline PenaltyFactors & operator=(const PenaltyFactors & p)
      {
        pos = p.pos;
        lWidth = p.lWidth;
        rWidth = p.rWidth;

        return *this;
      }

      ~PenaltyFactors(){}

      /// Penalty factor for the peak shape's position
      double pos;
      /// Penalty factor for the peak shape's left width parameter
      double lWidth;
      /// Penalty factor for the peak shape's right width parameter
      double rWidth;
    };

    /// Evaluation of the target function for nonlinear optimization.
    int residual(const gsl_vector * x, void * params, gsl_vector * f);

    /// Compute the Jacobian of the residual, where each row of the matrix corresponds to a point in the data.
    int jacobian(const gsl_vector * x, void * params, gsl_matrix * J);

    /// Driver function for the evaluation of function and jacobian.
    int evaluate(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J);

    /// Print all peak shapes
    void printSignal(const gsl_vector * x, void * param, float resolution = 0.25);
  }


  /**
    @brief This class provides the non-linear optimization of the peak parameters.

    Given a vector of peak shapes, this class optimizes all peak shapes parameters using a non-linear optimization.
    For the non-linear optimization we use the Levenberg-Marquardt algorithm provided by the gsl.
  */
  class OPENMS_DLLAPI OptimizePick
  {
public:

    struct Data
    {
      /// Positions and intensity values of the raw data
      std::vector<double> positions;
      std::vector<double> signal;
      /// This container contains the peak shapes to be optimized
      std::vector<PeakShape> peaks;

      OptimizationFunctions::PenaltyFactors penalties;

    };


    /// Raw data vector type
    typedef std::vector<Peak1D> RawDataVector;
    /// Raw data iterator type
    typedef RawDataVector::iterator PeakIterator;


    /// Constructor
    OptimizePick() :
      max_iteration_(0),
      eps_abs_(0),
      eps_rel_(0) {}

    /// Constructor to set the penalty factors, the number of optimization iterations as well as the threshold for the absolute and the relative error.
    OptimizePick(const struct OptimizationFunctions::PenaltyFactors & penalties_,
                 const int max_iteration_,
                 const double eps_abs_,
                 const double eps_rel_);

    /// Destructor
    ~OptimizePick();

    /// Non-mutable access to the penalty factors
    inline const struct OptimizationFunctions::PenaltyFactors & getPenalties() const { return penalties_; }
    /// Mutable access to the penalty factors
    inline struct OptimizationFunctions::PenaltyFactors & getPenalties() { return penalties_; }
    /// Mutable access to the penalty factors
    inline void setPenalties(const struct OptimizationFunctions::PenaltyFactors & penalties) { penalties_ = penalties; }

    /// Non-mutable access to the number of iterations
    inline UInt getNumberIterations() const { return max_iteration_; }
    /// Mutable access to the number of iterations
    inline unsigned int & getNumberIterations() { return max_iteration_; }
    /// Mutable access to the number of iterations
    inline void setNumberIterations(const int max_iteration) { max_iteration_ = max_iteration; }

    /// Non-mutable access to the maximum absolute error
    inline DoubleReal getMaxAbsError() const { return eps_abs_; }
    /// Mutable access to the maximum absolute error
    inline double & getMaxAbsError() { return eps_abs_; }
    /// Mutable access to the maximum absolute error
    inline void setMaxAbsError(double eps_abs) { eps_abs_ = eps_abs; }

    /// Non-mutable access to the maximum relative error
    inline DoubleReal getMaxRelError() const { return eps_rel_; }
    /// Mutable access to the maximum relative error
    inline double & getMaxRelError() { return eps_rel_; }
    /// Mutable access to the maximum relative error
    inline void setMaxRelError(double eps_rel) { eps_rel_ = eps_rel; }

    /// Start the optimization of the peak shapes peaks. The original peak shapes will be substituted by the optimized peak shapes.
    void optimize(std::vector<PeakShape> & peaks, Data & data);


protected:
    /// Penalty factors
    struct OptimizationFunctions::PenaltyFactors penalties_;

    /// Maximum number of iterations during optimization
    unsigned int max_iteration_;

    /// Maximum absolute and relative error used in the optimization.
    double eps_abs_;
    double eps_rel_;

//       /** @brief Returns the squared pearson coefficient.

//         Computes the correlation of the peak and the original data given by the peak endpoints.
//         If the value is near 1, the fitted peakshape and the raw data are expected to be very similar.
//     */
//     double correlate_(const PeakShape& peak,
//                                          double left_endpoint,
//                                          double right_endpoint,Data& data);

  };
}

#endif
