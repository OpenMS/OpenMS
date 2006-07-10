// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: OptimizePick.h,v 1.11 2006/04/11 15:29:39 elange Exp $
// $Author: elange $
// $Maintainer: Andreas Hildebrandt $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_OPTIMIZEPICK_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_OPTIMIZEPICK_H

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/KERNEL/DRawDataPoint.h>

#ifdef GSL_DEF
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#endif


#include <iostream>
#include <fstream>
#include <vector>

namespace OpenMS
{
  /**
  	@brief This class contains the non-linear optimization of peak parameters.

  	For non-linear optimization we use the Levenberg-Marquardt of the gsl.
  	We have to use function pointers for the gsl and can't put them into
  	a class, so we provide an extra namespace.

  	@ingroup PeakPicking
  	
  	@todo write a test (Eva)
   */
  namespace OptimizationFunctions
  {
    /** @name Type definitions
     */
    //@{
    typedef std::vector<DRawDataPoint<1> > RawDataVector;
    typedef RawDataVector::iterator RawDataPointIterator;
    ///
    //@}
    ///
    struct PenaltyFactors
    {
      PenaltyFactors() : pos(0), lWidth(0), rWidth(0) {}
      PenaltyFactors(const PenaltyFactors& p) : pos(p.pos), lWidth(p.lWidth), rWidth(p.rWidth) {}
      inline PenaltyFactors& operator=(const PenaltyFactors& p)
      {
        pos=p.pos;
        lWidth=p.lWidth;
        rWidth=p.rWidth;

        return *this;
      }
      ~PenaltyFactors(){}

      double pos;
      double lWidth;
      double rWidth;
    };

    /** positions and signal values **/
    extern std::vector<double> positions_;
    extern std::vector<double> signal_;
    extern std::vector<PeakShape> peaks_;

    /** Evaluation of the target function for nonlinear optimization. **/
    int residual(const gsl_vector* x, void* /* params */, gsl_vector* f);

    /** Compute the Jacobian of the residual, where each row of the matrix corresponds to a
     *  point in the data.
     */
    int jacobian(const gsl_vector* x, void* /* params */, gsl_matrix* J);

    /** Driver function for the evaluation of function and jacobian. **/
    int evaluate(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);

  }

  class OptimizePick
  {
  public:
    /** @name Type definitions
     */
    //@{
    typedef std::vector<DRawDataPoint<1> > RawDataVector;
    typedef RawDataVector::iterator RawDataPointIterator;
    //@}


    /** @name Constructors and Destructor
     */
    //@{
    OptimizePick( )
        : max_iteration_(0),
        eps_abs_(0),
    eps_rel_(0) {}
    ///
    OptimizePick(const struct OpenMS::OptimizationFunctions::PenaltyFactors& penalties_,
                 const int max_iteration_,
                 const double eps_abs_,
                 const double eps_rel_ );
    ///
    OptimizePick(const OptimizePick& opt)
        : penalties_(opt.penalties_),
        max_iteration_(opt.max_iteration_),
    eps_rel_(opt.eps_rel_){}
    ///
    ~OptimizePick();
    //@}

    /**	@name Assignment
     */
    //@{
    inline OptimizePick& operator=(const OptimizePick& opt)
    {
      penalties_=opt.penalties_;
      max_iteration_=opt.max_iteration_;
      eps_rel_=opt.eps_rel_;

      return *this;
    }
    //@}


    /**	Accessors
     */
    //@{
    /// Non-mutable access to the penalty parameter
    inline const struct OptimizationFunctions::PenaltyFactors& getPenalties() const { return penalties_; }
    /// Mutable access to the penalty parameter
    inline struct OptimizationFunctions::PenaltyFactors& getPenalties() { return penalties_; }
    /// Mutable access to the penalty parameter
    inline void setPenalties(const struct OptimizationFunctions::PenaltyFactors& penalties) { penalties_ = penalties; }

    /// Non-mutable access to the number of iterations
    inline const int getNumberIterations() const { return max_iteration_; }
    /// Mutable access to the number of iterations
    inline int getNumberIterations() { return max_iteration_; }
    /// Mutable access to the number of iterations
    inline void setNumberIterations(const int max_iteration) { max_iteration_ = max_iteration; }

    /// Non-mutable access to the maximum absolute error
    inline const double getMaxAbsError() const { return eps_abs_; }
    /// Mutable access to the maximum absolute error
    inline double getMaxAbsError() { return eps_abs_; }
    /// Mutable access to the maximum absolute error
    inline void setMaxAbsError(const double eps_abs) { eps_abs_ = eps_abs; }

    /// Non-mutable access to the maximum relative error
    inline const double getMaxRelError() const { return eps_rel_; }
    /// Mutable access to the maximum relative error
    inline double getMaxRelError() { return eps_rel_; }
    /// Mutable access to the maximum relative error
    inline void setMaxRelError(const double eps_rel) { eps_rel_ = eps_rel; }
    //@}


    /** Perform a nonlinear optimization of the peaks that have been determined for the
       *  current split array.
       */
    void optimize(std::vector<PeakShape>& peaks);

    double correlate(const PeakShape& peak,
                     double left_endpoint,
                     double right_endpoint) ;

  protected:
    // Penalty factors for some paramter in the optimization
    struct OptimizationFunctions::PenaltyFactors penalties_;

    /** Maximum number of iterations **/
    unsigned int max_iteration_;

    /** Test for the convergence of the sequence by comparing the last iteration step dx
    		with the absolute error epsabs and relative error epsrel to the current position x **/
    double eps_abs_;
    double eps_rel_;
  };
}

#endif
