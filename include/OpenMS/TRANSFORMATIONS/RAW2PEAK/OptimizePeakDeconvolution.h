// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_OPTIMIZEPEAKDECONVOLUTION_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_OPTIMIZEPEAKDECONVOLUTION_H

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/FORMAT/Param.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>

//#define DEBUG_DECONV
#include <iostream>
#ifdef DEBUG_DECONV
#include <iostream>
#include <fstream>
#endif
#include <vector>

namespace OpenMS
{

  namespace OptimizationFunctions
  {
    extern std::vector<PeakShape> peaks_DC_;
    extern std::vector<double> positions_DC_;
    extern std::vector<double> signal_DC_;

		/**
			 @brief Class for the penalty factors used during the optimization.
			
			 A great deviation (squared deviation) of a peak shape's position or its left or right width parameter can be penalised.
			 During the optimization negative heights may occur, they are penalised as well.
		
			 @todo use DefaultParamHandler (Alexandra)
		*/
    struct PenaltyFactorsInt : public PenaltyFactors
    {
      PenaltyFactorsInt():PenaltyFactors(),height(0){}
      PenaltyFactorsInt(const PenaltyFactorsInt& p) : PenaltyFactors(p), height(p.height) {}
      inline PenaltyFactorsInt& operator=(const PenaltyFactorsInt& p)
      {
				height=p.height;
				pos=p.pos;
				lWidth=p.lWidth;
				rWidth=p.rWidth;
				
				return *this;
      }
      ~PenaltyFactorsInt(){}

      double height;


    };
    

    
  }//namespace OptimizationFunctions

  /**
		 @brief This class provides the deconvolution of peak regions using non-linear optimization.

		 Given a vector of peak shapes, this class optimizes all peak shapes parameters using a non-linear optimization.
		 For the non-linear optimization we use the Levenberg-Marquardt algorithm provided by the gsl.
		 There are a few constraints for the parameters: the positions are equidistant according to the peptide
		 mass rule, e.g. two consecutive isotopic peaks are 1.003/charge away from each other. Besides the
		 peaks have all the same left and right width, respectively.

		 @ingroup PeakPicking
       
	*/
  class OptimizePeakDeconvolution
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
    ///Constructor
    OptimizePeakDeconvolution( )
			: max_iteration_(0),
				eps_abs_(0),
				eps_rel_(0),
				charge_(1){}
    
    ///Constructor
    OptimizePeakDeconvolution(const OptimizationFunctions::PenaltyFactorsInt& penalties,
															const int max_iteration,
															const double eps_abs,
															const double eps_rel,
															const int charge);
    /// Copy-Constructor
    OptimizePeakDeconvolution(const OptimizePeakDeconvolution& opt)
      : penalties_(opt.penalties_),
				max_iteration_(opt.max_iteration_),
				eps_abs_(opt.eps_abs_),
				eps_rel_(opt.eps_rel_),
				charge_(opt.charge_){}

    ///Destructor
    virtual ~OptimizePeakDeconvolution(){}
    //@}

    /**	@name Assignment
     */
    //@{
    inline OptimizePeakDeconvolution& operator=(const OptimizePeakDeconvolution& opt)
    {
      penalties_=opt.penalties_;
      max_iteration_=opt.max_iteration_;
      eps_rel_=opt.eps_rel_;
      eps_abs_=opt.eps_abs_;
      charge_=opt.charge_;

      return *this;
    }
    //@}


    /**	Accessors
     */
    //@{
		/// Non-mutable access to the penalty parameter
    inline const OptimizationFunctions::PenaltyFactorsInt& getPenalties() const { return penalties_; }
    /// Mutable access to the penalty parameter
    inline OptimizationFunctions::PenaltyFactorsInt& getPenalties() { return penalties_; }
    /// Mutable access to the penalty parameter
    inline void setPenalties(const OptimizationFunctions::PenaltyFactorsInt& penalties) { penalties_ = penalties; }
    
    /// Non-mutable access to the number of iterations
    inline const int& getNumberIterations() const { return max_iteration_; }
    /// Mutable access to the number of iterations
    inline int& getNumberIterations() { return max_iteration_; }
    /// Mutable access to the number of iterations
    inline void setNumberIterations(const int max_iteration) { max_iteration_ = max_iteration; }

    /// Non-mutable access to the maximum absolute error
    inline const double& getMaxAbsError() const { return eps_abs_; }
    /// Mutable access to the maximum absolute error
    inline double& getMaxAbsError() { return eps_abs_; }
    /// Mutable access to the maximum absolute error
    inline void setMaxAbsError(const double eps_abs) { eps_abs_ = eps_abs; }

    /// Non-mutable access to the maximum relative error
    inline const double& getMaxRelError() const { return eps_rel_; }
    /// Mutable access to the maximum relative error
    inline double& getMaxRelError() { return eps_rel_; }
    /// Mutable access to the maximum relative error
    inline void setMaxRelError(const double eps_rel) { eps_rel_ = eps_rel; }

    /// Non-mutable access to the charge state
    inline const int& getCharge() const { return charge_; }
    /// Mutable access to the charge state
    inline int& getCharge() { return charge_; }
    /// Mutable access to the charge
    inline void setCharge(const int charge) { charge_ = charge; }
    //@}


    /// Performs a nonlinear optimization of the peaks that belong to the current isotope pattern
    bool optimize(std::vector<PeakShape>& peaks,Param& param,int failure);

  protected:
    // Penalty factors for some paramter in the optimization
    OptimizationFunctions::PenaltyFactorsInt penalties_;

    /// Maximum number of iterations
    int max_iteration_;

    /// Test for the convergence of the sequence by comparing the last iteration step dx with the absolute error epsabs and relative error epsrel to the current position x
    double eps_abs_;
    double eps_rel_;

    /// Charge state of the current isotope pattern
    int charge_;

    /// distance between two isotopic peaks
    static const double dist_;

    /// A function to determine the number of peaks that lie in the current m/z interval given the distance between the peaks by the current charge state.
    int getNumberOfPeaks_(int charge, std::vector<PeakShape>& temp_shapes);

    // After each iteration the fwhm of all peaks is checked whether it isn't too large
    bool checkFWHM_(std::vector<PeakShape>& peaks, Param& param,gsl_multifit_fdfsolver *& fit);
  };// class
  
}// namespace OpenMS


#endif
