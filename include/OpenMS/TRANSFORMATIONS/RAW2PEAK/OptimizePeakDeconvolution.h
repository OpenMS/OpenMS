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
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

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
		*/
    struct PenaltyFactorsIntensity : public PenaltyFactors
    {
      PenaltyFactorsIntensity():PenaltyFactors(),height(0){}
      PenaltyFactorsIntensity(const PenaltyFactorsIntensity& p) : PenaltyFactors(p), height(p.height) {}
      inline PenaltyFactorsIntensity& operator=(const PenaltyFactorsIntensity& p)
      {
				height=p.height;
				pos=p.pos;
				lWidth=p.lWidth;
				rWidth=p.rWidth;
				
				return *this;
      }
      ~PenaltyFactorsIntensity(){}

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
		
		@ref OptimizePeakDeconvolution_Parameters are explained on a separate page.
		
		@ingroup PeakPicking  
	*/
  class OptimizePeakDeconvolution : public DefaultParamHandler
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
    OptimizePeakDeconvolution( );

    /// Copy-Constructor
    OptimizePeakDeconvolution(const OptimizePeakDeconvolution& opt)
      :DefaultParamHandler(opt),
			 penalties_(opt.penalties_),
			 charge_(opt.charge_){}

    ///Destructor
    virtual ~OptimizePeakDeconvolution(){}
    //@}

    /**	@name Assignment
     */
    //@{
    inline OptimizePeakDeconvolution& operator=(const OptimizePeakDeconvolution& opt)
    {
			DefaultParamHandler::operator=(opt);
      penalties_=opt.penalties_;
			charge_=opt.charge_;

      return *this;
    }
    //@}


    /**	Accessors
     */
    //@{
		/// Non-mutable access to the penalty parameter
    inline const OptimizationFunctions::PenaltyFactorsIntensity& getPenalties() const { return penalties_; }
		/// Mutable access to the penalty parameter
    inline void setPenalties(const OptimizationFunctions::PenaltyFactorsIntensity& penalties)
		{
			penalties_ = penalties;
			param_.setValue("penalties:left_width",penalties_.lWidth);
			param_.setValue("penalties:right_width",penalties_.rWidth);
			param_.setValue("penalties:height",penalties_.height);
			param_.setValue("penalties:position",penalties_.pos);
		}
		
		/// Non-mutable access to the charge
    inline const int getCharge() const { return charge_; }
    /// Mutable access to the charge
    inline void setCharge(const int charge) { charge_ = charge; }
    //@}


    /// Performs a nonlinear optimization of the peaks that belong to the current isotope pattern
    bool optimize(std::vector<PeakShape>& peaks,int failure);

  protected:
    // Penalty factors for some paramter in the optimization
    OptimizationFunctions::PenaltyFactorsIntensity penalties_;

    /// Charge state of the current isotope pattern
    int charge_;

    /// distance between two isotopic peaks
    static const double dist_;

    /// A function to determine the number of peaks that lie in the current m/z interval given the distance between the peaks by the current charge state.
    int getNumberOfPeaks_(int charge, std::vector<PeakShape>& temp_shapes);

    // After each iteration the fwhm of all peaks is checked whether it isn't too large
    bool checkFWHM_(std::vector<PeakShape>& peaks,gsl_multifit_fdfsolver *& fit);

		void updateMembers_();
  };// class
  
}// namespace OpenMS


#endif
