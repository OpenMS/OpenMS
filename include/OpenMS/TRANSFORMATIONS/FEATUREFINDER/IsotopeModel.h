// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

namespace OpenMS
{
  class EmpiricalFormula;

  /** 
		@brief Isotope distribution approximated using linear interpolation.

    This models a smoothed (widened) distribution, i.e. can be used to sample actual raw peaks (depending on the points you query).
    If you only want the distribution (no widening), use either
    EmpiricalFormula::getIsotopeDistribution() // for a certain sum formula
    or
    IsotopeDistribution::estimateFromPeptideWeight (double average_weight)  // for averagine

    Peak widening is achieved by either a Gaussian or Lorentzian shape.

		@htmlinclude OpenMS_IsotopeModel.parameters
	*/
	class OPENMS_DLLAPI IsotopeModel
  : public InterpolationModel
  {

		public:
		typedef InterpolationModel::CoordinateType CoordinateType;
		typedef InterpolationModel::CoordinateType IntensityType;

		enum Averagines{C=0,H,N,O,S,AVERAGINE_NUM};

    /// Default constructor
    IsotopeModel();

    ///  copy constructor
  	IsotopeModel(const IsotopeModel& source);

    /// destructor
    virtual ~IsotopeModel();

    /// assignment operator
    virtual IsotopeModel& operator = (const IsotopeModel& source);

		UInt getCharge();

		/// create new IsotopeModel object (needed by Factory)
		static BaseModel<1>* create()
    {
	     return new IsotopeModel();
  	}

		/// name of the model (needed by Factory)
    static const String getProductName()
    {
	     return "IsotopeModel";
  	}

		/** @brief set the offset of the model

			The whole model will be shifted to the new offset without being computing all over.
			This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
			standard deviations) but can get significant otherwise. In that case use setParameters()
			which enforces a recomputation of the model.
		*/
		void setOffset(CoordinateType offset);

		CoordinateType getOffset();

    /// return the Averagine peptide formula (mass calculated from mean mass and charge -- use .setParameters() to set them)
    EmpiricalFormula getFormula();

		/// set sample/supporting points of interpolation
		void setSamples(const EmpiricalFormula& formula);

		/** @brief get the center of the Isotope model

			 This is a m/z-value not necessarily the monoisotopic mass.
		*/
		CoordinateType getCenter() const;

    /** @brief the Isotope distribution (without widening) from the last setSamples() call

      Useful to determine the number of isotopes that the model contains and their position

    */
    const IsotopeDistribution& getIsotopeDistribution() const;


		protected:
			CoordinateType isotope_stdev_;
      CoordinateType isotope_lorentz_fwhm_;

			UInt charge_;
			CoordinateType mean_;
			CoordinateType monoisotopic_mz_;
			DoubleReal averagine_[AVERAGINE_NUM];
			Int max_isotope_;
			DoubleReal trim_right_cutoff_;
			DoubleReal isotope_distance_;
      IsotopeDistribution isotope_distribution_;

  		void updateMembers_();

  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H
