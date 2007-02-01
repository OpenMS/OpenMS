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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>

namespace OpenMS
{
  /** @brief Isotope distribution approximated using linear interpolation.

		Parameters:
		<table>
			<tr><td>interpolation_step</td>
					<td>step size used to interpolate model</td></tr>
			<tr><td>intensity_scaling</td>
					<td>factor used to scale the calculated intensities</td></tr>
			<tr><td>cutoff</td>
					<td>peak with intensity below cutoff is not considered
							 to be part of the model</td></tr>
			<tr><td>charge</td>
					<td>charge of the isotope distribution.</td></tr>
			<tr><td>isotope:stdev</td>
					<td>standard deviation of the isotope distribution,
							used to account for different data resolutions.</td></tr>
			<tr><td>isotope:distance</td>
					<td>distance between two isotopes of charge +1</td></tr>
			<tr><td>isotope:trim_right_cutoff</td>
					<td>use only isotopes with abundancies above this cutoff</td></tr>
			<tr><td>isotope:maximum</td>
					<td>maximum number of isotopes being used for the IsotopeModel</td></tr>
			<tr><td>statistics:mean</td>
					<td>mean of the data used to fit the model.</td></tr>
			<tr><td>avergines:C, N, H, O, S</td>
					<td>averagines are used to approximate the number of atoms of a given element
					 (C,H,N,O,S) given a mass</td></tr>
		</table>
		
		@todo Remove setParam method and use setParameters instead (Ole)
		
		@ingroup FeatureFinder
	*/
	class IsotopeModel
  : public InterpolationModel<>
  {

		public:
		typedef InterpolationModel<>::CoordinateType CoordinateType;
		typedef InterpolationModel<>::CoordinateType IntensityType;

		enum Averagines{C=0,H,N,O,S,AVERAGINE_NUM};

    /// Default constructor
    IsotopeModel();

    ///  copy constructor
  	IsotopeModel(const IsotopeModel& source);

    /// destructor 
    virtual ~IsotopeModel();

    /// assignment operator 
    virtual IsotopeModel& operator = (const IsotopeModel& source);

		void setParam(CoordinateType mean, UnsignedInt charge, CoordinateType isotope_stdev);

		UnsignedInt getCharge();

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
			standard deviations) but can get significant otherwise. In that case use setParam()
			which enforces a recomputation of the model.
		*/
		void setOffset(CoordinateType offset);

		const CoordinateType& getOffset();

		/// set sample/supporting points of interpolation
		void setSamples();

		/** @brief get the center of the Isotope model i.e. the position of the monoisotopic peak

			 This is a m/z-value not necessarily the monoisotopic mass.
		*/
		const CoordinateType getCenter() const;

		protected:
			CoordinateType isotope_stdev_;
			UnsignedInt charge_;
			CoordinateType mean_;
			CoordinateType monoisotopic_mz_;
			double averagine_[AVERAGINE_NUM];
			int max_isotope_;
			double trim_right_cutoff_;
			double isotope_distance_;
  
  		void updateMembers_();
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEMODEL_H
