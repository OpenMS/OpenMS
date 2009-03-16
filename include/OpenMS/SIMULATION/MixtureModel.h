// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------
#ifndef OPENMS_SIMULATION_MIXTUREMODEL_H
#define OPENMS_SIMULATION_MIXTUREMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

#include <numeric>

namespace OpenMS
{

  /** @brief A mixture model consisting of two Gaussian distributions.

		Can be used to model imperfect elution profiles.

		Parameters:
		<table>
			<tr><td>interpolation_step</td>
					<td>step size used to interpolate model</td></tr>
			<tr><td>intensity_scaling</td>
					<td>factor used to scale the calculated intensities</td></tr>
			<tr><td>cutoff</td>
					<td>peak with intensity below cutoff is not considered
							 to be part of the model</td></tr>
			<tr><td>bounding_box: min, max</td>
					<td>minimum and maximum coordinate value of bounding box enclosing
							 the data used to fit the model</td></tr>
			<tr><td>statistics: mean1, mean2, variance1, variance2</td>
					<td>means and variances  of the mixture</td></tr>
		</table>
	*/
	class OPENMS_DLLAPI MixtureModel
  : public InterpolationModel
  {
		public:
		typedef InterpolationModel::CoordinateType CoordinateType;
    typedef LinearInterpolation::container_type ContainerType;

    /// Default constructor
    MixtureModel();

    /// copy constructor
  	MixtureModel(const MixtureModel& source);

    /// destructor
    virtual ~MixtureModel();

    /// assignment operator
    virtual MixtureModel& operator = (const MixtureModel& source);

		/// create new MixtureModel object (function needed by Factory)
		static BaseModel<1>* create()
    {
	     return new MixtureModel();
  	}

		/// name of the model (needed by Factory)
    static const String getProductName()
    {
	     return "MixtureModel";
  	}

		/** @brief set the offset of the model

			The whole model will be shifted to the new offset without being computing all over.
			and without any discrepancy.
		*/
		void setOffset(double offset);

		/// set sample/supporting points of interpolation
		void setSamples();

		/// get the center of the BiGaussian model i.e. the position of the maximum
		CoordinateType getCenter() const;

		protected:

			CoordinateType min_;

			CoordinateType max_;

			/// Mixing proportion
			CoordinateType mix_prop_;

			/// First Gaussian
			Math::BasicStatistics<> statistics1_;

			/// Second Gaussian
			Math::BasicStatistics<> statistics2_;

			void updateMembers_();
 };

} // namespace OpenMS

#endif // OPENMS_SIMULATION_MIXTUREMODEL_H
