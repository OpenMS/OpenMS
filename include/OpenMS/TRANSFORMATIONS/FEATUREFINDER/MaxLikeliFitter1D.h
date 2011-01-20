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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MAXLIKELIFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MAXLIKELIFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>

namespace OpenMS
{

	/**
	@brief Abstract base class for all 1D-model fitters using maximum likelihood optimization.
	*/
	class OPENMS_DLLAPI MaxLikeliFitter1D
		: public Fitter1D
	{

	 public:

		/// default constructor
		MaxLikeliFitter1D()
			: Fitter1D()
		{
		}

		/// copy constructor
		MaxLikeliFitter1D(const MaxLikeliFitter1D& source)
			: Fitter1D(source)
		{
		}

		/// destructor
		virtual ~MaxLikeliFitter1D()
		{
		}

		/// assignment operator
		virtual MaxLikeliFitter1D& operator = (const MaxLikeliFitter1D& source)
		{
			if (&source ==this) return *this;

			Fitter1D::operator = (source);

			return *this;
		}

	 protected:

		/// fit an offset on the basis of the pearson correlation coefficient
		QualityType fitOffset_(InterpolationModel* model, const RawDataArrayType& set, const CoordinateType stdev1, const CoordinateType stdev2, const CoordinateType offset_step)
		{
			const CoordinateType offset_min = model->getInterpolation().supportMin() - stdev1;
			const CoordinateType offset_max = model->getInterpolation().supportMin() + stdev2;

			CoordinateType offset;
			QualityType correlation;

			//test model with default offset
			std::vector<Real> real_data;
			real_data.reserve(set.size());
			std::vector<Real> model_data;
			model_data.reserve(set.size());

			for (Size i=0; i < set.size(); ++i)
			{
				real_data.push_back(set[i].getIntensity());
				model_data.push_back( model->getIntensity( DPosition<1>(set[i].getPosition()) ) );
			}

			CoordinateType max_offset = model->getInterpolation().getOffset();
			QualityType max_correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());

			//test different offsets
			for ( offset = offset_min; offset <= offset_max; offset += offset_step )
			{
				// set offset
				model->setOffset( offset );

				// get samples
				model_data.clear();
				for (Size i=0; i < set.size(); ++i)
				{
					model_data.push_back( model->getIntensity( DPosition<1>(set[i].getPosition()) ) );
				}

				correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());

				if ( correlation > max_correlation )
				{
					max_correlation = correlation;
					max_offset = offset;
				}
			}

			model->setOffset( max_offset );

			return max_correlation;
		}

		void updateMembers_()
		{
			Fitter1D::updateMembers_();
		}

	};
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_MAXLIKELIFITTER1D_H
