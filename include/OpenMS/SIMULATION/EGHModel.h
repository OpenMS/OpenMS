// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_EGHMODEL_H
#define OPENMS_SIMULATION_EGHMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

#include <boost/math/tr1.hpp>

namespace OpenMS
{
	/**
		@brief Exponential-Gaussian hybrid distribution model for elution profiles.

    Lan K, Jorgenson JW.
    A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks.
    Journal of Chromatography A. 2001;915(1-2):1-13.
    Available at: http://linkinghub.elsevier.com/retrieve/pii/S0021967301005945

    @htmlinclude OpenMS_EGHModel.parameters

	*/
	class OPENMS_DLLAPI EGHModel
		: public InterpolationModel
	{

	 public:
		typedef InterpolationModel::CoordinateType CoordinateType;
		typedef Math::BasicStatistics<CoordinateType > BasicStatistics;
    typedef LinearInterpolation::container_type ContainerType;

		/// Default constructor
    EGHModel();

		/// copy constructor
    EGHModel(const EGHModel& source);

		/// destructor
		virtual ~EGHModel();

		/// assignment operator
		virtual EGHModel& operator = (const EGHModel& source);

		/// create new ElutionModel object (needed by Factory)
		static BaseModel<1>* create()
		{
			return new EGHModel();
  	}

		/// name of the model (needed by Factory)
		static const String getProductName()
		{
			return "EGHModel";
		}

		/// set offset without being computing all over and without any discrepancy
    void setOffset(CoordinateType offset);

		/// set sample/supporting points of interpolation
		void setSamples();

		/// get the center of the Gaussian model i.e. the position of the maximum
		CoordinateType getCenter() const;

	 protected:
		CoordinateType  min_;
		CoordinateType  max_;
		BasicStatistics statistics_;
		CoordinateType  height_; // H in paper
		CoordinateType  apex_rt_;

		CoordinateType  A_;
		CoordinateType  B_;

		CoordinateType  tau_;
		CoordinateType  sigma_square_;
		CoordinateType  sigma_square_2_;


		void updateMembers_();

    /// Computes a left & right boundary for the EGH Profile and sets the internal parameters accordingly
    void computeBoundaries_();

    /**
     * @brief Evaluate the EGH function at position rt
     *
     * @param rt        The position where the EGH function should be evaluated. Note that this is the position without the RT offset, meaning that the EGH apex is at position 0
     * @param egh_value The computed value
     */
    inline void evaluateEGH_(CoordinateType & rt, CoordinateType & egh_value)
    {
      if((sigma_square_2_ + tau_ * rt) > 0)
      {
        // evaluate egh ->
        egh_value = height_ * exp(
            (-1 * rt * rt)
            /
            (sigma_square_2_ + tau_ * rt)
            );
      }
      else
      {
        egh_value = 0.0;
      }
    }

	};

} // namespace OpenMS

#endif // OPENMS_SIMULATION_ELUTIONMODEL_H
