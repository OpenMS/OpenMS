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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

namespace OpenMS
{
/** @brief Normal distribution approximated using linear interpolation

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
	<tr><td>statistics: mean, variance</td>
			<td>mean and variance of the data used to fit the model</td></tr>
	</table>
	
	@todo Remove setParam method and use setParameters instead (Ole)
	
	@ingroup FeatureFinder
*/
class GaussModel
	: public InterpolationModel<>
{

public:
    typedef InterpolationModel<>::CoordinateType CoordinateType;
    typedef Math::BasicStatistics<CoordinateType > BasicStatistics;
    typedef InterpolationModel<> InterpolationModel;

    /// Default constructor
    GaussModel();

    /// copy constructor
    GaussModel(const GaussModel& source);

    /// destructor
    virtual ~GaussModel();

    /// assignment operator
    virtual GaussModel& operator = (const GaussModel& source);
    
    void setParam(const BasicStatistics& statistics, CoordinateType min, CoordinateType max);

    /// create new GaussModel object (needed by Factory)
    static BaseModel<1>* create()
    {
        return new GaussModel();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
        return "GaussModel";
    }

    /** @brief set the offset of the model

    	The whole model will be shifted to the new offset without being computing all over.
    	and without any discrepancy.
    */
    void setOffset(CoordinateType  offset);

    /// set sample/supporting points of interpolation
    void setSamples();

    /// get the center of the Gaussian model i.e. the position of the maximum
    const CoordinateType getCenter() const;

protected:
    CoordinateType  min_;
    CoordinateType  max_;
    BasicStatistics statistics_;

		void updateMembers_();
};
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSMODEL_H
