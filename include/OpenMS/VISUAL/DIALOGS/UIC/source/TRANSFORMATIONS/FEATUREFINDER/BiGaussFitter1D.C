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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <boost/math/special_functions/fpclassify.hpp>
namespace OpenMS
{
    BiGaussFitter1D::BiGaussFitter1D()			
    : MaxLikeliFitter1D()
    {
        setName(getProductName());
    
        defaults_.setValue("statistics:variance1",1.0,"Variance of the first gaussian, used for the lower half of the model.", StringList::create("advanced"));
        defaults_.setValue("statistics:variance2",1.0,"Variance of the second gaussian, used for the upper half of the model.", StringList::create("advanced"));
    
        defaultsToParam_();
    }
		
    BiGaussFitter1D::BiGaussFitter1D(const BiGaussFitter1D& source)
    : MaxLikeliFitter1D(source)
    {
        updateMembers_();
    }

    BiGaussFitter1D::~BiGaussFitter1D()
    {
    }

    BiGaussFitter1D& BiGaussFitter1D::operator = (const BiGaussFitter1D& source)
    {
        if (&source == this) return *this;
    
        MaxLikeliFitter1D::operator = (source);
        updateMembers_();
    
        return *this;
    }
		
    BiGaussFitter1D::QualityType BiGaussFitter1D::fit1d(const RawDataArrayType& set, InterpolationModel*& model)
    {
        // Calculate bounding box
        min_ = max_ = set[0].getPos();
        for ( UInt pos=1; pos < set.size(); ++pos)
        {
            CoordinateType tmp = set[pos].getPos();
            if ( min_ > tmp ) min_ = tmp;
            if ( max_ < tmp ) max_ = tmp;
        }
    
        // Enlarge the bounding box by a few multiples of the standard deviation
        {
            stdev1_ = sqrt ( statistics1_.variance() ) * tolerance_stdev_box_;
            stdev2_ = sqrt ( statistics2_.variance() ) * tolerance_stdev_box_;
            min_ -= stdev1_;
            max_ += stdev2_;
        }
        
        // build model
        model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("BiGaussModel"));
        model->setInterpolationStep( interpolation_step_ );
        Param tmp;
        tmp.setValue( "bounding_box:min", min_ );
        tmp.setValue( "bounding_box:max", max_ );
        tmp.setValue( "statistics:mean", statistics1_.mean() );
        tmp.setValue( "statistics:variance1", statistics1_.variance() );
        tmp.setValue( "statistics:variance2", statistics2_.variance() );
        model->setParameters( tmp );
        
        // fit offset
        QualityType quality;
        quality = fitOffset_(model, set, stdev1_, stdev2_, interpolation_step_);
        if (boost::math::isnan(quality) ) quality = -1.0;
        
        return quality;
    }
		
    void BiGaussFitter1D::updateMembers_()
    {
        MaxLikeliFitter1D::updateMembers_();
        statistics1_.setMean(param_.getValue("statistics:mean"));
        statistics1_.setVariance(param_.getValue("statistics:variance1"));
        statistics2_.setMean(param_.getValue("statistics:mean"));
        statistics2_.setVariance(param_.getValue("statistics:variance2"));
    }

}
