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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <numeric>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
    EmgModel::EmgModel()
		: InterpolationModel()
    {
      setName(getProductName());
      
      defaults_.setValue("bounding_box:min",0.0f,"Lower end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
      defaults_.setValue("bounding_box:max",1.0f,"Upper end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:mean",0.0f,"Centroid position of the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:variance",1.0f,"The variance of the model.", StringList::create("advanced"));
      defaults_.setValue("emg:height",100000.0f,"Height of the exponentially modified Gaussian.", StringList::create("advanced"));
      defaults_.setValue("emg:width",5.0f,"Width of the exponentially modified Gaussian.", StringList::create("advanced"));
      defaults_.setValue("emg:symmetry",5.0f,"Symmetry of the exponentially modified Gaussian.", StringList::create("advanced"));
      defaults_.setValue("emg:retention",1200.0f,"Retention time of the exponentially modified Gaussian.", StringList::create("advanced"));
      
      defaultsToParam_();
    }

    EmgModel::EmgModel(const EmgModel& source)
    : InterpolationModel(source)
    {
      setParameters( source.getParameters() );
      updateMembers_();
    }
  
    EmgModel::~EmgModel()
    {	
    }
  
    EmgModel& EmgModel::operator = (const EmgModel& source)
    {
      if (&source == this) return *this;
      
      InterpolationModel::operator = (source);
      setParameters( source.getParameters() );
      updateMembers_();
      
      return *this;
    }

    void EmgModel::setSamples()
    {
      LinearInterpolation::container_type& data = interpolation_.getData();
      data.clear();
      if (max_==min_) return;
      data.reserve( UInt ( (max_-min_) / interpolation_step_ + 1 ) );
      CoordinateType pos = min_;
  
      DoubleReal sqrt_2pi = sqrt(2*Constants::PI);
      DoubleReal term_sq2 = (-2.4055/sqrt(2.0));
      DoubleReal part1    = (height_*width_/symmetry_);
      DoubleReal part2    = pow(width_,2)/(2*pow(symmetry_,2));
      DoubleReal part3    = width_/symmetry_;
  
      for ( UInt i = 0; pos< max_; ++i)
      {
        pos = min_ + i * interpolation_step_;
        DoubleReal tmp = pos - retention_;
        
        // data.push_back (Simplified EMG)
        data.push_back((part1*sqrt_2pi*exp(part2-(tmp/symmetry_))/(1+exp(term_sq2*((tmp/width_)-part3)))));
      }
  
      interpolation_.setScale  ( interpolation_step_ );
      interpolation_.setOffset ( min_ );
    }

    void EmgModel::setOffset(CoordinateType offset)
    {
      DoubleReal diff = offset - getInterpolation().getOffset();
      min_ += diff;
      max_ += diff;
      statistics_.setMean(statistics_.mean() + diff);
  
      InterpolationModel::setOffset(offset);
  
      param_.setValue("bounding_box:min", min_);
      param_.setValue("bounding_box:max", max_);
      param_.setValue("statistics:mean", statistics_.mean());
    }

    EmgModel::CoordinateType EmgModel::getCenter() const
    {
      return statistics_.mean();
    }
    
    void EmgModel::updateMembers_()
    {
      InterpolationModel::updateMembers_();
  
      min_ = param_.getValue("bounding_box:min");
      max_ = param_.getValue("bounding_box:max");
      statistics_.setMean( param_.getValue("statistics:mean") );
      statistics_.setVariance(param_.getValue("statistics:variance"));
      height_ = param_.getValue("emg:height");
      width_ = param_.getValue("emg:width");
      symmetry_ = param_.getValue("emg:symmetry");
      retention_ = param_.getValue("emg:retention");
      
      setSamples();
    }
	
}
