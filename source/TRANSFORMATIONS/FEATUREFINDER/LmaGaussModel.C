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


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussModel.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <numeric>

namespace OpenMS
{
    LmaGaussModel::LmaGaussModel()
      : InterpolationModel()
    {
      setName(getProductName());
      
      defaults_.setValue("bounding_box:min",0.0f,"Lower end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
      defaults_.setValue("bounding_box:max",1.0f,"Upper end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:mean",0.0f,"Centroid position of the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:variance",1.0f,"The variance of the model.", StringList::create("advanced"));
      defaults_.setValue("lma:scale_factor",1000000.0f,"Scale factor for the intensity of the model.", StringList::create("advanced"));
      defaults_.setValue("lma:standard_deviation",5.0f,"The standard deviation (variance) of the model.", StringList::create("advanced"));
      defaults_.setValue("lma:expected_value",1200.0f,"The expected value (centroid position) of the model.", StringList::create("advanced"));
  
      defaultsToParam_();
    }
  
    LmaGaussModel::LmaGaussModel(const LmaGaussModel& source)
    : InterpolationModel(source)
    {
      setParameters( source.getParameters() );
      updateMembers_();
    }
  
    LmaGaussModel::~LmaGaussModel()
    {
    }
  
    LmaGaussModel& LmaGaussModel::operator = (const LmaGaussModel& source)
    {
      if (&source == this) return *this;
      
      InterpolationModel::operator = (source);
      setParameters( source.getParameters() );
      updateMembers_();
      
      return *this;
    }
  
    void LmaGaussModel::setSamples()
    {
      LinearInterpolation::container_type& data = interpolation_.getData();
      data.clear();
      if (max_==min_) return;
      data.reserve( UInt ( (max_ - min_) / interpolation_step_ + 1 ) );
      CoordinateType pos = min_;
  
      DoubleReal part1 = 1/(sqrt(2*Constants::PI)*standard_deviation_);
      DoubleReal part2 = (2*pow(standard_deviation_,2));
  
      for ( UInt i = 0; pos< max_; ++i)
      {
        pos = min_ + i * interpolation_step_;
        DoubleReal tmp = pos - expected_value_;
  
        // data.push_back(Gauss)
        data.push_back((part1*exp(-(pow(tmp,2))/part2)*scale_factor_));
      }
  
      interpolation_.setScale  ( interpolation_step_ );
      interpolation_.setOffset ( min_ );
    }
  
    void LmaGaussModel::setOffset(CoordinateType offset)
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
  
    LmaGaussModel::CoordinateType LmaGaussModel::getCenter() const
    {
      return statistics_.mean();
    }
  
    void LmaGaussModel::updateMembers_()
    {
      InterpolationModel::updateMembers_();
      
      min_ = param_.getValue("bounding_box:min");
      max_ = param_.getValue("bounding_box:max");
      statistics_.setMean( param_.getValue("statistics:mean") );
      statistics_.setVariance(param_.getValue("statistics:variance"));
      scale_factor_= param_.getValue("lma:scale_factor");
      standard_deviation_ = param_.getValue("lma:standard_deviation");
      expected_value_= param_.getValue("lma:expected_value");
      
      setSamples();
    }

}
