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
// $Authors:  Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/MixtureModel.h>

namespace OpenMS
{

  MixtureModel::MixtureModel()
  : InterpolationModel(), mix_prop_(0.7), statistics1_(), statistics2_()
  {
    setName(getProductName());

    defaults_.setValue("bounding_box:min",0.0,"Lower end of bounding box enclosing the data used to fit the model");
    defaults_.setValue("bounding_box:max",1.0,"Upper end of bounding box enclosing the data used to fit the model");
    defaults_.setValue("statistics:mean1",0.0,"Mean of the first Gaussian.");
    defaults_.setValue("statistics:mean2",0.0,"Mean of the second Gaussian.");
    defaults_.setValue("statistics:variance1",1.0,"Variance of the first Gaussian.");
    defaults_.setValue("statistics:variance2",1.0,"Variance of the second Gaussian.");
    defaults_.setValue("mix",0.7,"Mixing proportions (less or equal to 1)");

    defaultsToParam_();
  }

  MixtureModel::MixtureModel(const MixtureModel& source)
  : InterpolationModel(source)
  {
    setParameters( source.getParameters() );
    updateMembers_();
  }

  MixtureModel::~MixtureModel()
  {
  }

  MixtureModel& MixtureModel::operator = (const MixtureModel& source)
  {
    if (&source == this) return *this;

    InterpolationModel::operator = (source);
    setParameters( source.getParameters() );
    updateMembers_();

    return *this;
  }

  void MixtureModel::setSamples()
  {
    ContainerType& data = interpolation_.getData();
    data.clear();
    if (max_==min_) return;
    data.reserve( UInt ( (max_-min_) / interpolation_step_ + 1 ) );
    CoordinateType pos = min_;

    for ( UInt i = 0; pos< max_; ++i)
    {
      pos = min_ + i * interpolation_step_;
      data.push_back( (mix_prop_ * statistics1_.normalDensity_sqrt2pi(pos)) + ( (1-mix_prop_)  * statistics2_.normalDensity_sqrt2pi(pos) ) );
    }
    // scale data so that integral over distribution equals one
    // multiply sum by interpolation_step_ -> rectangular approximation of integral
    IntensityType factor = scaling_ / interpolation_step_ / accumulate ( data.begin(), data.end(), IntensityType(0) );

    for (ContainerType::iterator it=data.begin();	it!=data.end();	++it)
    {
      *it *= factor;
    }

    interpolation_.setScale  ( interpolation_step_ );
    interpolation_.setOffset ( min_ );
  }

  void MixtureModel::updateMembers_()
  {
    InterpolationModel::updateMembers_();

    min_         = param_.getValue("bounding_box:min");
    max_        = param_.getValue("bounding_box:max");
    mix_prop_ =  param_.getValue("mix");

    statistics1_.setMean(param_.getValue("statistics:mean1"));
    statistics2_.setMean(param_.getValue("statistics:mean2"));
    statistics1_.setVariance(param_.getValue("statistics:variance1"));
    statistics2_.setVariance(param_.getValue("statistics:variance2"));

    setSamples();
  }

  void MixtureModel::setOffset(double offset)
  {
    double diff = offset - getInterpolation().getOffset();
    min_ += diff;
    max_ += diff;
    statistics1_.setMean(statistics1_.mean()+diff);
    statistics2_.setMean(statistics2_.mean()+diff);

    InterpolationModel::setOffset(offset);

    param_.setValue("bounding_box:min", min_);
    param_.setValue("bounding_box:max", max_);
    param_.setValue("statistics:mean1", statistics1_.mean());
    param_.setValue("statistics:mean2", statistics2_.mean());
  }

  MixtureModel::CoordinateType MixtureModel::getCenter() const
  {
    return ( (mix_prop_ * statistics1_.mean()) + ((1-mix_prop_) * statistics2_.mean()) / 2) ;
  }

}
