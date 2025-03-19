// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/FEATUREFINDER/BaseModel_impl.h>

using namespace OpenMS;

// Default constructor
BaseModel::BaseModel() : DefaultParamHandler("BaseModel")
{
  defaults_.setValue("cutoff", 0.0, "Low intensity cutoff of the model. Peaks below this intensity are not considered part of the model.");
}

// Copy constructor
 BaseModel(const BaseModel& source) = default;

// Destructor
~BaseModel() override = default;

// Assignment operator
BaseModel& BaseModel::operator=(const BaseModel& source)
{
  if (&source == this)
    return *this;

  DefaultParamHandler::operator=(source);
  cut_off_ = source.cut_off_;

  return *this;
}

// Check if position @p pos is part of the model regarding the models cut-off.
bool BaseModel::isContained(const PositionType& pos) const
{
  return getIntensity(pos) >= cut_off_;
}

// Get cutoff value
BaseModel::IntensityType BaseModel::getCutOff() const
{
  return cut_off_;
}

// Set cutoff value
void BaseModel::setCutOff(IntensityType cut_off)
{
  cut_off_ = cut_off;
  param_.setValue("cutoff", cut_off_);
}

// Fill stream with reasonable set of samples from the model (i.e., for printing)
void BaseModel::getSamples(std::ostream& os)
{
  SamplesType samples;
  getSamples(samples);
  for (typename SamplesType::const_iterator it = samples.begin(); it != samples.end(); ++it)
  {
    os << *it << std::endl;
  }
}

// Update members
void BaseModel::updateMembers_()
{
  cut_off_ = (double)param_.getValue("cutoff");
}
