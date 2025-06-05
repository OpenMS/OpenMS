// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Sample.h>

using namespace std;

namespace OpenMS
{

  const std::string Sample::NamesOfSampleState[] = {"Unknown", "solid", "liquid", "gas", "solution", "emulsion", "suspension"};

  Sample::Sample() :
    MetaInfoInterface(),
    state_(SAMPLENULL),
    mass_(0.0),
    volume_(0.0),
    concentration_(0.0)
  {
  }

  Sample::Sample(const Sample & source) :
    MetaInfoInterface(source),
    name_(source.name_),
    number_(source.number_),
    comment_(source.comment_),
    organism_(source.organism_),
    state_(source.state_),
    mass_(source.mass_),
    volume_(source.volume_),
    concentration_(source.concentration_),
    subsamples_(source.subsamples_)
  {
  }

  Sample::~Sample() = default;

  Sample & Sample::operator=(const Sample & source)
  {
    if (&source == this)
    {
      return *this;
    }

    name_ = source.name_;
    number_ = source.number_;
    comment_ = source.comment_;
    organism_ = source.organism_;
    state_ = source.state_;
    mass_ = source.mass_;
    volume_ = source.volume_;
    concentration_ = source.concentration_;
    subsamples_ = source.subsamples_;
    MetaInfoInterface::operator=(source);
    return *this;
  }

  bool Sample::operator==(const Sample & rhs) const
  {
    if (
      name_ != rhs.name_ ||
      number_ != rhs.number_ ||
      comment_ != rhs.comment_ ||
      organism_ != rhs.organism_ ||
      state_ != rhs.state_ ||
      mass_ != rhs.mass_ ||
      volume_ != rhs.volume_ ||
      concentration_ != rhs.concentration_ ||
      subsamples_ != rhs.subsamples_ ||
      MetaInfoInterface::operator!=(rhs))
    {
      return false;
    }
    return true;
  }

  const String & Sample::getName() const
  {
    return name_;
  }

  void Sample::setName(const String & name)
  {
    name_ = name;
  }

  const String & Sample::getOrganism() const
  {
    return organism_;
  }

  void Sample::setOrganism(const String & organism)
  {
    organism_ = organism;
  }

  const String & Sample::getNumber() const
  {
    return number_;
  }

  void Sample::setNumber(const String & number)
  {
    number_ = number;
  }

  const String & Sample::getComment() const
  {
    return comment_;
  }

  void Sample::setComment(const String & comment)
  {
    comment_ = comment;
  }

  Sample::SampleState Sample::getState() const
  {
    return state_;
  }

  void Sample::setState(Sample::SampleState state)
  {
    state_ = state;
  }

  double Sample::getMass() const
  {
    return mass_;
  }

  void Sample::setMass(double mass)
  {
    mass_ = mass;
  }

  double Sample::getVolume() const
  {
    return volume_;
  }

  void Sample::setVolume(double volume)
  {
    volume_ = volume;
  }

  double Sample::getConcentration() const
  {
    return concentration_;
  }

  void Sample::setConcentration(double concentration)
  {
    concentration_ = concentration;
  }

  const vector<Sample> & Sample::getSubsamples() const
  {
    return subsamples_;
  }

  vector<Sample> & Sample::getSubsamples()
  {
    return subsamples_;
  }

  void Sample::setSubsamples(const vector<Sample> & subsamples)
  {
    subsamples_ = subsamples;
  }
}

