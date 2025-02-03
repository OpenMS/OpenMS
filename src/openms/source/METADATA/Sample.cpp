// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Sample.h>

#include <OpenMS/METADATA/SampleTreatment.h>

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
    // delete old treatments
    for (auto& it : treatments_)
    {
      delete it;
    }
    treatments_.clear();
    // copy treatments
    for (const auto & it : source.treatments_)
    {
      treatments_.push_back(it->clone());
    }
  }

  Sample::~Sample()
  {
    for (auto& it : treatments_)
    {
      delete it;
    }
  }

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

    // delete old treatments
    for (auto& it : treatments_)
    {
      delete it;
    }
    treatments_.clear();
    // copy treatments
    for (const auto& it : source.treatments_)
    {
      treatments_.push_back(it->clone());
    }
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
      MetaInfoInterface::operator!=(rhs) ||
      treatments_.size() != rhs.treatments_.size()
      )
    {
      return false;
    }

    // treatments
    auto it2 = rhs.treatments_.begin();
    for (auto it = treatments_.begin(); it != treatments_.end(); ++it, ++it2)
    {
      if (*it != *it2)
      {
        return false;
      }
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

  void Sample::addTreatment(const SampleTreatment & treatment, Int before_position)
  {
    if (before_position > Int(treatments_.size()))
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, before_position, treatments_.size());
    }
    list<SampleTreatment *>::iterator it;
    if (before_position >= 0)
    {
      it = treatments_.begin();
      for (Int i = 0; i < before_position; ++i)
      {
        ++it;
      }
    }
    else
    {
      it = treatments_.end();
    }
    SampleTreatment * tmp = treatment.clone();
    treatments_.insert(it, tmp);
  }

  const SampleTreatment & Sample::getTreatment(UInt position) const
  {
    if (position >= treatments_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, position, treatments_.size());
    }
    list<SampleTreatment *>::const_iterator it = treatments_.begin();
    for (Size i = 0; i < position; ++i)
    {
      ++it;
    }
    return **it;
  }

  SampleTreatment & Sample::getTreatment(UInt position)
  {
    if (position >= treatments_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, position, treatments_.size());
    }
    list<SampleTreatment *>::const_iterator it = treatments_.begin();
    for (Size i = 0; i < position; ++i)
    {
      ++it;
    }
    return **it;
  }

  void Sample::removeTreatment(UInt position)
  {
    if (position >= treatments_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, position, treatments_.size());
    }
    list<SampleTreatment *>::iterator it = treatments_.begin();
    for (Size i = 0; i < position; ++i)
    {
      ++it;
    }
    delete(*it);
    treatments_.erase(it);
  }

  Int Sample::countTreatments() const
  {
    return (Int)treatments_.size();
  }

}

