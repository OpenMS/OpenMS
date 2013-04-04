// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/SampleTreatment.h>

#include <iostream>

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
    //delete old treatments
    for (list<SampleTreatment *>::iterator it = treatments_.begin(); it != treatments_.end(); ++it)
    {
      delete(*it);
    }
    treatments_.clear();
    //copy treatments
    for (list<SampleTreatment *>::const_iterator it = source.treatments_.begin(); it != source.treatments_.end(); ++it)
    {
      treatments_.push_back((*it)->clone());
    }
  }

  Sample::~Sample()
  {
    for (list<SampleTreatment *>::iterator it = treatments_.begin(); it != treatments_.end(); ++it)
    {
      delete(*it);
    }
  }

  Sample & Sample::operator=(const Sample & source)
  {
    if (&source == this)
      return *this;

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
    //delete old treatments
    for (list<SampleTreatment *>::iterator it = treatments_.begin(); it != treatments_.end(); ++it)
    {
      delete(*it);
    }
    treatments_.clear();
    //copy treatments
    for (list<SampleTreatment *>::const_iterator it = source.treatments_.begin(); it != source.treatments_.end(); ++it)
    {
      treatments_.push_back((*it)->clone());
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

    //treatments
    list<SampleTreatment *>::const_iterator it2 = rhs.treatments_.begin();
    for (list<SampleTreatment *>::const_iterator it = treatments_.begin(); it != treatments_.end(); ++it, ++it2)
    {
      if (*it != *it2)
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

  DoubleReal Sample::getMass() const
  {
    return mass_;
  }

  void Sample::setMass(DoubleReal mass)
  {
    mass_ = mass;
  }

  DoubleReal Sample::getVolume() const
  {
    return volume_;
  }

  void Sample::setVolume(DoubleReal volume)
  {
    volume_ = volume;
  }

  DoubleReal Sample::getConcentration() const
  {
    return concentration_;
  }

  void Sample::setConcentration(DoubleReal concentration)
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, before_position, treatments_.size());
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, position, treatments_.size());
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, position, treatments_.size());
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, position, treatments_.size());
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
