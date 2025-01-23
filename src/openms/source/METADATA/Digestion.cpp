// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Digestion.h>

using namespace std;

namespace OpenMS
{


  Digestion::Digestion() :
    SampleTreatment("Digestion"),
    enzyme_(""),
    digestion_time_(0.0),
    temperature_(0.0),
    ph_(0.0)
  {

  }

  Digestion::~Digestion() = default;

  SampleTreatment * Digestion::clone() const
  {
    SampleTreatment * tmp = new Digestion(*this);
    return tmp;
  }

  bool Digestion::operator==(const SampleTreatment & rhs) const
  {
    if (type_ != rhs.getType())
    {
      return false;
    }
    const Digestion * tmp = dynamic_cast<const Digestion *>(&rhs);
    return SampleTreatment::operator==(* tmp) &&
           enzyme_ == tmp->enzyme_ &&
           digestion_time_ == tmp->digestion_time_ &&
           temperature_ == tmp->temperature_ &&
           ph_ == tmp->ph_;
  }

  const String & Digestion::getEnzyme() const
  {
    return enzyme_;
  }

  void Digestion::setEnzyme(const String & enzyme)
  {
    enzyme_ = enzyme;
  }

  double Digestion::getDigestionTime() const
  {
    return digestion_time_;
  }

  void Digestion::setDigestionTime(double digestion_time)
  {
    digestion_time_ = digestion_time;
  }

  double Digestion::getTemperature() const
  {
    return temperature_;
  }

  void Digestion::setTemperature(double temperature)
  {
    temperature_ = temperature;
  }

  double Digestion::getPh() const
  {
    return ph_;
  }

  void Digestion::setPh(double ph)
  {
    ph_ = ph;
  }

}

