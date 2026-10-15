// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/HPLC.h>

#include <utility>

using namespace std;

namespace OpenMS
{

  HPLC::HPLC() :
    instrument_(),
    column_(),
    temperature_(21),
    pressure_(0),
    flux_(0),
    comment_(),
    gradient_()
  {
  }

  HPLC::~HPLC() = default;

  bool HPLC::operator==(const HPLC & rhs) const
  {
    return (instrument_ == rhs.instrument_) &&
           (column_ == rhs.column_) &&
           (temperature_ == rhs.temperature_) &&
           (pressure_ == rhs.pressure_) &&
           (flux_ == rhs.flux_) &&
           (comment_ == rhs.comment_) &&
           (gradient_ == rhs.gradient_);
  }

  bool HPLC::operator!=(const HPLC & rhs) const
  {
    return !(operator==(rhs));
  }

  const std::string & HPLC::getInstrument() const
  {
    return instrument_;
  }

  void HPLC::setInstrument(const std::string & instrument)
  {
    instrument_ = instrument;
  }

  const std::string & HPLC::getColumn() const
  {
    return column_;
  }

  void HPLC::setColumn(const std::string & column)
  {
    column_ = column;
  }

  Int HPLC::getTemperature() const
  {
    return temperature_;
  }

  void HPLC::setTemperature(Int temperature)
  {
    temperature_ = temperature;
  }

  UInt HPLC::getPressure() const
  {
    return pressure_;
  }

  void HPLC::setPressure(UInt pressure)
  {
    pressure_ = pressure;
  }

  UInt HPLC::getFlux() const
  {
    return flux_;
  }

  void HPLC::setFlux(UInt flux)
  {
    flux_ = flux;
  }

  std::string HPLC::getComment() const
  {
    return comment_;
  }

  void HPLC::setComment(std::string comment)
  {
    comment_ = std::move(comment);
  }

  Gradient & HPLC::getGradient()
  {
    return gradient_;
  }

  const Gradient & HPLC::getGradient() const
  {
    return gradient_;
  }

  void HPLC::setGradient(const Gradient & gradient)
  {
    gradient_ = gradient;
  }

} // namespace OpenMS

