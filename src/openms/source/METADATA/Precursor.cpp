// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Precursor.h>

using namespace std;

namespace OpenMS
{

  const std::string Precursor::NamesOfActivationMethod[] = {
    "Unknown",
    "Collision-induced dissociation", 
    "Post-source decay", 
    "Plasma desorption", 
    "Surface-induced dissociation", 
    "Blackbody infrared radiative dissociation", 
    "Electron capture dissociation", 
    "Infrared multiphoton dissociation", 
    "Sustained off-resonance irradiation", 
    "High-energy collision-induced dissociation", 
    "Low-energy collision-induced dissociation", 
    "Photodissociation", 
    "Electron transfer dissociation", 
    "Electron transfer and collision-induced dissociation",
    "Electron transfer and higher-energy collision dissociation",
    "Pulsed q dissociation",
    "trap-type collision-induced dissociation",
    "beam-type collision-induced dissociation", // == HCD
    "in-source collision-induced dissociation",
    "Bruker proprietary method"
    };
  
  const std::string Precursor::NamesOfActivationMethodShort[] = {
    "UNKNOWN",
    "CID", 
    "PSD", 
    "PD", 
    "SID", 
    "BIRD", 
    "ECD", 
    "IMD", 
    "SORI", 
    "HCID", 
    "LCID", 
    "PHD", 
    "ETD", 
    "ETciD",
    "EThcD",
    "PQD",
    "TRAP",
    "HCD",
    "INSOURCE",
    "LIFT"
  };

  Precursor::Precursor(Precursor&& rhs) noexcept :
      CVTermList(std::move(rhs)),
      Peak1D(std::move(rhs)),
      activation_method_(rhs.activation_method_),
      activation_energy_(rhs.activation_energy_),
      window_low_(rhs.window_low_),
      window_up_(rhs.window_up_),
      drift_time_(rhs.drift_time_),
      drift_window_low_(rhs.drift_window_low_),
      drift_window_up_(rhs.drift_window_up_),
      drift_time_unit_(rhs.drift_time_unit_),
      charge_(rhs.charge_)
  {
  }

  bool Precursor::operator==(const Precursor& rhs) const
  {
    return activation_method_ == rhs.activation_method_ &&
           activation_energy_ == rhs.activation_energy_ &&
           window_low_ == rhs.window_low_ &&
           window_up_ == rhs.window_up_ &&
           drift_time_ == rhs.drift_time_ &&
           drift_window_up_ == rhs.drift_window_up_ &&
           drift_window_low_ == rhs.drift_window_low_ &&
           drift_time_unit_ == rhs.drift_time_unit_ &&
           charge_ == rhs.charge_ &&
           Peak1D::operator==(rhs) &&
           CVTermList::operator==(rhs);
  }

  bool Precursor::operator!=(const Precursor& rhs) const
  {
    return !(operator==(rhs));
  }

  Precursor::ActivationMethod Precursor::getActivationMethod() const
  {
    return activation_method_;
  }

  void Precursor::setActivationMethod(Precursor::ActivationMethod activation_method)
  {
    activation_method_ = activation_method;
  }

  StringList Precursor::getActivationMethodsAsString() const
  {
    StringList am;
    // Only add the activation method if it's valid (not SIZE_OF_ACTIVATIONMETHOD)
    if (activation_method_ != UNKNOWN)
    {
      am.push_back(NamesOfActivationMethod[activation_method_]);
    }
    return am;
  }


  double Precursor::getActivationEnergy() const
  {
    return activation_energy_;
  }

  void Precursor::setActivationEnergy(double activation_energy)
  {
    activation_energy_ = activation_energy;
  }

  double Precursor::getIsolationWindowLowerOffset() const
  {
    return window_low_;
  }

  void Precursor::setIsolationWindowLowerOffset(double bound)
  {
    if (bound < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor::setIsolationWindowLowerOffset() received a negative lower offset", String(bound));
    }
    window_low_ = bound;
  }

  double Precursor::getIsolationWindowUpperOffset() const
  {
    return window_up_;
  }

  void Precursor::setIsolationWindowUpperOffset(double bound)
  {
    if (bound < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor::setIsolationWindowUpperOffset() received a negative lower offset", String(bound));
    }
    window_up_ = bound;
  }

  double Precursor::getDriftTime() const
  {
    return drift_time_;
  }

  void Precursor::setDriftTime(double drift_time)
  {
    drift_time_ = drift_time;
  }

  DriftTimeUnit Precursor::getDriftTimeUnit() const
  {
    return drift_time_unit_;
  }

  void Precursor::setDriftTimeUnit(DriftTimeUnit dt)
  {
    drift_time_unit_ = dt;
  }

  double Precursor::getDriftTimeWindowLowerOffset() const
  {
    return drift_window_low_;
  }

  void Precursor::setDriftTimeWindowLowerOffset(double bound)
  {
    OPENMS_PRECONDITION(bound >= 0, "Relative drift time offset needs to be positive.")
    drift_window_low_ = bound;
  }

  double Precursor::getDriftTimeWindowUpperOffset() const
  {
    return drift_window_up_;
  }

  void Precursor::setDriftTimeWindowUpperOffset(double bound)
  {
    OPENMS_PRECONDITION(bound >= 0, "Relative drift time offset needs to be positive.")
    drift_window_up_ = bound;
  }

  Int Precursor::getCharge() const
  {
    return charge_;
  }

  void Precursor::setCharge(Int charge)
  {
    charge_ = charge;
  }


}

