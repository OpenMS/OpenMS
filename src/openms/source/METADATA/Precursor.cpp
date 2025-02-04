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
      activation_methods_(std::move(rhs.activation_methods_)),
      activation_energy_(rhs.activation_energy_),
      window_low_(rhs.window_low_),
      window_up_(rhs.window_up_),
      drift_time_(rhs.drift_time_),
      drift_window_low_(rhs.drift_window_low_),
      drift_window_up_(rhs.drift_window_up_),
      drift_time_unit_(rhs.drift_time_unit_),
      charge_(rhs.charge_),
      possible_charge_states_(std::move(rhs.possible_charge_states_))
  {
  }

  bool Precursor::operator==(const Precursor& rhs) const
  {
    return activation_methods_ == rhs.activation_methods_ &&
           activation_energy_ == rhs.activation_energy_ &&
           window_low_ == rhs.window_low_ &&
           window_up_ == rhs.window_up_ &&
           drift_time_ == rhs.drift_time_ &&
           drift_window_up_ == rhs.drift_window_up_ &&
           drift_window_low_ == rhs.drift_window_low_ &&
           drift_time_unit_ == rhs.drift_time_unit_ &&
           charge_ == rhs.charge_ &&
           possible_charge_states_ == rhs.possible_charge_states_ &&
           Peak1D::operator==(rhs) &&
           CVTermList::operator==(rhs);
  }

  bool Precursor::operator!=(const Precursor& rhs) const
  {
    return !(operator==(rhs));
  }

  const set<Precursor::ActivationMethod>& Precursor::getActivationMethods() const
  {
    return activation_methods_;
  }

  set<Precursor::ActivationMethod>& Precursor::getActivationMethods()
  {
    return activation_methods_;
  }

  StringList Precursor::getActivationMethodsAsString() const
  {
    StringList am;
    am.reserve(activation_methods_.size());
    for (const auto& m : activation_methods_)
    {
      am.push_back(NamesOfActivationMethod[m]);
    }
    return am;
  }

  void Precursor::setActivationMethods(const set<Precursor::ActivationMethod> & activation_methods)
  {
    activation_methods_ = activation_methods;
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

  std::vector<Int> & Precursor::getPossibleChargeStates()
  {
    return possible_charge_states_;
  }

  const std::vector<Int> & Precursor::getPossibleChargeStates() const
  {
    return possible_charge_states_;
  }

  void Precursor::setPossibleChargeStates(const std::vector<Int> & possible_charge_states)
  {
    possible_charge_states_ = possible_charge_states;
  }

}

