// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Mathias Walzer $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Precursor.h>

using namespace std;

namespace OpenMS
{

  const std::string Precursor::NamesOfActivationMethod[] = {"Collision-induced dissociation", "Post-source decay", "Plasma desorption", "Surface-induced dissociation", "Blackbody infrared radiative dissociation", "Electron capture dissociation", "Infrared multiphoton dissociation", "Sustained off-resonance irradiation", "High-energy collision-induced dissociation", "Low-energy collision-induced dissociation", "Photodissociation", "Electron transfer dissociation", "Pulsed q dissociation"};
  const std::string Precursor::NamesOfActivationMethodShort[] = { "CID", "PSD", "PD", "SID", "BIRD", "ECD", "IMD", "SORI", "HCID", "LCID", "PHD", "ETD", "PQD" };

  Precursor::Precursor() :
    CVTermList(),
    Peak1D(),
    activation_methods_(),
    activation_energy_(0.0),
    window_low_(0.0),
    window_up_(0.0),
    drift_time_(-1),
    drift_window_low_(0.0),
    drift_window_up_(0.0),
    drift_time_unit_(Precursor::DriftTimeUnit::NONE),
    charge_(0),
    possible_charge_states_()
  {
  }

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

  Precursor::~Precursor()
  {
  }

  bool Precursor::operator==(const Precursor & rhs) const
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

  bool Precursor::operator!=(const Precursor & rhs) const
  {
    return !(operator==(rhs));
  }

  const set<Precursor::ActivationMethod> & Precursor::getActivationMethods() const
  {
    return activation_methods_;
  }

  set<Precursor::ActivationMethod> & Precursor::getActivationMethods()
  {
    return activation_methods_;
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
    if (bound < 0) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor::setIsolationWindowLowerOffset() received a negative lower offset", String(bound));
    window_low_ = bound;
  }

  double Precursor::getIsolationWindowUpperOffset() const
  {
    return window_up_;
  }

  void Precursor::setIsolationWindowUpperOffset(double bound)
  {
    if (bound < 0) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor::setIsolationWindowUpperOffset() received a negative lower offset", String(bound));
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

  Precursor::DriftTimeUnit Precursor::getDriftTimeUnit() const
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

