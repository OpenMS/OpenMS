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
// $Maintainer: Mathias Walzer $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Precursor.h>

using namespace std;

namespace OpenMS
{

  const std::string Precursor::NamesOfActivationMethod[] = {"Collision-induced dissociation", "Post-source decay", "Plasma desorption", "Surface-induced dissociation", "Blackbody infrared radiative dissociation", "Electron capture dissociation", "Infrared multiphoton dissociation", "Sustained off-resonance irradiation", "High-energy collision-induced dissociation", "Low-energy collision-induced dissociation", "Photodissociation", "Electron transfer dissociation", "Pulsed q dissociation"};

  Precursor::Precursor() :
    CVTermList(),
    Peak1D(),
    activation_methods_(),
    activation_energy_(0.0),
    window_low_(0.0),
    window_up_(0.0),
    charge_(0),
    possible_charge_states_()
  {
  }

  Precursor::Precursor(const Precursor & source) :
    CVTermList(source),
    Peak1D(source),
    activation_methods_(source.activation_methods_),
    activation_energy_(source.activation_energy_),
    window_low_(source.window_low_),
    window_up_(source.window_up_),
    charge_(source.charge_),
    possible_charge_states_(source.possible_charge_states_)
  {
  }

  Precursor::~Precursor()
  {
  }

  Precursor & Precursor::operator=(const Precursor & source)
  {
    if (&source == this)
      return *this;

    CVTermList::operator=(source);
    Peak1D::operator=(source);
    activation_methods_ = source.activation_methods_;
    activation_energy_ = source.activation_energy_;
    window_low_ = source.window_low_;
    window_up_ = source.window_up_;
    charge_ = source.charge_;
    possible_charge_states_ = source.possible_charge_states_;

    return *this;
  }

  bool Precursor::operator==(const Precursor & rhs) const
  {
    return activation_methods_ == rhs.activation_methods_ &&
           activation_energy_ == rhs.activation_energy_ &&
           window_low_ == rhs.window_low_ &&
           window_up_ == rhs.window_up_ &&
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

  DoubleReal Precursor::getActivationEnergy() const
  {
    return activation_energy_;
  }

  void Precursor::setActivationEnergy(DoubleReal activation_energy)
  {
    activation_energy_ = activation_energy;
  }

  DoubleReal Precursor::getIsolationWindowLowerOffset() const
  {
    return window_low_;
  }

  void Precursor::setIsolationWindowLowerOffset(DoubleReal bound)
  {
    window_low_ = bound;
  }

  DoubleReal Precursor::getIsolationWindowUpperOffset() const
  {
    return window_up_;
  }

  void Precursor::setIsolationWindowUpperOffset(DoubleReal bound)
  {
    window_up_ = bound;
  }

  Int Precursor::getCharge() const
  {
    return charge_;
  }

  void Precursor::setCharge(Int charge)
  {
    charge_ = charge;
    return;
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
