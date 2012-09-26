// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#include <OpenMS/METADATA/IonSource.h>

using namespace std;

namespace OpenMS
{

  const std::string IonSource::NamesOfInletType[] = {"Unknown", "Direct", "Batch", "Chromatography", "Particle beam", "Membrane sparator", "Open split", "Jet separator", "Septum", "Reservoir", "Moving belt", "Moving wire", "Flow injection analysis", "Electro spray", "Thermo spray", "Infusion", "Continuous flow fast atom bombardment", "Inductively coupled plasma", "Membrane inlet", "Nanospray inlet"};

  const std::string IonSource::NamesOfIonizationMethod[] = {"Unknown", "Electrospray ionisation", "Electron ionization", "Chemical ionisation", "Fast atom bombardment", "Thermospray", "Laser desorption", "Field desorption", "Flame ionization", "Plasma desorption", "Secondary ion MS", "Thermal ionization", "Atmospheric pressure ionisation", "ISI", "Collsion induced decomposition", "Collsiona activated decomposition", "HN", "Atmospheric pressure chemical ionization", "Atmospheric pressure photo ionization", "Inductively coupled plasma", "Nano electrospray ionization", "Micro electrospray ionization", "Surface enhanced laser desorption ionization", "Surface enhanced neat desorption", "Fast ion bombardment", "Matrix-assisted laser desorption ionization", "Multiphoton ionization", "Desorption ionization", "Flowing afterglow", "Field ionization", "Glow discharge ionization", "Negative ion chemical ionization", "Neutralization reionization mass spectrometry", "Photoionization", "Pyrolysis mass spectrometry", "Resonance enhanced multiphoton ionization", "Adiabatic ionization", "Associative ionization", "Autodetachment", "Autoionization", "Charge exchange ionization", "Chemi-ionization", "Dissociative ionization", "Liquid secondary ionization", "Penning ionization", "Soft ionization", "Spark ionization", "Surface ionization", "Vertical ionization", "Atmospheric pressure matrix-assisted laser desorption ionization", "Desorption/ionization on silicon", "Surface-assisted laser desorption ionization"};

  const std::string IonSource::NamesOfPolarity[] = {"Unknown", "Positive", "Negative"};

  IonSource::IonSource() :
    MetaInfoInterface(),
    inlet_type_(INLETNULL),
    ionization_method_(IONMETHODNULL),
    polarity_(POLNULL),
    order_(0)
  {
  }

  IonSource::IonSource(const IonSource & source) :
    MetaInfoInterface(source),
    inlet_type_(source.inlet_type_),
    ionization_method_(source.ionization_method_),
    polarity_(source.polarity_),
    order_(source.order_)
  {
  }

  IonSource::~IonSource()
  {
  }

  IonSource & IonSource::operator=(const IonSource & source)
  {
    if (&source == this)
      return *this;

    order_ = source.order_;
    inlet_type_ = source.inlet_type_;
    ionization_method_ = source.ionization_method_;
    polarity_ = source.polarity_;
    MetaInfoInterface::operator=(source);

    return *this;
  }

  bool IonSource::operator==(const IonSource & rhs) const
  {
    return order_ == rhs.order_ &&
           inlet_type_ == rhs.inlet_type_ &&
           ionization_method_ == rhs.ionization_method_ &&
           polarity_ == rhs.polarity_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool IonSource::operator!=(const IonSource & rhs) const
  {
    return !(operator==(rhs));
  }

  IonSource::InletType IonSource::getInletType() const
  {
    return inlet_type_;
  }

  void IonSource::setInletType(IonSource::InletType inlet_type)
  {
    inlet_type_ = inlet_type;
  }

  IonSource::IonizationMethod IonSource::getIonizationMethod() const
  {
    return ionization_method_;
  }

  void IonSource::setIonizationMethod(IonSource::IonizationMethod ionization_type)
  {
    ionization_method_ = ionization_type;
  }

  IonSource::Polarity IonSource::getPolarity() const
  {
    return polarity_;
  }

  void IonSource::setPolarity(IonSource::Polarity polarity)
  {
    polarity_ = polarity;
  }

  Int IonSource::getOrder() const
  {
    return order_;
  }

  void IonSource::setOrder(Int order)
  {
    order_ = order;
  }

}
