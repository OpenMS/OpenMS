// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IonSource.h>

using namespace std;

namespace OpenMS
{

  const std::string IonSource::NamesOfInletType[] = {"Unknown", "Direct", "Batch", "Chromatography", "Particle beam", "Membrane sparator", "Open split", "Jet separator", "Septum", "Reservoir", "Moving belt", "Moving wire", "Flow injection analysis", "Electro spray", "Thermo spray", "Infusion", "Continuous flow fast atom bombardment", "Inductively coupled plasma", "Membrane inlet", "Nanospray inlet"};

  const std::string IonSource::NamesOfIonizationMethod[] = {"Unknown", "Electrospray ionisation", "Electron ionization", "Chemical ionisation", "Fast atom bombardment", "Thermospray", "Laser desorption", "Field desorption", "Flame ionization", "Plasma desorption", "Secondary ion MS", "Thermal ionization", "Atmospheric pressure ionisation", "ISI", "Collsion induced decomposition", "Collsiona activated decomposition", "HN", "Atmospheric pressure chemical ionization", "Atmospheric pressure photo ionization", "Inductively coupled plasma", "Nano electrospray ionization", "Micro electrospray ionization", "Surface enhanced laser desorption ionization", "Surface enhanced neat desorption", "Fast ion bombardment", "Matrix-assisted laser desorption ionization", "Multiphoton ionization", "Desorption ionization", "Flowing afterglow", "Field ionization", "Glow discharge ionization", "Negative ion chemical ionization", "Neutralization reionization mass spectrometry", "Photoionization", "Pyrolysis mass spectrometry", "Resonance enhanced multiphoton ionization", "Adiabatic ionization", "Associative ionization", "Autodetachment", "Autoionization", "Charge exchange ionization", "Chemi-ionization", "Dissociative ionization", "Liquid secondary ionization", "Penning ionization", "Soft ionization", "Spark ionization", "Surface ionization", "Vertical ionization", "Atmospheric pressure matrix-assisted laser desorption ionization", "Desorption/ionization on silicon", "Surface-assisted laser desorption ionization"};

  const std::string IonSource::NamesOfPolarity[] = {"unknown", "positive", "negative"};

  IonSource::IonSource() :
    MetaInfoInterface(),
    inlet_type_(INLETNULL),
    ionization_method_(IONMETHODNULL),
    polarity_(POLNULL),
    order_(0)
  {
  }

  IonSource::~IonSource() = default;

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

