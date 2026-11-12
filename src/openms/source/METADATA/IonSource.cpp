// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{

  const std::string IonSource::NamesOfInletType[] = {"Unknown", "Direct", "Batch", "Chromatography", "Particle beam", "Membrane sparator", "Open split", "Jet separator", "Septum", "Reservoir", "Moving belt", "Moving wire", "Flow injection analysis", "Electro spray", "Thermo spray", "Infusion", "Continuous flow fast atom bombardment", "Inductively coupled plasma", "Membrane inlet", "Nanospray inlet"};

  const std::string IonSource::NamesOfIonizationMethod[] = {"Unknown", "Electrospray ionisation", "Electron ionization", "Chemical ionisation", "Fast atom bombardment", "Thermospray", "Laser desorption", "Field desorption", "Flame ionization", "Plasma desorption", "Secondary ion MS", "Thermal ionization", "Atmospheric pressure ionisation", "ISI", "Collsion induced decomposition", "Collsiona activated decomposition", "HN", "Atmospheric pressure chemical ionization", "Atmospheric pressure photo ionization", "Inductively coupled plasma", "Nano electrospray ionization", "Micro electrospray ionization", "Surface enhanced laser desorption ionization", "Surface enhanced neat desorption", "Fast ion bombardment", "Matrix-assisted laser desorption ionization", "Multiphoton ionization", "Desorption ionization", "Flowing afterglow", "Field ionization", "Glow discharge ionization", "Negative ion chemical ionization", "Neutralization reionization mass spectrometry", "Photoionization", "Pyrolysis mass spectrometry", "Resonance enhanced multiphoton ionization", "Adiabatic ionization", "Associative ionization", "Autodetachment", "Autoionization", "Charge exchange ionization", "Chemi-ionization", "Dissociative ionization", "Liquid secondary ionization", "Penning ionization", "Soft ionization", "Spark ionization", "Surface ionization", "Vertical ionization", "Atmospheric pressure matrix-assisted laser desorption ionization", "Desorption/ionization on silicon", "Surface-assisted laser desorption ionization"};

  const std::string IonSource::NamesOfPolarity[] = {"unknown", "positive", "negative"};

  IonSource::IonSource() :
    MetaInfoInterface(),
    inlet_type_(InletType::INLETNULL),
    ionization_method_(IonizationMethod::IONMETHODNULL),
    polarity_(Polarity::POLNULL),
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

  StringList IonSource::getAllNamesOfInletType()
  {
    return Helpers::enumNames(NamesOfInletType, InletType::SIZE_OF_INLETTYPE);
  }

  StringList IonSource::getAllNamesOfIonizationMethod()
  {
    return Helpers::enumNames(NamesOfIonizationMethod, IonizationMethod::SIZE_OF_IONIZATIONMETHOD);
  }

  StringList IonSource::getAllNamesOfPolarity()
  {
    return Helpers::enumNames(NamesOfPolarity, Polarity::SIZE_OF_POLARITY);
  }

  const std::string& IonSource::inletTypeToString(InletType type)
  {
    return Helpers::enumToName(NamesOfInletType, type, InletType::SIZE_OF_INLETTYPE, "SIZE_OF_INLETTYPE", OPENMS_PRETTY_FUNCTION);
  }

  IonSource::InletType IonSource::toInletType(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfInletType, InletType::SIZE_OF_INLETTYPE, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

  const std::string& IonSource::ionizationMethodToString(IonizationMethod method)
  {
    return Helpers::enumToName(NamesOfIonizationMethod, method, IonizationMethod::SIZE_OF_IONIZATIONMETHOD, "SIZE_OF_IONIZATIONMETHOD", OPENMS_PRETTY_FUNCTION);
  }

  IonSource::IonizationMethod IonSource::toIonizationMethod(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfIonizationMethod, IonizationMethod::SIZE_OF_IONIZATIONMETHOD, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

  const std::string& IonSource::polarityToString(Polarity polarity)
  {
    return Helpers::enumToName(NamesOfPolarity, polarity, Polarity::SIZE_OF_POLARITY, "SIZE_OF_POLARITY", OPENMS_PRETTY_FUNCTION);
  }

  IonSource::Polarity IonSource::toPolarity(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfPolarity, Polarity::SIZE_OF_POLARITY, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

}

