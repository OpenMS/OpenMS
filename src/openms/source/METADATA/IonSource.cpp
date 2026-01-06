// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IonSource.h>
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

  StringList IonSource::getAllNamesOfInletType()
  {
    StringList names;
    names.reserve(SIZE_OF_INLETTYPE);
    for (size_t i = 0; i < SIZE_OF_INLETTYPE; ++i)
    {
      names.push_back(NamesOfInletType[i]);
    }
    return names;
  }

  StringList IonSource::getAllNamesOfIonizationMethod()
  {
    StringList names;
    names.reserve(SIZE_OF_IONIZATIONMETHOD);
    for (size_t i = 0; i < SIZE_OF_IONIZATIONMETHOD; ++i)
    {
      names.push_back(NamesOfIonizationMethod[i]);
    }
    return names;
  }

  StringList IonSource::getAllNamesOfPolarity()
  {
    StringList names;
    names.reserve(SIZE_OF_POLARITY);
    for (size_t i = 0; i < SIZE_OF_POLARITY; ++i)
    {
      names.push_back(NamesOfPolarity[i]);
    }
    return names;
  }

  const std::string& IonSource::inletTypeToString(InletType type)
  {
    if (type == SIZE_OF_INLETTYPE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_INLETTYPE");
    }
    return NamesOfInletType[static_cast<size_t>(type)];
  }

  IonSource::InletType IonSource::toInletType(const std::string& name)
  {
    auto first = &NamesOfInletType[0];
    auto last = &NamesOfInletType[SIZE_OF_INLETTYPE];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<InletType>(it - first);
  }

  const std::string& IonSource::ionizationMethodToString(IonizationMethod method)
  {
    if (method == SIZE_OF_IONIZATIONMETHOD)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_IONIZATIONMETHOD");
    }
    return NamesOfIonizationMethod[static_cast<size_t>(method)];
  }

  IonSource::IonizationMethod IonSource::toIonizationMethod(const std::string& name)
  {
    auto first = &NamesOfIonizationMethod[0];
    auto last = &NamesOfIonizationMethod[SIZE_OF_IONIZATIONMETHOD];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<IonizationMethod>(it - first);
  }

  const std::string& IonSource::polarityToString(Polarity polarity)
  {
    if (polarity == SIZE_OF_POLARITY)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_POLARITY");
    }
    return NamesOfPolarity[static_cast<size_t>(polarity)];
  }

  IonSource::Polarity IonSource::toPolarity(const std::string& name)
  {
    auto first = &NamesOfPolarity[0];
    auto last = &NamesOfPolarity[SIZE_OF_POLARITY];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<Polarity>(it - first);
  }

}

