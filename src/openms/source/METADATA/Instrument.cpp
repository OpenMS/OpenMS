// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Instrument.h>

using namespace std;

namespace OpenMS
{

  const std::string Instrument::NamesOfIonOpticsType[] = {"Unknown", "magnetic deflection", "delayed extraction", "collision quadrupole", "selected ion flow tube", "time lag focusing", "reflectron", "einzel lens", "first stability region", "fringing field", "kinetic energy analyzer", "static field"};

  Instrument::Instrument() :
    MetaInfoInterface(),
    ion_optics_(UNKNOWN)
  {

  }

  Instrument::~Instrument() = default;

  bool Instrument::operator==(const Instrument & rhs) const
  {
    return software_ == rhs.software_ &&
           name_ == rhs.name_ &&
           vendor_ == rhs.vendor_ &&
           model_ == rhs.model_ &&
           customizations_ == rhs.customizations_ &&
           ion_sources_ == rhs.ion_sources_ &&
           mass_analyzers_ == rhs.mass_analyzers_ &&
           ion_detectors_ == rhs.ion_detectors_ &&
           ion_optics_ == rhs.ion_optics_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool Instrument::operator!=(const Instrument & rhs) const
  {
    return !(operator==(rhs));
  }

  const String & Instrument::getName() const
  {
    return name_;
  }

  void Instrument::setName(const String & name)
  {
    name_ = name;
  }

  const String & Instrument::getVendor() const
  {
    return vendor_;
  }

  void Instrument::setVendor(const String & vendor)
  {
    vendor_ = vendor;
  }

  const String & Instrument::getModel() const
  {
    return model_;
  }

  void Instrument::setModel(const String & model)
  {
    model_ = model;
  }

  const String & Instrument::getCustomizations() const
  {
    return customizations_;
  }

  void Instrument::setCustomizations(const String & customizations)
  {
    customizations_ = customizations;
  }

  const std::vector<IonSource> & Instrument::getIonSources() const
  {
    return ion_sources_;
  }

  std::vector<IonSource> & Instrument::getIonSources()
  {
    return ion_sources_;
  }

  void Instrument::setIonSources(const std::vector<IonSource> & ion_sources)
  {
    ion_sources_ = ion_sources;
  }

  const std::vector<MassAnalyzer> & Instrument::getMassAnalyzers() const
  {
    return mass_analyzers_;
  }

  std::vector<MassAnalyzer> & Instrument::getMassAnalyzers()
  {
    return mass_analyzers_;
  }

  void Instrument::setMassAnalyzers(const std::vector<MassAnalyzer> & mass_analyzers)
  {
    mass_analyzers_ = mass_analyzers;
  }

  const std::vector<IonDetector> & Instrument::getIonDetectors() const
  {
    return ion_detectors_;
  }

  std::vector<IonDetector> & Instrument::getIonDetectors()
  {
    return ion_detectors_;
  }

  void Instrument::setIonDetectors(const std::vector<IonDetector> & ion_detectors)
  {
    ion_detectors_ = ion_detectors;
  }

  const Software & Instrument::getSoftware() const
  {
    return software_;
  }

  Software & Instrument::getSoftware()
  {
    return software_;
  }

  void Instrument::setSoftware(const Software & software)
  {
    software_ = software;
  }

  Instrument::IonOpticsType Instrument::getIonOptics() const
  {
    return ion_optics_;
  }

  void Instrument::setIonOptics(Instrument::IonOpticsType ion_optics)
  {
    ion_optics_ = ion_optics;
  }

}

