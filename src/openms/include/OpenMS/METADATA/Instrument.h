// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/Software.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief Description of a MS instrument

      It contains information like vendor, model, ion source(s), mass analyzer(s), and ion detector(s).

      The parts (IonSource, MassAnalyzer, IonDetector) all have a @em order member.
      The order can be ignored, as long the instrument has this default setup:
      - one ion source
      - one or many mass analyzers
      - one ion detector

      For more complex instruments, the order should be defined.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Instrument :
    public MetaInfoInterface
  {

public:

    /// ion optics type
    enum IonOpticsType
    {
      UNKNOWN,                  ///< unknown
      MAGNETIC_DEFLECTION,      ///< magnetic deflection
      DELAYED_EXTRACTION,       ///< delayed extraction
      COLLISION_QUADRUPOLE,     ///< collision quadrupole
      SELECTED_ION_FLOW_TUBE,   ///< selected ion flow tube
      TIME_LAG_FOCUSING,        ///< time lag focusing
      REFLECTRON,               ///< reflectron
      EINZEL_LENS,              ///< einzel lens
      FIRST_STABILITY_REGION,   ///< first stability region
      FRINGING_FIELD,           ///< fringing field
      KINETIC_ENERGY_ANALYZER,  ///< kinetic energy analyzer
      STATIC_FIELD,             ///< static field
      SIZE_OF_IONOPTICSTYPE
    };

    /// Names of inlet types
    static const std::string NamesOfIonOpticsType[SIZE_OF_IONOPTICSTYPE];

    /// Constructor
    Instrument();
    /// Copy constructor
    Instrument(const Instrument &) = default;
    /// Move constructor
    Instrument(Instrument&&) = default;
    /// Destructor
    ~Instrument();

    /// Assignment operator
    Instrument & operator=(const Instrument &) = default;
    /// Move assignment operator
    Instrument& operator=(Instrument&&) & = default;

    /// Equality operator
    bool operator==(const Instrument & rhs) const;
    /// Equality operator
    bool operator!=(const Instrument & rhs) const;

    /// returns the name of the instrument
    const String & getName() const;
    /// sets the name of the instrument
    void setName(const String & name);

    /// returns the instrument vendor
    const String & getVendor() const;
    /// sets the instrument vendor
    void setVendor(const String & vendor);

    /// returns the instrument model
    const String & getModel() const;
    /// sets the instrument model
    void setModel(const String & model);

    /// returns a description of customizations
    const String & getCustomizations() const;
    /// sets a description of customizations
    void setCustomizations(const String & customizations);

    /// returns a const reference to the ion source list
    const std::vector<IonSource> & getIonSources() const;
    /// returns a mutable reference to the ion source list
    std::vector<IonSource> & getIonSources();
    /// sets the ion source list
    void setIonSources(const std::vector<IonSource> & ion_sources);

    /// returns a const reference to the mass analyzer list
    const std::vector<MassAnalyzer> & getMassAnalyzers() const;
    /// returns a mutable reference to the mass analyzer list
    std::vector<MassAnalyzer> & getMassAnalyzers();
    /// sets the mass analyzer list
    void setMassAnalyzers(const std::vector<MassAnalyzer> & mass_analyzers);

    /// returns a const reference to the ion detector list
    const std::vector<IonDetector> & getIonDetectors() const;
    /// returns a mutable reference to the ion detector list
    std::vector<IonDetector> & getIonDetectors();
    /// sets the ion detector list
    void setIonDetectors(const std::vector<IonDetector> & ion_detectors);

    /// returns a const reference to the instrument software
    const Software & getSoftware() const;
    /// returns a mutable reference to the instrument software
    Software & getSoftware();
    /// sets the instrument software
    void setSoftware(const Software & software);

    /// returns the ion optics type
    IonOpticsType getIonOptics() const;
    /// sets the ion optics type
    void setIonOptics(IonOpticsType ion_optics);

protected:

    String name_;
    String vendor_;
    String model_;
    String customizations_;
    std::vector<IonSource> ion_sources_;
    std::vector<MassAnalyzer> mass_analyzers_;
    std::vector<IonDetector> ion_detectors_;
    Software software_;
    IonOpticsType ion_optics_;
  };
} // namespace OpenMS

