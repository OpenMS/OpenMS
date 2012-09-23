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

#ifndef OPENMS_METADATA_INSTRUMENT_H
#define OPENMS_METADATA_INSTRUMENT_H

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

      For more complex instuments, the order should be defined.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Instrument :
    public MetaInfoInterface
  {

public:

    /// ion optics type
    enum IonOpticsType
    {
      UNKNOWN,                                          ///< unknown
      MAGNETIC_DEFLECTION,                      ///< magnetic deflection
      DELAYED_EXTRACTION,                       ///< delayed extraction
      COLLISION_QUADRUPOLE,                 ///< collision quadrupole
      SELECTED_ION_FLOW_TUBE,               ///< selected ion flow tube
      TIME_LAG_FOCUSING,                        ///< time lag focusing
      REFLECTRON,                                       ///< reflectron
      EINZEL_LENS,                                      ///< einzel lens
      FIRST_STABILITY_REGION,               ///< first stability region
      FRINGING_FIELD,                               ///< fringing field
      KINETIC_ENERGY_ANALYZER,              ///< kinetic energy analyzer
      STATIC_FIELD,                                     ///< static field
      SIZE_OF_IONOPTICSTYPE
    };
    /// Names of inlet types
    static const std::string NamesOfIonOpticsType[SIZE_OF_IONOPTICSTYPE];

    /// Constructor
    Instrument();
    /// Copy constructor
    Instrument(const Instrument & source);
    /// Destructor
    ~Instrument();

    /// Assignement operator
    Instrument & operator=(const Instrument & source);

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
    /// sets the a description of customizations
    void setCustomizations(const String & customizations);

    /// returns a const reference to the ion source list
    const std::vector<IonSource> & getIonSources() const;
    /// returns a mutable reference to the ion source list
    std::vector<IonSource> & getIonSources();
    /// sets the ion source list
    void setIonSources(const std::vector<IonSource> & ion_sources);

    /// returns a const reference to the mass analyer list
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

#endif // OPENMS_METADATA_INSTRUMENT_H
