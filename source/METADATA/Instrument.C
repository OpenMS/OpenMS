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

  Instrument::Instrument(const Instrument & source) :
    MetaInfoInterface(source),
    name_(source.name_),
    vendor_(source.vendor_),
    model_(source.model_),
    customizations_(source.customizations_),
    ion_sources_(source.ion_sources_),
    mass_analyzers_(source.mass_analyzers_),
    ion_detectors_(source.ion_detectors_),
    software_(source.software_),
    ion_optics_(source.ion_optics_)
  {

  }

  Instrument::~Instrument()
  {

  }

  Instrument & Instrument::operator=(const Instrument & source)
  {
    if (&source == this)
      return *this;

    MetaInfoInterface::operator=(source);
    software_ = source.software_;
    name_ = source.name_;
    vendor_ = source.vendor_;
    model_ = source.model_;
    customizations_ = source.customizations_;
    ion_sources_ = source.ion_sources_;
    mass_analyzers_ = source.mass_analyzers_;
    ion_detectors_ = source.ion_detectors_;
    ion_optics_ = source.ion_optics_;

    return *this;
  }

  bool Instrument::operator==(const Instrument & rhs) const
  {
//if (software_ != rhs.software_) cout << "Instrument - " << __LINE__ << endl;
//if (name_ != rhs.name_) cout << "Instrument - " << __LINE__ << endl;
//if (vendor_ != rhs.vendor_) cout << "Instrument - " << __LINE__ << endl;
//if (model_ != rhs.model_) cout << "Instrument - " << __LINE__ << endl;
//if (customizations_ != rhs.customizations_) cout << "Instrument - " << __LINE__ << endl;
//if (ion_sources_ != rhs.ion_sources_) cout << "Instrument - " << __LINE__ << endl;
//if (mass_analyzers_ != rhs.mass_analyzers_) cout << "Instrument - " << __LINE__ << endl;
//if (ion_detectors_ != rhs.ion_detectors_) cout << "Instrument - " << __LINE__ << endl;
//if (ion_optics_ != rhs.ion_optics_) cout << "Instrument - " << __LINE__ << endl;
//if (MetaInfoInterface::operator!=(rhs)) cout << "Instrument - " << __LINE__ << endl;

    return software_ == rhs.software_ &&
           name_ == rhs.name_ &&
           vendor_ == rhs.vendor_ &&
           model_ == rhs.model_ &&
           customizations_ == rhs.customizations_ &&
           ion_sources_ == rhs.ion_sources_ &&
           mass_analyzers_ == rhs.mass_analyzers_ &&
           ion_detectors_ == rhs.ion_detectors_ &&
           ion_optics_ == rhs.ion_optics_ &&
           MetaInfoInterface::operator==(rhs)
    ;
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
