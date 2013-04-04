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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ExperimentalSettings.h>

using namespace std;

namespace OpenMS
{

  ExperimentalSettings::ExperimentalSettings() :
    MetaInfoInterface(),
    DocumentIdentifier(),
    sample_(),
    source_files_(),
    contacts_(),
    instrument_(),
    hplc_(),
    datetime_(),
    comment_(),
    protein_identifications_(),
    fraction_identifier_()
  {
  }

  ExperimentalSettings::ExperimentalSettings(const ExperimentalSettings & source) :
    MetaInfoInterface(source),
    DocumentIdentifier(source),
    sample_(source.sample_),
    source_files_(source.source_files_),
    contacts_(source.contacts_),
    instrument_(source.instrument_),
    hplc_(source.hplc_),
    datetime_(source.datetime_),
    comment_(source.comment_),
    protein_identifications_(source.protein_identifications_),
    fraction_identifier_(source.fraction_identifier_)
  {
  }

  ExperimentalSettings::~ExperimentalSettings()
  {
  }

  ExperimentalSettings & ExperimentalSettings::operator=(const ExperimentalSettings & source)
  {
    if (&source == this)
      return *this;

    sample_ = source.sample_;
    source_files_ = source.source_files_;
    contacts_ = source.contacts_;
    instrument_ = source.instrument_;
    hplc_ = source.hplc_;
    datetime_ = source.datetime_;
    comment_ = source.comment_;
    protein_identifications_ = source.protein_identifications_;
    fraction_identifier_ = source.fraction_identifier_;
    MetaInfoInterface::operator=(source);
    DocumentIdentifier::operator=(source);

    return *this;
  }

  bool ExperimentalSettings::operator==(const ExperimentalSettings & rhs) const
  {
    return sample_ == rhs.sample_ &&
           source_files_ == rhs.source_files_ &&
           contacts_ == rhs.contacts_ &&
           instrument_ == rhs.instrument_ &&
           hplc_ == rhs.hplc_ &&
           datetime_ == rhs.datetime_ &&
           protein_identifications_ == rhs.protein_identifications_ &&
           comment_ == rhs.comment_ &&
           fraction_identifier_ == rhs.fraction_identifier_ &&
           MetaInfoInterface::operator==(rhs) &&
           DocumentIdentifier::operator==(rhs);
  }

  bool ExperimentalSettings::operator!=(const ExperimentalSettings & rhs) const
  {
    return !(operator==(rhs));
  }

  const Sample & ExperimentalSettings::getSample() const
  {
    return sample_;
  }

  Sample & ExperimentalSettings::getSample()
  {
    return sample_;
  }

  void ExperimentalSettings::setSample(const Sample & sample)
  {
    sample_ = sample;
  }

  const vector<SourceFile> & ExperimentalSettings::getSourceFiles() const
  {
    return source_files_;
  }

  vector<SourceFile> & ExperimentalSettings::getSourceFiles()
  {
    return source_files_;
  }

  void ExperimentalSettings::setSourceFiles(const vector<SourceFile> & source_file)
  {
    source_files_ = source_file;
  }

  const vector<ContactPerson> & ExperimentalSettings::getContacts() const
  {
    return contacts_;
  }

  vector<ContactPerson> & ExperimentalSettings::getContacts()
  {
    return contacts_;
  }

  void ExperimentalSettings::setContacts(const std::vector<ContactPerson> & contacts)
  {
    contacts_ = contacts;
  }

  const Instrument & ExperimentalSettings::getInstrument() const
  {
    return instrument_;
  }

  Instrument & ExperimentalSettings::getInstrument()
  {
    return instrument_;
  }

  void ExperimentalSettings::setInstrument(const Instrument & instrument)
  {
    instrument_ = instrument;
  }

  const DateTime & ExperimentalSettings::getDateTime() const
  {
    return datetime_;
  }

  void ExperimentalSettings::setDateTime(const DateTime & date)
  {
    datetime_ = date;
  }

  const HPLC & ExperimentalSettings::getHPLC() const
  {
    return hplc_;
  }

  HPLC & ExperimentalSettings::getHPLC()
  {
    return hplc_;
  }

  void ExperimentalSettings::setHPLC(const HPLC & hplc)
  {
    hplc_ = hplc;
  }

  std::ostream & operator<<(std::ostream & os, const ExperimentalSettings & /*exp*/)
  {
    os << "-- EXPERIMENTALSETTINGS BEGIN --" << std::endl;
    os << "-- EXPERIMENTALSETTINGS END --" << std::endl;
    return os;
  }

  const vector<ProteinIdentification> & ExperimentalSettings::getProteinIdentifications() const
  {
    return protein_identifications_;
  }

  vector<ProteinIdentification> & ExperimentalSettings::getProteinIdentifications()
  {
    return protein_identifications_;
  }

  void ExperimentalSettings::setProteinIdentifications(const vector<ProteinIdentification> & protein_identifications)
  {
    protein_identifications_ = protein_identifications;
  }

  const String & ExperimentalSettings::getComment() const
  {
    return comment_;
  }

  void ExperimentalSettings::setComment(const String & comment)
  {
    comment_ = comment;
  }

  const String & ExperimentalSettings::getFractionIdentifier() const
  {
    return fraction_identifier_;
  }

  void ExperimentalSettings::setFractionIdentifier(const String & fraction_identifier)
  {
    fraction_identifier_ = fraction_identifier;
  }

}
