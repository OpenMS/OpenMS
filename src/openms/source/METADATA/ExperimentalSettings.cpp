// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <ostream>

using namespace std;

namespace OpenMS
{

  ExperimentalSettings::~ExperimentalSettings() = default;

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
    os << "-- EXPERIMENTALSETTINGS BEGIN --\n";
    os << "-- EXPERIMENTALSETTINGS END --\n";
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
