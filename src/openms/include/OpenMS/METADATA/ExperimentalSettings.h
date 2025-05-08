// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief Description of the experimental settings

      These settings are valid for the whole experiment.
      See SpectrumSettings for settings which are specific to an MSSpectrum.
      See ChromatogramSettings for settings which are specific to an MSChromatogram.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI ExperimentalSettings :
    public MetaInfoInterface,
    public DocumentIdentifier
  {
public:

    /// Constructor
    ExperimentalSettings() = default;
    /// Copy constructor
    ExperimentalSettings(const ExperimentalSettings &) = default;
    /// Move constructor
    ExperimentalSettings(ExperimentalSettings &&) = default;
    /// Destructor
    ~ExperimentalSettings() override;

    /// Assignment operator
    ExperimentalSettings & operator=(const ExperimentalSettings &) = default;
    /// Move assignment operator
    ExperimentalSettings & operator=(ExperimentalSettings &&) = default;

    /// Equality operator
    bool operator==(const ExperimentalSettings & rhs) const;
    /// Equality operator
    bool operator!=(const ExperimentalSettings & rhs) const;

    /// returns a const reference to the sample description
    const Sample & getSample() const;
    /// returns a mutable reference to the sample description
    Sample & getSample();
    /// sets the sample description
    void setSample(const Sample & sample);

    /// returns a const reference to the source data file
    const std::vector<SourceFile> & getSourceFiles() const;
    /// returns a mutable reference to the source data file
    std::vector<SourceFile> & getSourceFiles();
    /// sets the source data file
    void setSourceFiles(const std::vector<SourceFile> & source_files);

    /// returns a const reference to the list of contact persons
    const std::vector<ContactPerson> & getContacts() const;
    /// returns a mutable reference to the list of contact persons
    std::vector<ContactPerson> & getContacts();
    /// sets the list of contact persons
    void setContacts(const std::vector<ContactPerson> & contacts);

    /// returns a const reference to the MS instrument description
    const Instrument & getInstrument() const;
    /// returns a mutable reference to the MS instrument description
    Instrument & getInstrument();
    /// sets the MS instrument description
    void setInstrument(const Instrument & instrument);

    /// returns a const reference to the description of the HPLC run
    const HPLC & getHPLC() const;
    /// returns a mutable reference to the description of the HPLC run
    HPLC & getHPLC();
    /// sets the description of the HPLC run
    void setHPLC(const HPLC & hplc);

    /// returns the date the experiment was performed
    const DateTime & getDateTime() const;
    /// sets the date the experiment was performed
    void setDateTime(const DateTime & date);

    /// returns the free-text comment
    const String & getComment() const;
    /// sets the free-text comment
    void setComment(const String & comment);

    /// returns a const reference to the protein ProteinIdentification vector
    const std::vector<ProteinIdentification> & getProteinIdentifications() const;
    /// returns a mutable reference to the protein ProteinIdentification vector
    std::vector<ProteinIdentification> & getProteinIdentifications();
    /// sets the protein ProteinIdentification vector
    void setProteinIdentifications(const std::vector<ProteinIdentification> & protein_identifications);

    /// returns fraction identifier
    const String & getFractionIdentifier() const;
    /// sets the fraction identifier
    void setFractionIdentifier(const String & fraction_identifier);

protected:
    Sample sample_;
    std::vector<SourceFile> source_files_;
    std::vector<ContactPerson> contacts_;
    Instrument instrument_;
    HPLC hplc_;
    DateTime datetime_;
    String comment_;
    std::vector<ProteinIdentification> protein_identifications_;
    String fraction_identifier_;
  };

  ///Print the contents to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const ExperimentalSettings & exp);

} // namespace OpenMS

