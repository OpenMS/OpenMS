// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_EXPERIMENTALSETTINGS_H
#define OPENMS_METADATA_EXPERIMENTALSETTINGS_H

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
      See SpectrumSettings for settings which are spectrum specific.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI ExperimentalSettings :
    public MetaInfoInterface,
    public DocumentIdentifier
  {
public:
    ///Constructor
    ExperimentalSettings();
    ///Copy constructor
    ExperimentalSettings(const ExperimentalSettings & source);
    ///Destructor
    ~ExperimentalSettings() override;

    ///Assignment operator
    ExperimentalSettings & operator=(const ExperimentalSettings & source);

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

#endif // OPENMS_METADATA_EXPERIMENTALSETTINGS_H
