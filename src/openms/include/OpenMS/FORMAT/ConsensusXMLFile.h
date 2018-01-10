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
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CONSENSUSXMLFILE_H
#define OPENMS_FORMAT_CONSENSUSXMLFILE_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

namespace OpenMS
{
  /**
    @brief This class provides Input functionality for ConsensusMaps and Output functionality for
    alignments and quantitation.

    This class can be used to load the content of a consensusXML file into a ConsensusMap
    or to save the content of an ConsensusMap object into an XML file.

    A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/.

  @todo Take care that unique ids are assigned properly by TOPP tools before calling ConsensusXMLFile::store().  There will be a message on LOG_INFO but we will make no attempt to fix the problem in this class.  (all developers)

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ConsensusXMLFile :
    public Internal::XMLHandler,
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    ConsensusXMLFile();
    ///Destructor
    ~ConsensusXMLFile() override;


    /**
    @brief Loads a consensus map from file and calls updateRanges

    @exception Exception::FileNotFound is thrown if the file could not be opened
    @exception Exception::ParseError is thrown if an error occurs during parsing
    @exception Exception::MissingInformation is thrown if source files are missing/duplicated or map-IDs are referencing non-existing maps
    */
    void load(const String& filename, ConsensusMap& map);

    /**
    @brief Stores a consensus map to file

    @exception Exception::UnableToCreateFile is thrown if the file name is not writable
    @exception Exception::IllegalArgument is thrown if the consensus map is not valid
    @exception Exception::MissingInformation is thrown if source files are missing/duplicated or map-IDs are referencing non-existing maps
    */
    void store(const String& filename, const ConsensusMap& consensus_map);

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

protected:

    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    // Docu in base class
    void characters(const XMLCh* const chars, const XMLSize_t length) override;


    /// Writes a peptide identification to a stream (for assigned/unassigned peptide identifications)
    void writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name, UInt indentation_level);


    /// Options that can be set
    PeakFileOptions options_;

    ///@name Temporary variables for parsing
    //@{
    ConsensusMap* consensus_map_;
    ConsensusFeature act_cons_element_;
    DPosition<2> pos_;
    double it_;
    //@}

    /// Pointer to last read object as a MetaInfoInterface, or null.
    MetaInfoInterface* last_meta_;
    /// Temporary protein ProteinIdentification
    ProteinIdentification prot_id_;
    /// Temporary peptide ProteinIdentification
    PeptideIdentification pep_id_;
    /// Temporary protein hit
    ProteinHit prot_hit_;
    /// Temporary peptide hit
    PeptideHit pep_hit_;
    /// Temporary peptide evidences
    std::vector<PeptideEvidence> peptide_evidences_;
    /// Map from protein id to accession
    Map<String, String> proteinid_to_accession_;
    /// Map from search identifier concatenated with protein accession to id
    Map<String, Size> accession_to_id_;
    /// Map from identification run identifier to file xs:id (for linking peptide identifications to the corresponding run)
    Map<String, String> identifier_id_;
    /// Map from file xs:id to identification run identifier (for linking peptide identifications to the corresponding run)
    Map<String, String> id_identifier_;
    /// Temporary search parameters file
    ProteinIdentification::SearchParameters search_param_;

    UInt progress_;

  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_CONSENSUSXMLFILE_H

