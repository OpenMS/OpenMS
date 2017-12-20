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

#ifndef OPENMS_FORMAT_IDXMLFILE_H
#define OPENMS_FORMAT_IDXMLFILE_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Used to load and store idXML files

    This class is used to load and store documents that implement
    the schema of idXML files.

    A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/.

    One file can contain several ProteinIdentification runs. Each run consists of peptide hits stored in
    PeptideIdentification and (optional) protein hits stored in Identification. Peptide and protein
    hits are connected via a string identifier. We use the search engine and the date as identifier.

    @note This format will eventually be replaced by the HUPO-PSI (mzIdentML and mzQuantML)) AnalysisXML formats!

    @ingroup FileIO
  */
  class OPENMS_DLLAPI IdXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    // both ConsensusXMLFile and FeatureXMLFile use some protected IdXML helper functions to parse identifications without code duplication
    friend class ConsensusXMLFile;
    friend class FeatureXMLFile;

    /// Constructor
    IdXMLFile();

    /**
        @brief Loads the identifications of an idXML file without identifier

        The information is read in and the information is stored in the
        corresponding variables

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids);

    /**
        @brief Loads the identifications of an idXML file

        The information is read in and the information is stored in the
        corresponding variables

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, String& document_id);

    /**
        @brief Stores the data in an idXML file

        The data is read in and stored in the file 'filename'. PeptideHits are sorted by score.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(String filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids, const String& document_id = "");
  

protected:
    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    /// Add data from ProteinGroups to a MetaInfoInterface
    void addProteinGroups_(MetaInfoInterface& meta, const std::vector<ProteinIdentification::ProteinGroup>& groups, const String& group_name, const std::map<String, UInt>& accession_to_id);

    /// Read and store ProteinGroup data
    void getProteinGroups_(std::vector<ProteinIdentification::ProteinGroup>& groups, const String& group_name);

    /**
      * Helper function to create the XML string for the amino acids before and after the peptide position in a protein.
      * Can be reused by e.g. ConsensusXML, FeatureXML to write PeptideHit elements  
      */
    static String createFlankingAAXMLString_(const std::vector<PeptideEvidence> & pes);

    /**
      * Helper function to create the XML string for the position of the peptide in a protein.
      * Can be reused by e.g. ConsensusXML, FeatureXML to write PeptideHit elements  
      */
    static String createPositionXMLString_(const std::vector<PeptideEvidence> & pes);


    /**
      * Helper function to write out fragment annotations as user param fragment_annotation
      */  
    static void writeFragmentAnnotations_(const String & tag_name, std::ostream & os, 
                                          std::vector<PeptideHit::PeakAnnotation> annotations, UInt indent); 

    /**
      * Helper function to parse fragment annotations from string
      */  
    static void parseFragmentAnnotation_(const String& s, std::vector<PeptideHit::PeakAnnotation> & annotations);
    

    /// @name members for loading data
    //@{
    /// Pointer to fill in protein identifications
    std::vector<ProteinIdentification>* prot_ids_;
    /// Pointer to fill in peptide identifications
    std::vector<PeptideIdentification>* pep_ids_;
    /// Pointer to last read object with MetaInfoInterface
    MetaInfoInterface* last_meta_;
    /// Search parameters map (key is the "id")
    std::map<String, ProteinIdentification::SearchParameters> parameters_;
    /// Temporary search parameters variable
    ProteinIdentification::SearchParameters param_;
    /// Temporary id
    String id_;
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
    std::map<String, String> proteinid_to_accession_;
    /// Document identifier
    String* document_id_;
    /// true if a prot id is contained in the current run
    bool prot_id_in_run_;
    //@}
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_IDXMLFILE_H
