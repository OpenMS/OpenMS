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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_OMSSAXMLFILE_H
#define OPENMS_FORMAT_OMSSAXMLFILE_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

#include <vector>

namespace OpenMS
{
  class String;
  class ModificationsDB;

  /**
    @brief Used to load OMSSAXML files

    This class is used to load documents that implement
    the schema of OMSSAXML files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI OMSSAXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Default constructor
    OMSSAXMLFile();

    /// Destructor
    ~OMSSAXMLFile() override;
    /**
      @brief loads data from a OMSSAXML file

      @param filename The file to be loaded
      @param protein_identification Protein identifications belonging to the whole experiment
      @param id_data The identifications with m/z and RT
      @param load_proteins If this flag is set to false, the protein identifications are not loaded
      @param load_empty_hits Many spectra will not return a hit. Report empty peptide identifications?

      This class serves to read in a OMSSAXML file. The information can be
      retrieved via the load function.

          @exception FileNotFound
          @exception ParseError

      @ingroup FileIO
    */
    void load(const String& filename,
              ProteinIdentification& protein_identification,
              std::vector<PeptideIdentification>& id_data,
              bool load_proteins = true,
              bool load_empty_hits = true);

    /// sets the valid modifications
    void setModificationDefinitionsSet(const ModificationDefinitionsSet& rhs);

protected:
    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    // Docu in base class
    void characters(const XMLCh* const chars, const XMLSize_t /*length*/) override;

private:

    OMSSAXMLFile(const OMSSAXMLFile& rhs);

    OMSSAXMLFile& operator=(const OMSSAXMLFile& rhs);

    /// reads the mapping file needed for modifications
    void readMappingFile_();

    /// the identifications (storing the peptide hits)
    std::vector<PeptideIdentification>* peptide_identifications_;

    ProteinHit actual_protein_hit_;

    PeptideHit actual_peptide_hit_;

    PeptideEvidence actual_peptide_evidence_;

    std::vector<PeptideEvidence> actual_peptide_evidences_;

    PeptideIdentification actual_peptide_id_;

    ProteinIdentification actual_protein_id_;

    String tag_;

    /// site of the actual modification (simple position in the peptide)
    UInt actual_mod_site_;

    /// type of the modification
    String actual_mod_type_;

    /// modifications of the peptide defined by site and type
    std::vector<std::pair<UInt, String> > modifications_;

    /// should protein hits be read from the file?
    bool load_proteins_;

    /// should empty peptide identifications be loaded or skipped?
    bool load_empty_hits_;

    /// modifications mapping file from OMSSA mod num to UniMod accession
    Map<UInt, std::vector<ResidueModification> > mods_map_;

    /// modification mapping reverse, from the modification to the mod_num
    Map<String, UInt> mods_to_num_;

    /// modification definitions set of the search, needed to annotate fixed modifications
    ModificationDefinitionsSet mod_def_set_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_OMSSAXMLFILE_H
