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

#ifndef OPENMS_FORMAT_XTANDEMXMLFILE_H
#define OPENMS_FORMAT_XTANDEMXMLFILE_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

namespace OpenMS
{
  class String;
  class ProteinIdentification;

  /**
    @brief Used to load XTandemXML files

    This class is used to load documents that implement
    the schema of XTandemXML files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI XTandemXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Default constructor
    XTandemXMLFile();

    /// Destructor
    ~XTandemXMLFile() override;
    /**
      @brief loads data from an X! Tandem XML file

      @param filename the file to be loaded
      @param protein_identification protein identifications belonging to the whole experiment
      @param id_data the identifications with m/z and RT
      @param mod_def_set Fixed and variable modifications defined for the search. May be extended with additional (X! Tandem default) modifications if those are found in the file.

      This class serves to read in an X! Tandem XML file. The information can be
      retrieved via the load function.

      @ingroup FileIO
    */
    void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, ModificationDefinitionsSet& mod_def_set);


protected:

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void characters(const XMLCh* const chars, const XMLSize_t /*length*/) override;

    XTandemXMLFile(const XTandemXMLFile& rhs);

    XTandemXMLFile& operator=(const XTandemXMLFile& rhs);

private:

    ProteinIdentification* protein_identification_;

    // true during "note" element containing protein accession
    bool is_protein_note_;

    // true during "note" element containing spectrum ID
    bool is_spectrum_note_;

    // peptide hits per spectrum
    std::map<UInt, std::vector<PeptideHit> > peptide_hits_;

    // protein hits
    std::vector<ProteinHit> protein_hits_;

    // protein unique IDs (assigned by X! Tandem), to keep track of which proteins were already seen
    std::set<UInt> protein_uids_;

    // accession of the current protein
    String current_protein_;

    // charge of current peptide
    Int current_charge_;

    // X! Tandem ID of current peptide
    UInt current_id_;

    // tag
    String tag_;

    // start position of current peptide in protein sequence
    UInt current_start_;

    // stop position of current peptide in protein sequence
    UInt current_stop_;

    // previous peptide sequence
    String previous_seq_;

    // mapping from X! Tandem ID to spectrum ID
    std::map<UInt, String> spectrum_ids_;

    // modification definitions
    ModificationDefinitionsSet mod_def_set_;

    // modifications used by X! Tandem by default
    ModificationDefinitionsSet default_nterm_mods_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_XTANDEMXMLFILE_H
