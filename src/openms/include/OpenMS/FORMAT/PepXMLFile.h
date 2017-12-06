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
// $Maintainer: Chris Bielow, Hendrik Weisser $
// $Authors: Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEPXMLFILE_H
#define OPENMS_FORMAT_PEPXMLFILE_H

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

#include <vector>
#include <map>
#include <set>


namespace OpenMS
{
  /**
    @brief Used to load and store PepXML files

    This class is used to load and store documents that implement the schema of PepXML files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI PepXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Constructor
    PepXMLFile();

    /// Destructor
    ~PepXMLFile() override;

    /**
        @brief Loads peptide sequences with modifications out of a PepXML file

        @param filename PepXML file to load
        @param proteins Protein identification output
        @param peptides Peptide identification output
        @param experiment_name Experiment file name, which is used to extract the corresponding search results from the PepXML file.
        @param lookup Helper for looking up retention times (PepXML may contain only scan numbers).

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename,
              std::vector<ProteinIdentification>& proteins,
              std::vector<PeptideIdentification>& peptides,
              const String& experiment_name,
              const SpectrumMetaDataLookup& lookup);

    /**
        @brief @a load function with empty defaults for some parameters (see above)

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename,
              std::vector<ProteinIdentification>& proteins,
              std::vector<PeptideIdentification>& peptides,
              const String& experiment_name = "");

    /**
        @brief Stores idXML as PepXML file

        @exception Exception::UnableToCreateFile is thrown if the file could not be opened for writing
    */
    void store(const String& filename, std::vector<ProteinIdentification>& protein_ids, 
               std::vector<PeptideIdentification>& peptide_ids, const String& mz_file = "",
               const String& mz_name = "", bool peptideprophet_analyzed = false);

    /**
        @brief Whether we should keep the native spectrum name of the pepXML

        @note This will lead to a "pepxml_spectrum_name" meta value being added
        to each PeptideIdentification containing the original name of the
        spectrum in TPP format.
    */
    void keepNativeSpectrumName(bool keep) 
    {
      keep_native_name_ = keep;
    }

protected:

    /// Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    /// Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

private:

    /// Fill @p scan_map_
    void makeScanMap_();

    /// Read RT, m/z, charge information from attributes of "spectrum_query"
    void readRTMZCharge_(const xercesc::Attributes& attributes);

    /**
        @brief find modification name given a modified AA mass

        Matches a mass of a modified AA to a mod in our modification db
        For ambiguous mods, the first (arbitrary) is returned
        If no mod is found an error is issued and the return string is empty
        @note A duplicate of this function is also used in ProtXMLFile

        @param mass Modified AA's mass
        @param origin AA one letter code
        @param modification_description [out] Name of the modification, e.g. 'Carboxymethyl (C)'
    */
    void matchModification_(const double mass, const String& origin, String& modification_description);

    struct AminoAcidModification
    {
      String aminoacid;
      String massdiff;
      double mass;
      bool variable;
      String description;
      String terminus;

      AminoAcidModification() :
        mass(0),
        variable(false)
      {
      }

      AminoAcidModification(const AminoAcidModification& rhs) :
        aminoacid(rhs.aminoacid),
        massdiff(rhs.massdiff),
        mass(rhs.mass),
        variable(rhs.variable),
        description(rhs.description),
        terminus(rhs.terminus)
      {
      }

      virtual ~AminoAcidModification()
      {
      }

      AminoAcidModification& operator=(const AminoAcidModification& rhs)
      {
        if (this != &rhs)
        {
          aminoacid = rhs.aminoacid;
          massdiff = rhs.massdiff;
          mass = rhs.mass;
          variable = rhs.variable;
          description = rhs.description;
          terminus = rhs.terminus;
        }
        return *this;
      }

    };

    /// Pointer to the list of identified proteins
    std::vector<ProteinIdentification>* proteins_;

    /// Pointer to the list of identified peptides
    std::vector<PeptideIdentification>* peptides_;

    /// Pointer to wrapper for looking up spectrum meta data
    const SpectrumMetaDataLookup* lookup_;

    /// Name of the associated experiment (filename of the data file, extension will be removed)
    String exp_name_;

    /// Set name of search engine
    String search_engine_;

    /// Several optional attributes of spectrum_query
    String native_spectrum_name_;
    String experiment_label_;
    String swath_assay_;
    String status_;

    /// Get RT and m/z for peptide ID from precursor scan (should only matter for RT)?
    bool use_precursor_data_;

    /// Mapping between scan number in the pepXML file and index in the corresponding MSExperiment
    std::map<Size, Size> scan_map_;

    /// Hydrogen data (for mass types)
    Element hydrogen_;

    /// Are we currently in an "analysis_summary" element (should be skipped)?
    bool analysis_summary_;

    /// Whether we should keep the native spectrum name of the pepXML
    bool keep_native_name_;

    /// Are we currently in an "search_score_summary" element (should be skipped)?
    bool search_score_summary_;

    /// Are we currently in an "search_summary" element (should be skipped)?
    bool search_summary_;

    /// Do current entries belong to the experiment of interest (for pepXML files that bundle results from different experiments)?
    bool wrong_experiment_;

    /// Have we seen the experiment of interest at all?
    bool seen_experiment_;

    /// Have we checked the "base_name" attribute in the "msms_run_summary" element?
    bool checked_base_name_;

    /// current base name
    String current_base_name_;

    /// References to currently active ProteinIdentifications
    std::vector<std::vector<ProteinIdentification>::iterator> current_proteins_;

    /// Search parameters of the current identification run
    ProteinIdentification::SearchParameters params_;

    /// Enzyme name associated with the current identification run
    String enzyme_;

    /// PeptideIdentification instance currently being processed
    PeptideIdentification current_peptide_;

    /// Analysis result instance currently being processed
    PeptideHit::PepXMLAnalysisResult current_analysis_result_;

    /// PeptideHit instance currently being processed
    PeptideHit peptide_hit_;

    /// Sequence of the current peptide hit
    String current_sequence_;

    /// RT and m/z of current PeptideIdentification
    double rt_, mz_;

    /// Precursor ion charge
    Int charge_;

    /// ID of current search result
    UInt search_id_;

    /// Identifier linking PeptideIdentifications and ProteinIdentifications
    String prot_id_;

    /// Date the pepXML file was generated
    DateTime date_;

    /// Mass of a hydrogen atom (monoisotopic/average depending on case)
    double hydrogen_mass_;

    /// The modifications of the current peptide hit (position is 1-based)
    std::vector<std::pair<String, Size> > current_modifications_;

    /// Fixed aminoacid modifications
    std::vector<AminoAcidModification> fixed_modifications_;

    /// Variable aminoacid modifications
    std::vector<AminoAcidModification> variable_modifications_;

    //@}

    static const double mod_tol_;
    static const double xtandem_artificial_mod_tol_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEPXMLFILE_H
