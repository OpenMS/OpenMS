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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZTABFILE_H
#define OPENMS_FORMAT_MZTABFILE_H

#include <OpenMS/FORMAT/MzTab.h>

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include <vector>
#include <algorithm>

namespace OpenMS
{
  class String;
  class SVOutStream;
/**
    @brief File adapter for MzTab files

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MzTabFile
  {
public:
    ///Default constructor
    MzTabFile();
    ///Destructor
    ~MzTabFile();

    typedef std::map<std::pair<String, String>, std::vector<PeptideHit> > MapAccPepType;

    // store MzTab file
    void store(const String& filename, const MzTab& mz_tab) const;

    // Set store behaviour of optional "reliability" and "uri" columns (default=no)
    void storeProteinReliabilityColumn(bool store);
    void storePeptideReliabilityColumn(bool store);
    void storePSMReliabilityColumn(bool store);
    void storeSmallMoleculeReliabilityColumn(bool store);
    void storeProteinUriColumn(bool store);
    void storePeptideUriColumn(bool store);
    void storePSMUriColumn(bool store);
    void storeSmallMoleculeUriColumn(bool store);
    void storeProteinGoTerms(bool store);

    // load MzTab file
    void load(const String& filename, MzTab& mz_tab);

protected:
    bool store_protein_reliability_;
    bool store_peptide_reliability_;
    bool store_psm_reliability_;
    bool store_smallmolecule_reliability_;
    bool store_protein_uri_;
    bool store_peptide_uri_;
    bool store_psm_uri_;
    bool store_smallmolecule_uri_;
    bool store_protein_goterms_;

    void generateMzTabMetaDataSection_(const MzTabMetaData& map, StringList& sl) const;

    void generateMzTabProteinSection_(const MzTabProteinSectionRows& rows, StringList& sl, const std::vector<String>& optional_columns) const;

    String generateMzTabProteinHeader_(const MzTabProteinSectionRow& reference_row, const Size n_best_search_engine_scores, const std::vector<String>& optional_columns) const;

    String generateMzTabProteinSectionRow_(const MzTabProteinSectionRow& row, const std::vector<String>& optional_columns) const;

    void generateMzTabPeptideSection_(const MzTabPeptideSectionRows& rows, StringList& sl, const std::vector<String>& optional_columns) const;

    String generateMzTabPeptideHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_score, Size assays, Size study_variables, const std::vector<String>& optional_columns) const;

    String generateMzTabPeptideSectionRow_(const MzTabPeptideSectionRow& row, const std::vector<String>& optional_columns) const;

    void generateMzTabPSMSection_(const MzTabPSMSectionRows& rows, StringList& sl, const std::vector<String>& optional_columns) const;

    String generateMzTabPSMHeader_(Size n_search_engine_scores, const std::vector<String>& optional_columns) const;

    String generateMzTabPSMSectionRow_(const MzTabPSMSectionRow& row, const std::vector<String>& optional_columns) const;

    void generateMzTabSmallMoleculeSection_(const MzTabSmallMoleculeSectionRows& map, StringList& sl, const std::vector<String>& optional_columns) const;

    String generateMzTabSmallMoleculeHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_score, Size assays, Size study_variables, const std::vector<String>& optional_columns) const;

    String generateMzTabSmallMoleculeSectionRow_(const MzTabSmallMoleculeSectionRow& row, const std::vector<String>& optional_columns) const;

    // auxiliary functions
    // extract two integers from string (e.g. search_engine_score[1]_ms_run[2] -> 1,2)
    static std::pair<int, int> extractIndexPairsFromBrackets_(const String& s);

    static void sortPSM_(std::vector<PeptideIdentification>::iterator begin, std::vector<PeptideIdentification>::iterator end);

    static void keepFirstPSM_(std::vector<PeptideIdentification>::iterator begin, std::vector<PeptideIdentification>::iterator end);

    /// Extract protein and peptide identifications for each run. maps are assumed empty.
    static void partitionIntoRuns_(const std::vector<PeptideIdentification>& pep_ids,
                                   const std::vector<ProteinIdentification>& pro_ids,
                                   std::map<String, std::vector<PeptideIdentification> >& map_run_to_pepids,
                                   std::map<String, std::vector<ProteinIdentification> >& map_run_to_proids
                                   );


    /// create links from protein to peptides
    static void createProteinToPeptideLinks_(const std::map<String, std::vector<PeptideIdentification> >& map_run_to_pepids, MapAccPepType& map_run_accession_to_pephits);

    /// Extracts, if possible a unique protein accession for a peptide hit in mzTab format. Otherwise NA is returned
    static String extractProteinAccession_(const PeptideHit& peptide_hit);

    /// Extracts, modifications and positions of a peptide hit in mzTab format
    static String extractPeptideModifications_(const PeptideHit& peptide_hit);

    /// Map search engine identifier to CV, param etc.
    static String mapSearchEngineToCvParam_(const String& openms_search_engine_name);

    static String mapSearchEngineScoreToCvParam_(const String& openms_search_engine_name, double score, String score_type);

    static String extractNumPeptides_(const String& common_identifier, const String& protein_accession,
                                      const MapAccPepType& map_run_accession_to_peptides);

    // mzTab definition of distinct
    static String extractNumPeptidesDistinct_(String common_identifier, String protein_accession,
                                              const MapAccPepType& map_run_accession_to_peptides);

    // same as distinct but additional constraint of uniqueness (=maps to exactly one Protein)
    static String extractNumPeptidesUnambiguous_(String common_identifier, String protein_accession,
                                                 const MapAccPepType& map_run_accession_to_peptides);

    static std::map<String, Size> extractNumberOfSubSamples_(const std::map<String, std::vector<ProteinIdentification> >& map_run_to_proids);

    static void writePeptideHeader_(SVOutStream& output, std::map<String, Size> n_sub_samples);

    static void writeProteinHeader_(SVOutStream& output, std::map<String, Size> n_sub_samples);

    static void writeProteinData_(SVOutStream& output,
                                  const ProteinIdentification& prot_id,
                                  Size run_count,
                                  String input_filename,
                                  bool has_coverage,
                                  const MapAccPepType& map_run_accession_to_peptides,
                                  const std::map<String, Size>& map_run_to_num_sub
                                  );

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MZTABFILE_H
