// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>

#include <vector>

namespace OpenMS
{

class OPENMS_DLLAPI SimpleSearchEngineAlgorithm :
  public DefaultParamHandler,
  public ProgressLogger
{
  public:
    SimpleSearchEngineAlgorithm(); 

    /// Exit codes
    enum class ExitCodes
    {
      EXECUTION_OK,
      INPUT_FILE_EMPTY,
      UNEXPECTED_RESULT,
      UNKNOWN_ERROR,
      ILLEGAL_PARAMETERS
    };

    /// @brief search spectra against database
    ExitCodes search(const String& in_mzML, 
      const String& in_db, 
      std::vector<ProteinIdentification>& prot_ids,
      std::vector<PeptideIdentification>& pep_ids) const;
  protected:
    void updateMembers_() override;

    /// Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
    struct AnnotatedHit_
    {
      StringView sequence;
      SignedSize peptide_mod_index; ///< enumeration index of the non-RNA peptide modification
      double score = 0; ///< main score
      std::vector<PeptideHit::PeakAnnotation> fragment_annotations;      
      double prefix_fraction = 0; ///< fraction of annotated b-ions
      double suffix_fraction = 0; ///< fraction of annotated y-ions
      double mean_error = 0.0; ///< mean absolute fragment mass error

      static bool hasBetterScore(const AnnotatedHit_& a, const AnnotatedHit_& b)
      {
        if (a.score != b.score) return a.score > b.score;
        // compare the mod_index first, as it is cheaper than the strncmp() of the sequences
        // there doesn't have to be a certain ordering (that makes sense), we just need it to be thread-safe
        if (b.peptide_mod_index != a.peptide_mod_index) return a.peptide_mod_index < b.peptide_mod_index;
        return a.sequence < b.sequence;
      }
    };

    /// @brief filter, deisotope, decharge spectra
    static void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm);

    /// @brief filter and annotate search results
    /// most of the parameters are used to properly add meta data to the id objects
    void postProcessHits_(const PeakMap& exp, 
      std::vector<std::vector<SimpleSearchEngineAlgorithm::AnnotatedHit_> >& annotated_hits, 
      std::vector<ProteinIdentification>& protein_ids, 
      std::vector<PeptideIdentification>& peptide_ids, 
      Size top_hits,
      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
      Size max_variable_mods_per_peptide,
      const StringList& modifications_fixed,
      const StringList& modifications_variable,
      Int peptide_missed_cleavages,
      double precursor_mass_tolerance,
      double fragment_mass_tolerance,
      const String& precursor_mass_tolerance_unit_ppm,
      const String& fragment_mass_tolerance_unit_ppm,
      const Int precursor_min_charge,
      const Int precursor_max_charge,
      const String& enzyme,
      const String& database_name) const;

    double precursor_mass_tolerance_;
    String precursor_mass_tolerance_unit_;

    Size precursor_min_charge_;
    Size precursor_max_charge_;

    IntList precursor_isotopes_;

    double fragment_mass_tolerance_;

    String fragment_mass_tolerance_unit_;

    StringList modifications_fixed_;

    StringList modifications_variable_;

    Size modifications_max_variable_mods_per_peptide_;

    String enzyme_;

    bool decoys_;

    StringList annotate_psm_;

    Size peptide_min_size_;
    Size peptide_max_size_;
    Size peptide_missed_cleavages_;

    String peptide_motif_;

    Size report_top_hits_;
};

} // namespace

