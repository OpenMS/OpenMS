// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/EntrapmentEvaluation.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <unordered_set>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_EntrapmentEvaluator EntrapmentEvaluator

@brief Measures the false discovery proportion of identifications with an entrapment database.

Search a database that contains entrapment sequences (for example made by DecoyDatabase with
@p -entrapment_ratio 1), estimate the FDR as usual (FalseDiscoveryRate, PercolatorAdapter or a search engine's
own rescoring), then run this tool on the result. At each q-value threshold it reports how many discoveries were
accepted, how many of them match only entrapment proteins, and estimates of the true false discovery proportion
(FDP), following Wen et al., "Assessment of false discovery rate control in tandem mass spectrometry analysis using
entrapment", Nat. Methods (2025). With N accepted discoveries, N_E of them entrapment:

  - lower_bound_fdp = N_E / N: a lower bound. Above the nominal FDR it shows that FDR control fails; below it, it
    shows nothing.
  - combined_fdp = N_E / N * (1 + T/E): an estimated upper bound. Below the nominal FDR it is evidence of control.
  - paired_fdp = (N_E + N_{E>=s>T} + 2 N_{E>T>=s}) / N: a tighter estimated upper bound, at the peptide level, for a
    database with one entrapment protein per target as DecoyDatabase -entrapment_ratio 1 writes it (needs
    @p -database). Every entrapment peptide is paired with the target peptide it was shuffled from;
    N_{E>=s>T} counts accepted entrapment peptides whose paired target is not accepted, N_{E>T>=s} those whose paired
    target is accepted but does not rank above them (q-value, then @p -score_meta_value; a tie counts against the FDP).

Two levels are reported: "psm" (the top hit of each spectrum) and "peptide" (distinct unmodified sequences, I = L).
The paper evaluates peptides, not PSMs. A peptide is accepted at its peptide-level q-value when the PSMs carry one
(FalseDiscoveryRate -peptide true writes it as "peptide q-value"), otherwise at the best q-value of its PSMs, which is
not a peptide-level FDR (the tool warns).

T/E is the ratio of searchable target to entrapment peptides. It is computed from @p -database when given
(by in-silico digestion, I and L taken as one residue), otherwise taken from @p -te_ratio. A discovery counts as
entrapment only when all its proteins are entrapment proteins and, with @p -database, its sequence is not also a
target peptide of the database with I and L taken as one residue (I and L have the same mass, so such a PSM cannot be
told apart from the target peptide; search engines usually map peptides to proteins without this equivalence).

Use it to validate changes to search and rescoring: a change that yields more identifications at the same nominal
FDR is only an improvement if the FDP stays controlled.

With @p -class_meta_value the result is also broken down by the value of a peptide hit meta value, for
example to check the FDR within a modification class.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_EntrapmentEvaluator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_EntrapmentEvaluator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPEntrapmentEvaluator :
  public TOPPBase
{
public:
  TOPPEntrapmentEvaluator() :
    TOPPBase("EntrapmentEvaluator", "Measures the false discovery proportion of identifications with an entrapment database.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Identifications with q-values (as main score or 'q-value' meta value). Each file is evaluated on its own and all files together.");
    setValidFormats_("in", {"idXML", "idparquet"});
    registerOutputFile_("out", "<file>", "", "Optional tab-separated output (one row per file name, class and threshold).", false);
    setValidFormats_("out", {"tsv"});
    registerInputFile_("database", "<file>", "", "Searched database (target + entrapment, decoys are skipped). Used to compute T/E (if not given, 'te_ratio' is used), to count an entrapment PSM whose sequence is also a target peptide with I = L as target, and to pair entrapment with target peptides for 'paired_fdp'.", false);
    setValidFormats_("database", {"fasta"});
    registerDoubleOption_("te_ratio", "<ratio>", 1.0, "Ratio of searchable target to entrapment peptides (used if no 'database' is given).", false);
    setMinFloat_("te_ratio", 0.0);
    registerStringOption_("entrapment_string", "<string>", "ENTRAP_", "Accession prefix of entrapment proteins.", false);
    registerStringOption_("decoy_string", "<string>", "DECOY_", "Decoy string of the database. Decoy proteins are skipped when computing T/E, and a hit counts as entrapment if all its non-decoy proteins are entrapment.", false);
    registerStringOption_("decoy_string_position", "<choice>", "prefix", "Is 'decoy_string' a prefix or a suffix?", false);
    setValidStrings_("decoy_string_position", {"prefix", "suffix"});
    registerDoubleList_("q_value", "<thresholds>", {0.01}, "q-value thresholds at which to report.", false);
    registerStringOption_("class_meta_value", "<name>", "", "Also report per value of this peptide hit meta value.", false);
    registerStringOption_("enzyme", "<enzyme>", "Trypsin", "Enzyme for the in-silico digestion used to compute T/E.", false, true);
    registerIntOption_("missed_cleavages", "<num>", 2, "Missed cleavages for the T/E digestion.", false, true);
    setMinInt_("missed_cleavages", 0);
    registerIntOption_("min_length", "<num>", 7, "Minimum peptide length for the T/E digestion.", false, true);
    setMinInt_("min_length", 1);
    registerIntOption_("max_length", "<num>", 40, "Maximum peptide length for the T/E digestion.", false, true);
    setMinInt_("max_length", 1);
    registerStringOption_("score_meta_value", "<name>", "", "Peptide hit meta value (higher is better) that ranks peptides with equal q-values against their paired target for 'paired_fdp'. Empty: ties count against the FDP.", false, true);
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList in = getStringList_("in");
    const std::string out = getStringOption_("out");
    const std::string entrapment_prefix = getStringOption_("entrapment_string");
    const DoubleList thresholds = getDoubleList_("q_value");
    const std::string class_meta_value = getStringOption_("class_meta_value");

    const std::string score_meta_value = getStringOption_("score_meta_value");
    double te_ratio = getDoubleOption_("te_ratio");
    const std::string database = getStringOption_("database");
    const std::string decoy_string = getStringOption_("decoy_string");
    const bool decoy_is_prefix = getStringOption_("decoy_string_position") == "prefix";
    std::unordered_set<std::string> target_peptides;
    EntrapmentEvaluation::PeptidePairs pairs;
    if (!database.empty())
    {
      std::vector<FASTAFile::FASTAEntry> db;
      FASTAFile().load(database, db);
      const std::string enzyme = getStringOption_("enzyme");
      const Size missed_cleavages = getIntOption_("missed_cleavages");
      const Size min_length = getIntOption_("min_length");
      const Size max_length = getIntOption_("max_length");
      te_ratio = EntrapmentEvaluation::estimateTERatio(db, entrapment_prefix, decoy_string, decoy_is_prefix, enzyme,
                                                       missed_cleavages, min_length, max_length);
      target_peptides = EntrapmentEvaluation::targetPeptides(db, entrapment_prefix, decoy_string, decoy_is_prefix, enzyme,
                                                             missed_cleavages, min_length, max_length);
      Size mismatched = 0;
      pairs = EntrapmentEvaluation::pairPeptides(db, entrapment_prefix, decoy_string, decoy_is_prefix, enzyme, missed_cleavages, min_length, max_length, &mismatched);
      OPENMS_LOG_INFO << "T/E from '" << database << "': " << te_ratio << '\n';
      OPENMS_LOG_INFO << "Entrapment peptides paired with their target peptide: " << pairs.size() << '\n';
      if (mismatched > 0)
      {
        OPENMS_LOG_WARN << "Warning: " << mismatched << " entrapment peptides or proteins do not match their target in "
                        << "composition and were not paired. 'paired_fdp' expects a database made by DecoyDatabase -entrapment_ratio 1." << endl;
      }
      if (pairs.empty())
      {
        OPENMS_LOG_WARN << "Warning: no entrapment peptide could be paired with a target peptide, so 'paired_fdp' is not "
                        << "reported. It needs one entrapment protein per target, named '" << entrapment_prefix
                        << "<target accession>' (DecoyDatabase -entrapment_ratio 1)." << endl;
      }
    }
    const EntrapmentEvaluation::PeptidePairs* pairs_ptr = pairs.empty() ? nullptr : &pairs;

    std::vector<std::pair<std::string, PeptideIdentificationList>> runs;
    PeptideIdentificationList all;
    // label each file by its name, or by its path if two inputs share a name
    std::map<std::string, Size> name_count;
    for (const auto& file : in) ++name_count[std::filesystem::path(file).filename().string()];
    bool has_peptide_q = false;
    for (const auto& file : in)
    {
      std::vector<ProteinIdentification> proteins;
      PeptideIdentificationList peptides;
      FileHandler().loadIdentifications(file, proteins, peptides, {FileTypes::IDXML, FileTypes::IDPARQUET});
      for (auto& pid : peptides)
      {
        pid.sort();
        if (!has_peptide_q && !pid.getHits().empty()) has_peptide_q = pid.getHits()[0].metaValueExists(Constants::UserParam::PEPTIDE_Q_VALUE);
      }
      all.insert(all.end(), peptides.begin(), peptides.end());
      const std::string name = std::filesystem::path(file).filename().string();
      runs.emplace_back(name_count[name] > 1 ? std::string(file) : name, std::move(peptides));
    }
    if (runs.size() > 1) runs.emplace_back("all files", std::move(all));
    if (!has_peptide_q)
    {
      OPENMS_LOG_WARN << "Warning: the hits carry no peptide-level q-value ('" << Constants::UserParam::PEPTIDE_Q_VALUE
                      << "', e.g. from FalseDiscoveryRate -peptide true). Peptide-level rows accept a peptide at the best "
                      << "q-value of its PSMs, which is not a peptide-level FDR." << endl;
    }

    auto fmt = [](double v)
    {
      if (std::isnan(v)) return std::string("NA");
      std::ostringstream os;
      os << std::setprecision(4) << v;
      return os.str();
    };
    std::ostringstream table;
    table << "file\tlevel\tclass\tq_value\taccepted\tentrapment\tlower_bound_fdp\tcombined_fdp\tpaired_fdp"
             "\tentrapment_target_not_accepted\tentrapment_above_target\tte_ratio\n";
    const std::vector<std::pair<std::string, EntrapmentEvaluation::Level>> levels = {
      {"psm", EntrapmentEvaluation::Level::PSM}, {"peptide", EntrapmentEvaluation::Level::PEPTIDE}};
    for (const auto& [name, peptides] : runs)
    {
      for (const auto& [level_name, level] : levels)
      {
        for (double q : thresholds)
        {
          const auto results = EntrapmentEvaluation::evaluate(peptides, q, entrapment_prefix, te_ratio, level, class_meta_value,
                                                              database.empty() ? nullptr : &target_peptides, pairs_ptr,
                                                              score_meta_value, decoy_string, decoy_is_prefix);
          for (const auto& [cls, r] : results)
          {
            table << name << '\t' << level_name << '\t' << cls << '\t' << q << '\t' << r.accepted << '\t' << r.entrapment << '\t'
                  << fmt(r.lowerBoundFDP()) << '\t' << fmt(r.combinedFDP()) << '\t' << fmt(r.pairedFDP()) << '\t'
                  << (r.paired ? std::to_string(r.entrapment_target_not_accepted) : std::string("NA")) << '\t'
                  << (r.paired ? std::to_string(r.entrapment_above_target) : std::string("NA")) << '\t'
                  << fmt(r.te_ratio) << '\n';
          }
        }
      }
    }

    OPENMS_LOG_INFO << table.str();
    if (!out.empty())
    {
      std::ofstream os(out);
      if (!os) throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out);
      os << table.str();
    }
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPEntrapmentEvaluator tool;
  return tool.main(argc, argv);
}

/// @endcond
