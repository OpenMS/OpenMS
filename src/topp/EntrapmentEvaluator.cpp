// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/EntrapmentEvaluation.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

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
own rescoring), then run this tool on the result. It reports, at each q-value threshold, how many PSMs and
peptides were accepted, how many of them match only entrapment proteins, and the resulting estimate of the
true false discovery proportion (FDP):

  FDP = N_entrapment / N_accepted * (1 + T/E)

T/E is the ratio of searchable target to entrapment peptides. It is computed from @p -database when given
(by in-silico digestion, I and L taken as one residue), otherwise taken from @p -te_ratio. A PSM counts as
entrapment only when all its proteins are entrapment proteins and, with @p -database, its sequence is not also a
target peptide of the database with I and L taken as one residue (I and L have the same mass, so such a PSM cannot be
told apart from the target peptide; search engines usually map peptides to proteins without this equivalence).

A nominal FDR that holds gives an FDP at or below the threshold. Use it to validate changes to search and
rescoring: a change that yields more identifications at the same nominal FDR is only an improvement if the FDP
stays controlled.

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
    registerInputFile_("database", "<file>", "", "Searched database (target + entrapment, decoys are skipped). Used to compute T/E (if not given, 'te_ratio' is used) and to count an entrapment PSM whose sequence is also a target peptide with I = L as target.", false);
    setValidFormats_("database", {"fasta"});
    registerDoubleOption_("te_ratio", "<ratio>", 1.0, "Ratio of searchable target to entrapment peptides (used if no 'database' is given).", false);
    setMinFloat_("te_ratio", 0.0);
    registerStringOption_("entrapment_string", "<string>", "ENTRAP_", "Accession prefix of entrapment proteins.", false);
    registerStringOption_("decoy_string", "<string>", "DECOY_", "Decoy string of the database (to skip decoys when computing T/E).", false);
    registerStringOption_("decoy_string_position", "<choice>", "prefix", "Is 'decoy_string' a prefix or a suffix?", false);
    setValidStrings_("decoy_string_position", {"prefix", "suffix"});
    registerDoubleList_("q_value", "<thresholds>", {0.01}, "q-value thresholds at which to report.", false);
    registerStringOption_("class_meta_value", "<name>", "", "Also report per value of this peptide hit meta value.", false);
    registerStringOption_("enzyme", "<enzyme>", "Trypsin", "Enzyme for the in-silico digestion used to compute T/E.", false, true);
    registerIntOption_("missed_cleavages", "<num>", 2, "Missed cleavages for the T/E digestion.", false, true);
    registerIntOption_("min_length", "<num>", 7, "Minimum peptide length for the T/E digestion.", false, true);
    registerIntOption_("max_length", "<num>", 40, "Maximum peptide length for the T/E digestion.", false, true);
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList in = getStringList_("in");
    const std::string out = getStringOption_("out");
    const std::string entrapment_prefix = getStringOption_("entrapment_string");
    const DoubleList thresholds = getDoubleList_("q_value");
    const std::string class_meta_value = getStringOption_("class_meta_value");

    double te_ratio = getDoubleOption_("te_ratio");
    const std::string database = getStringOption_("database");
    std::unordered_set<std::string> target_peptides;
    if (!database.empty())
    {
      std::vector<FASTAFile::FASTAEntry> db;
      FASTAFile().load(database, db);
      const std::string decoy_string = getStringOption_("decoy_string");
      const bool decoy_is_prefix = getStringOption_("decoy_string_position") == "prefix";
      const std::string enzyme = getStringOption_("enzyme");
      te_ratio = EntrapmentEvaluation::estimateTERatio(db, entrapment_prefix, decoy_string, decoy_is_prefix, enzyme,
                                                       getIntOption_("missed_cleavages"), getIntOption_("min_length"),
                                                       getIntOption_("max_length"));
      target_peptides = EntrapmentEvaluation::targetPeptides(db, entrapment_prefix, decoy_string, decoy_is_prefix, enzyme,
                                                             getIntOption_("missed_cleavages"), getIntOption_("min_length"),
                                                             getIntOption_("max_length"));
      OPENMS_LOG_INFO << "T/E from '" << database << "': " << te_ratio << '\n';
    }

    std::vector<std::pair<std::string, PeptideIdentificationList>> runs;
    PeptideIdentificationList all;
    // label each file by its name, or by its path if two inputs share a name
    std::map<std::string, Size> name_count;
    for (const auto& file : in) ++name_count[std::filesystem::path(file).filename().string()];
    for (const auto& file : in)
    {
      std::vector<ProteinIdentification> proteins;
      PeptideIdentificationList peptides;
      FileHandler().loadIdentifications(file, proteins, peptides, {FileTypes::IDXML, FileTypes::IDPARQUET});
      for (auto& pid : peptides) pid.sort();
      all.insert(all.end(), peptides.begin(), peptides.end());
      const std::string name = std::filesystem::path(file).filename().string();
      runs.emplace_back(name_count[name] > 1 ? std::string(file) : name, std::move(peptides));
    }
    if (runs.size() > 1) runs.emplace_back("all files", std::move(all));

    std::ostringstream table;
    table << "file\tclass\tq_value\taccepted_psms\tentrapment_psms\tpsm_fdp\taccepted_peptides\tentrapment_peptides\tpeptide_fdp\tte_ratio\n";
    for (const auto& [name, peptides] : runs)
    {
      for (double q : thresholds)
      {
        const auto results = EntrapmentEvaluation::evaluate(peptides, q, entrapment_prefix, te_ratio, class_meta_value,
                                                            database.empty() ? nullptr : &target_peptides);
        for (const auto& [cls, r] : results)
        {
          table << name << '\t' << cls << '\t' << q << '\t' << r.accepted_psms << '\t' << r.entrapment_psms << '\t'
                << std::setprecision(4) << r.psmFDP() << '\t' << r.accepted_peptides << '\t' << r.entrapment_peptides << '\t'
                << r.peptideFDP() << '\t' << r.te_ratio << std::setprecision(6) << '\n';
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
