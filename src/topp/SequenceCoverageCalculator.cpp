// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/SequenceCoverage.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SequenceCoverageCalculator SequenceCoverageCalculator

@brief Prints information about idXML files.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref
TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SequenceCoverageCalculator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SequenceCoverageCalculator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSequenceCoverageCalculator : public TOPPBase
{
public:
  TOPPSequenceCoverageCalculator(): TOPPBase("SequenceCoverageCalculator", "Annotates coverage information to idXML files.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_database", "<file>", "", "input file containing the database in FASTA format");
    setValidFormats_("in_database", ListUtils::create<String>("fasta"));
    registerInputFile_("in_peptides", "<file>", "", "input file containing the identified peptides");
    setValidFormats_("in_peptides", ListUtils::create<String>("idXML"), true);
    registerOutputFile_("out", "<file>", "", "Optional text output file. If left out, the output is written to the command line.", false);
    setValidFormats_("out", ListUtils::create<String>("idXML"), true);
  }

  struct CoverageInfo
  {
    double coverage {};                  // fraction of sequence covered by peptides
    Size unique_peptide_sequences {};    // unique peptide sequences (mods ignored)
    Size unique_peptidoforms {};         // unique peptidoforms (sequence + mods; incl. unmodified)
    Size unique_modified_peptidoforms {}; // unique modified peptidoforms
  };

  ExitCodes outputTo_(ostream& os, String out)
  {
    vector<ProteinIdentification> protein_identifications;
    PeptideIdentificationList identifications;
    vector<FASTAFile::FASTAEntry> proteins;
    vector<double> statistics;
    vector<Size> counts;
    vector<Size> peptidoform_counts;
    vector<Size> mod_counts;
    vector<PeptideHit> temp_hits;

    protein_identifications.push_back(ProteinIdentification());
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_name = getStringOption_("in_peptides");
    String database_name = getStringOption_("in_database");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    String document_id;
    FileHandler().loadIdentifications(inputfile_name, protein_identifications, identifications, {FileTypes::IDXML});
    FASTAFile().load(database_name, proteins);

    statistics.resize(proteins.size(), 0.);
    counts.resize(proteins.size(), 0);
    peptidoform_counts.resize(proteins.size(), 0);
    mod_counts.resize(proteins.size(), 0);

    Size identified_spectra = 0;
    unordered_set<String> unique_peptide_sequences_overall;
    unordered_set<String> unique_peptidoforms_overall;
    unordered_set<String> unique_modified_peptidoforms_overall;
    for (const auto& pid : identifications)
    {
      const auto& hits = pid.getHits();
      if (hits.empty())
      {
        continue;
      }
      if (hits.size() > 1)
      {
        OPENMS_LOG_ERROR << "Spectrum with more than one identification found, which is not allowed.\n"
                         << "Use the IDFilter with the -best_hits option to filter for best hits."
                         << "\n";
        return ILLEGAL_PARAMETERS;
      }

      ++identified_spectra;
      const AASequence& seq = hits[0].getSequence();
      unique_peptide_sequences_overall.emplace(seq.toUnmodifiedString());
      unique_peptidoforms_overall.emplace(seq.toString());
      if (seq.isModified())
      {
        unique_modified_peptidoforms_overall.emplace(seq.toString());
      }
    }
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    unordered_map<String, CoverageInfo> prot2cov;
    os << "proteinID\tcoverage (%)\tunique_peptide_sequences\tunique_peptidoforms\tunique_modified_peptidoforms\n";
    for (Size j = 0; j < proteins.size(); ++j)
    {
      const set<String> accession {proteins[j].identifier};
      AASequence protein_seq = AASequence::fromString(proteins[j].sequence);
      std::vector<AASequence> peptides;
      unordered_set<String> unique_peptide_sequences;
      unordered_set<String> unique_peptidoforms;
      unordered_set<String> unique_modified_peptidoforms;

      for (Size i = 0; i < identifications.size(); ++i)
      {
        const auto& hits = identifications[i].getHits();
        if (hits.empty())
        {
          continue;
        }
        temp_hits = PeptideIdentification::getReferencingHits(hits, accession);

        if (temp_hits.size() == 1)
        {
          const AASequence& seq = temp_hits[0].getSequence();
          peptides.push_back(seq);
          unique_peptide_sequences.emplace(seq.toUnmodifiedString());
          unique_peptidoforms.emplace(seq.toString());
          if (seq.isModified())
          {
            unique_modified_peptidoforms.emplace(seq.toString());
          }
        }
      }
      double coverage_percent = SequenceCoverage::getCoverage(protein_seq, peptides);

      statistics[j] = coverage_percent / 100.0;

      counts[j] = unique_peptide_sequences.size();
      peptidoform_counts[j] = unique_peptidoforms.size();
      mod_counts[j] = unique_modified_peptidoforms.size();

      // details for this protein
      if (counts[j] > 0)
      {
        prot2cov[proteins[j].identifier] = {statistics[j], counts[j], peptidoform_counts[j], mod_counts[j]};
        os << proteins[j].identifier << "\t" << statistics[j] * 100 << "\t" << counts[j] << "\t" << peptidoform_counts[j] << "\t"
           << mod_counts[j] << "\n";
      }

      // os << statistics[j] << endl;
    }

    // os << "Sum of coverage is " << accumulate(statistics.begin(), statistics.end(), 0.) << endl;

    // update meta values in protein_identifications
    for (auto& prot_id : protein_identifications)
    {
      for (auto& prot_hit : prot_id.getHits())
      {
        auto it = prot2cov.find(prot_hit.getAccession());
        if (it != prot2cov.end())
        {
          const CoverageInfo& info = it->second;
          prot_hit.setMetaValue("coverage", info.coverage);
          prot_hit.setMetaValue("unique_peptides", info.unique_peptide_sequences);
          prot_hit.setMetaValue("unique_peptidoforms", info.unique_peptidoforms);
          prot_hit.setMetaValue("unique_modified_peptides", info.unique_modified_peptidoforms);
          prot_hit.setMetaValue("unique_modified_peptidoforms", info.unique_modified_peptidoforms);
        }
      }
    }
    FileHandler().storeIdentifications(out, protein_identifications, identifications, {FileTypes::IDXML});

    // safe average helpers to avoid division by zero
    auto safe_avg_double = [](const std::vector<double>& v) -> double { return v.empty() ? 0.0 : (accumulate(v.begin(), v.end(), 0.0) / v.size()); };
    auto safe_avg_size = [](const std::vector<Size>& v) -> double {
      return v.empty() ? 0.0 : (static_cast<double>(accumulate(v.begin(), v.end(), static_cast<Size>(0))) / v.size());
    };

    os << "Average coverage per protein is " << safe_avg_double(statistics) * 100 << "\n";
    os << "Average number of unique peptide sequences per protein is " << safe_avg_size(counts) << "\n";
    os << "Average number of unique peptidoforms per protein is " << safe_avg_size(peptidoform_counts) << "\n";
    os << "Average number of unique modified peptidoforms per protein is " << safe_avg_size(mod_counts) << "\n";
    os << "Number of identified spectra: " << identified_spectra << "\n";
    os << "Number of unique peptide sequences: " << unique_peptide_sequences_overall.size() << "\n";
    os << "Number of unique peptidoforms: " << unique_peptidoforms_overall.size() << "\n";
    os << "Number of unique modified peptidoforms: " << unique_modified_peptidoforms_overall.size() << "\n";

    vector<double>::iterator it = statistics.begin();
    vector<Size>::iterator it2 = counts.begin();
    vector<Size>::iterator it_pf = peptidoform_counts.begin();
    vector<Size>::iterator it3 = mod_counts.begin();
    while (it != statistics.end())
    {
      if (*it == 0.)
      {
        it = statistics.erase(it);
        it2 = counts.erase(it2);
        it_pf = peptidoform_counts.erase(it_pf);
        it3 = mod_counts.erase(it3);
      }
      else
      {
        ++it;
        ++it2;
        ++it_pf;
        ++it3;
      }
    }
    os << "Average coverage per found protein (" << statistics.size() << ") is " << safe_avg_double(statistics) * 100 << "\n";
    os << "Average number of unique peptide sequences per found protein is " << safe_avg_size(counts) << "\n";
    os << "Average number of unique peptidoforms per found protein is " << safe_avg_size(peptidoform_counts) << "\n";
    os << "Average number of unique modified peptidoforms per found protein is " << safe_avg_size(mod_counts) << "\n";

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    String out = getStringOption_("out");
    return outputTo_(getGlobalLogInfo(), out);
  }
};


int main(int argc, const char** argv)
{
  TOPPSequenceCoverageCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
