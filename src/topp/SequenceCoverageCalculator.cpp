// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <unordered_map>
#include <numeric>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SequenceCoverageCalculator SequenceCoverageCalculator

@brief Prints information about idXML files.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SequenceCoverageCalculator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SequenceCoverageCalculator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSequenceCoverageCalculator :
  public TOPPBase
{
public:
  TOPPSequenceCoverageCalculator() :
    TOPPBase("SequenceCoverageCalculator", "Annotates coverage information to idXML files.")
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

  void getStartAndEndIndex(const String& sequence, const String& substring, pair<Size, Size>& indices)
  {
    indices.first = 0;
    indices.second = 0;

    if (sequence.hasSubstring(substring))
    {
      for (Size i = 0; i <= sequence.size() - substring.size(); ++i)
      {
        Size temp_index = i;
        Size temp_count = 0;
        while (temp_index < sequence.size()
              && temp_count < substring.size()
              && sequence.at(temp_index) == substring.at(temp_index - i))
        {
          ++temp_index;
          ++temp_count;
        }
        if (temp_count == substring.size())
        {
          indices.first = i;
          indices.second = temp_index;
          i = sequence.size();
        }
      }
    }
  }

  struct CoverageInfo
  {
    double coverage{};  // fraction of sequence covered by peptides
    Size count{};       // number of unique peptides
    Size mod_count{};   // number of unique modified peptides
  };

  ExitCodes outputTo_(ostream& os, String out)
  {
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    vector<FASTAFile::FASTAEntry> proteins;
    vector<double> statistics;
    vector<Size> counts;
    vector<Size> mod_counts;
    vector<PeptideHit> temp_hits;
    vector<Size> coverage;
    Size spectrum_count = 0;
    unordered_map<String, Size> unique_peptides;
    unordered_map<String, Size> temp_unique_peptides;
    unordered_map<String, Size> temp_modified_unique_peptides;

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
    mod_counts.resize(proteins.size(), 0);
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    unordered_map<String, CoverageInfo> prot2cov;
    os << "proteinID\tcoverage (%)\tunique hits\n";
    for (Size j = 0; j < proteins.size(); ++j)
    {
      coverage.clear();
      coverage.resize(proteins[j].sequence.size(), 0);
      temp_unique_peptides.clear();
      temp_modified_unique_peptides.clear();

      for (Size i = 0; i < identifications.size(); ++i)
      {
        if (!identifications[i].empty())
        {
          if (identifications[i].getHits().size() > 1)
          {
            OPENMS_LOG_ERROR << "Spectrum with more than one identification found, which is not allowed.\n"
                      << "Use the IDFilter with the -best_hits option to filter for best hits." << endl;
            return ILLEGAL_PARAMETERS;
          }

          set<String> accession;
          accession.insert(proteins[j].identifier);
          temp_hits = PeptideIdentification::getReferencingHits(identifications[i].getHits(), accession);

          if (temp_hits.size() == 1)
          {
            pair<Size, Size> indices;
            getStartAndEndIndex(proteins[j].sequence, temp_hits[0].getSequence().toUnmodifiedString(), indices);
            for (Size k = indices.first; k < indices.second; ++k)
            {
              coverage[k] = 1;
            }
            if (indices.first != indices.second)
            {
              // os <<  temp_hits[0].getSequence().toUnmodifiedString() << endl;
            }
            ++spectrum_count;
            if (unique_peptides.find(temp_hits[0].getSequence().toString()) == unique_peptides.end())
            {
              unique_peptides.insert(make_pair(temp_hits[0].getSequence().toString(), 0));
            }
            if (temp_unique_peptides.find(temp_hits[0].getSequence().toUnmodifiedString()) == temp_unique_peptides.end())
            {
              temp_unique_peptides.insert(make_pair(temp_hits[0].getSequence().toUnmodifiedString(), 0));
            }
            if (temp_modified_unique_peptides.find(temp_hits[0].getSequence().toUnmodifiedString()) == temp_modified_unique_peptides.end())
            {
              temp_modified_unique_peptides.insert(make_pair(temp_hits[0].getSequence().toString(), 0));
            }
          }
        }
      }
/* << proteins[j].sequence << endl;
                for (Size k = 0; k < coverage.size(); ++k)
                {
                    os << coverage[k];
                }
                os << endl;
*/
      // statistics[j] = make_pair(,
      // accumulate(coverage.begin(), coverage.end(), 0) / proteins[j].sequence.size());
      statistics[j] = ((double) accumulate(coverage.begin(), coverage.end(), Size(0))) / proteins[j].sequence.size();
      counts[j] = temp_unique_peptides.size();
      mod_counts[j] = temp_modified_unique_peptides.size();

      // details for this protein
      if (counts[j] > 0)
      {
        prot2cov[proteins[j].identifier] = { statistics[j], counts[j], mod_counts[j] };
        os << proteins[j].identifier << "\t" << statistics[j] * 100 << "\t" << counts[j] << "\n";
      }

// os << statistics[j] << endl;
    }

// os << "Sum of coverage is " << accumulate(statistics.begin(), statistics.end(), 0.) << endl;

    // update meta values in protein_identifications
    for (auto& prot_id : protein_identifications)
    {
      for (auto& prot_hit : prot_id.getHits())
      {
        if (prot2cov.find(prot_hit.getAccession()) != prot2cov.end())
        {
          prot_hit.setMetaValue("coverage", prot2cov[prot_hit.getAccession()].coverage);
          prot_hit.setMetaValue("unique_peptides", prot2cov[prot_hit.getAccession()].count);
          prot_hit.setMetaValue("unique_modified_peptides", prot2cov[prot_hit.getAccession()].mod_count);
        }
      }
    }
    FileHandler().storeIdentifications(out, protein_identifications, identifications, {FileTypes::IDXML});

    os << "Average coverage per protein is " << (accumulate(statistics.begin(), statistics.end(), 0.) / statistics.size()) << endl;
    os << "Average number of peptides per protein is " << (((double) accumulate(counts.begin(), counts.end(), 0.)) / counts.size()) << endl;
    os << "Average number of un/modified peptides per protein is " << (((double) accumulate(mod_counts.begin(), mod_counts.end(), 0.)) / mod_counts.size()) << endl;
    os << "Number of identified spectra: " << spectrum_count << endl;
    os << "Number of unique identified peptides: " << unique_peptides.size() << endl;

    vector<double>::iterator it = statistics.begin();
    vector<Size>::iterator it2 = counts.begin();
    vector<Size>::iterator it3 = mod_counts.begin();
    while (it != statistics.end())
    {
      if (*it == 0.)
      {
        it = statistics.erase(it);
        it2 = counts.erase(it2);
        it3 = mod_counts.erase(it3);
      }
      else
      {
        ++it;
        ++it2;
        ++it3;
      }
    }
    os << "Average coverage per found protein (" << statistics.size() << ") is " << (accumulate(statistics.begin(), statistics.end(), 0.) / statistics.size()) << endl;
    os << "Average number of peptides per found protein is " << (((double) accumulate(counts.begin(), counts.end(), 0.)) / counts.size()) << endl;
    os << "Average number of un/modified peptides per protein is " << (((double) accumulate(mod_counts.begin(), mod_counts.end(), 0.)) / mod_counts.size()) << endl;

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    String out = getStringOption_("out");
    return outputTo_(OpenMS_Log_info, out);
  }

};


int main(int argc, const char** argv)
{
  TOPPSequenceCoverageCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
