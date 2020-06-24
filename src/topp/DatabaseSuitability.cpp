// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cfloat>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_DatabaseSuitability DatabaseSuitability

@brief Calculates a suitability for a database which was used a for peptide identification search. Also reports the quality of LC-MS spectra.

This tool uses the metrics and algorithms first presented in 'Assessing protein sequence database suitability using de novo sequencing' by Richard S. Johnson, Brian C. Searle, Brook L. Nunn, Jason M. Gilmore, Molly Phillips, Chris T. Amemiya, Michelle Heck, Michael J. MacCoss.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

Citation c = { "Richard S. Johnson, Brian C. Searle, Brook L. Nunn, Jason M. Gilmore, Molly Phillips, Chris T. Amemiya, Michelle Heck, Michael J. MacCoss",
                    "Assessing protein sequence database suitability using de novo sequencing",
                    "Molecular & Cellular Proteomics. January 1, 2020; 19, 1: 198-208",
                    "10.1074/mcp.TIR119.001752" };

class DatabaseSuitability :
  public TOPPBase
{
public:
  DatabaseSuitability() :
    TOPPBase("DatabaseSuitability", "Computes a suitability score for a database which was used for a peptide identification search. Also reports the quality of LC-MS spectra.", false, {c})
  {
  }
protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_id", "<file>", "", "Input idXML file from peptide search with combined database with added de novo peptide (after FDR)");
    setValidFormats_("in_id", { "idXML" });
    registerInputFile_("in_spec", "<file>", "", "Input MzML file used for the peptide identification");
    setValidFormats_("in_spec", { "mzML" });
    registerInputFile_("in_novo", "<file>", "", "Input idXML file containing de novo peptides");
    setValidFormats_("in_novo", { "idXML" });
    registerOutputFile_("out", "<file>", "", "Optional tsv output containing database suitability information as well as spectral quality.", false);
    setValidFormats_("out", { "tsv" });
    registerDoubleOption_("novor_fract", "<double>", 1, "Set the fraction of how many cases, where a de novo peptide scores just higher than the database peptide, you wish to re-rank.", false, true);
    setMinFloat_("novor_fract", 0);
    setMaxFloat_("novor_fract", 1);
    registerFlag_("force_no_re_rank", "Use this flag if you want to disable re-ranking. Cases, where a de novo peptide scores just higher than the database peptide, are overlooked and counted as a de novo hit. This might underestimate the database quality.", true);
  }


  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in_id = getStringOption_("in_id");
    String in_spec = getStringOption_("in_spec");
    String in_novo = getStringOption_("in_novo");
    String out = getStringOption_("out");
    double novo_fract = getDoubleOption_("novor_fract");
    bool no_re_rank = getFlag_("force_no_re_rank");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    IdXMLFile x;
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    x.load(in_id, prot_ids, pep_ids);

    vector<ProteinIdentification> novo_prots;
    vector<PeptideIdentification> novo_peps;
    x.load(in_novo, novo_prots, novo_peps);

    // load mzML file in scope because we only need the number of ms2 spectra and no data
    // this saves some memory
    UInt64 count_ms2_lvl;
    {
      MzMLFile m;
      PeakFileOptions op;
      op.setMSLevels({2}); // only ms2
      op.setFillData(false); // no data
      m.setOptions(op);
      PeakMap exp;
      m.load(in_spec, exp);
      count_ms2_lvl = exp.size();
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // db suitability

    double cut_off;
    if (!no_re_rank)
    {
      cut_off = getDecoyCutOff_(pep_ids, novo_fract);
      if (cut_off == DBL_MAX)
      {
        OPENMS_LOG_ERROR << "Could not compute decoy cut off. Re-ranking impossible. If you want to ignore this, set the 'force_no_re_rank' flag." << endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    UInt64 count_db = 0;
    UInt64 count_novo = 0;
    UInt64 count_re_ranked = 0;
    UInt64 count_interest = 0;

    for (const auto& pep_id : pep_ids)
    {
      const vector<PeptideHit>& hits = pep_id.getHits();

      if (hits.empty()) continue;

      const PeptideHit& top_hit = hits[0];

      // skip if the top hit is a decoy hit
      if (top_hit.getMetaValue("target_decoy") == "decoy") continue;

      // check if top hit is found in de novo protein
      if (!isNovoHit_(top_hit)) // top hit is db hit
      {
        ++count_db;
        continue;
      }
      else // top hit is novo hit
      {
        if (hits.size() == 1)
        {
          ++count_novo;
          continue;
        }

        // find the second target hit, skip all decoy hits inbetween
        const PeptideHit* second_hit = nullptr;
        String target = "target";
        for (UInt i = 1; i < hits.size(); ++i)
        {
          if (target.find(hits[i].getMetaValue("target_decoy"), 0) == 0) // also check for "target+decoy" value
          {
            second_hit = &hits[i];
          }
        }
        if (second_hit == nullptr)
        {
          ++count_novo;
          continue;
        }

        // check if second hit is db hit
        if (isNovoHit_(*second_hit)) // second hit is also de novo hit
        {
          ++count_novo;
        }
        else // second hit is db hit
        {
          ++count_interest;
          // check for re-ranking
          if (no_re_rank)
          {
            ++count_novo;
            continue;
          }

          // check for xcorr score
          if (!top_hit.metaValueExists("MS:1002252") || !(*second_hit).metaValueExists("MS:1002252"))
          {
            throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported right now."));
          }

          if (double(top_hit.getMetaValue("MS:1002252")) - double((*second_hit).getMetaValue("MS:1002252")) <= cut_off)
          {
            ++count_db;
            ++count_re_ranked;
          }
          else ++count_novo;
        }
      }
    }

    double suitability = double(count_db) / (count_db + count_novo); //db suitability

    // spectra quality

    UInt64 count_novo_seq = 0;
    set<AASequence> unique_novo;

    for (const auto& pep_id : novo_peps)
    {
      if (pep_id.getHits().empty()) continue;
      ++count_novo_seq;
      unique_novo.insert(pep_id.getHits()[0].getSequence());
    }

    double id_rate = double(count_novo_seq) / count_ms2_lvl; // spectral quality (id rate of novo seqs)

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    OPENMS_LOG_INFO << count_db << " / " << (count_db + count_novo) << " top hits were found in the database." << endl;
    OPENMS_LOG_INFO << count_novo << " / " << (count_db + count_novo) << " top hits were only found in the concatenated de novo peptide." << endl;
    OPENMS_LOG_INFO << count_interest << " times scored a de novo hit above a database hit. Of those times " << count_re_ranked << " top de novo hits where re-ranked using a decoy cut-off of " << cut_off << endl;
    OPENMS_LOG_INFO << "database suitability [0, 1]: " << suitability << endl << endl;
    OPENMS_LOG_INFO << unique_novo.size() << " / " << count_novo_seq << " de novo sequences are unique" << endl;
    OPENMS_LOG_INFO << count_ms2_lvl << " ms2 spectra found" << endl;
    OPENMS_LOG_INFO << "spectral quality (id rate of de novo sequences) [0, 1]: " << id_rate << endl << endl;

    if (!out.empty())
    {
      OPENMS_LOG_INFO << "Writing output to: " << out << endl << endl;

      std::ofstream os(out);
      os.precision(writtenDigits(double()));
      os << "key" << "value\n";
      os << "#top_db_hits\t" << count_db << "\n";
      os << "#top_novo_hits\t" << count_novo << "\n";
      os << "db_suitability\t" << suitability << "\n";
      os << "#total_novo_seqs\t" << count_novo_seq << "\n";
      os << "#unique_novo_seqs\t" << unique_novo.size() << "\n";
      os << "#ms2_spectra" << count_ms2_lvl << "\n";
      os << "spectral_quality\t" << id_rate << "\n";
      os.close();
    }

    return EXECUTION_OK;
  }
private:
  // Calculates the difference of the xcorr scores from the first two decoy hits in a peptide identification.
  // If there aren't at least two decoy hits in the top ten, DBL_MAX is returned.
  // Also checks for target-decoy information and if FDR was run.
  double getDecoyDiff_(const PeptideIdentification& pep_id)
  {
    double diff = DBL_MAX;

    // get the score of the first two decoy hits
    double decoy_1 = DBL_MAX;
    double decoy_2 = DBL_MAX;
    UInt curr_hit = 0;

    for (const auto& hit : pep_id.getHits())
    {
      if (curr_hit > 10) break;
      ++curr_hit;

      if (!hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
      }

      if (pep_id.getScoreType() != "q-value")
      {
        if (!hit.metaValueExists("q-value"))
        {
          throw(Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No q-value found at peptide identification nor at peptide hits. Make sure 'False Discovery Rate' is run beforehand."));
        }
      }

      if (!hit.metaValueExists("MS:1002252"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported right now."));
      }

      if (decoy_1 == DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_1 = hit.getMetaValue("MS:1002252");
        continue;
      }
      if (decoy_1 < DBL_MAX && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_2 = hit.getMetaValue("MS:1002252");
        break;
      }
    }

    if (decoy_2 < DBL_MAX) // if there are two decoy hits
    {
      diff = abs(decoy_1 - decoy_2) / pep_id.getMZ(); // normalized by mw
    }

    // if there aren't two decoy hits DBL_MAX is returned
    return diff;
  }

  // Calculates all decoy differences of N given peptide identifications.
  // Returns the the (1-novor_fract)*N highest one.
  double getDecoyCutOff_(const vector<PeptideIdentification>& pep_ids, double novor_fract)
  {
    // get all decoy diffs of peptide ids with at least two decoy hits
    vector<double> diffs;
    for (const auto& pep_id : pep_ids)
    {
      double diff = getDecoyDiff_(pep_id);
      if (diff < DBL_MAX)
      {
        diffs.push_back(diff);
      }
    }

    if (double(diffs.size()) / pep_ids.size() < 0.2)
    {
      throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Under 20 % of peptide identifications have two decoy hits. This is not enough for re-ranking. Use the 'force_no_re_rank' flag to still compute a suitability score."));
    }

    // sort the diffs decreasing and get the (1-novo_fract)*N one
    auto sort_end = diffs.begin() + (1 - novor_fract) * diffs.size();
    
    std::partial_sort(diffs.begin(), sort_end + 1, diffs.end(), std::greater<double>());

    return *sort_end;
  }

  // Checks if all protein accessions contain 'CONCAT_PEPTIDE' (appended by IDFileConverter) and the hit is therefore considered a de novo hit.
  // If at least one accession doesn't contain 'CONCAT_PEPTIDE' the hit is considered a database hit.
  bool isNovoHit_(const PeptideHit& hit)
  {
    set<String> accessions = hit.extractProteinAccessionsSet();
    for (const String& acc : accessions)
    {
      if (acc.find(Constants::UserParam::CONCAT_PEPTIDE) == String::npos)
      {
        return false;
      }
    }
    return true;
  }
};


// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  DatabaseSuitability tool;
  return tool.main(argc, argv);
}

/// @endcond
