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

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_DatabaseSuitability DatabaseSuitability

@brief Calculates a suitability for a database which was used a for peptide identification search. Also reports the quality of LC-MS spectra.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class DatabaseSuitability :
  public TOPPBase
{
public:
  DatabaseSuitability() :
    TOPPBase("DatabaseSuitability", "Computes a suitability score for a database which was used for a peptide identification search. Also reports the quality of LC-MS spectra.", false)
  {
  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_id", "<file>", "", "Input idXML file from peptide search (after FDR)");
    setValidFormats_("in_id", { "idXML" });
    registerInputFile_("in_spec", "<file>", "", "Input MzML file");
    setValidFormats_("in_spec", { "mzML" });
    registerInputFile_("in_novo", "<file>", "", "Input idXML file containing de novo peptides");
    setValidFormats_("in_novo", { "idXML" });
    registerOutputFile_("out", "<file>", "", "Optional tsv output", false); //tsv output
    setValidFormats_("out", { "tsv" });
    registerIntOption_("novor_fract", "<integer>", 1, "Set the percentage of de novo peptides to capture with a score higher than the fasta score.", false, true);
    registerFlag_("force_no_re_rank", "Use this flag if you want to disable re-ranking. This might yeild in underperformance.", true);
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
    Int novo_fract = getIntOption_("novor_fract");
    bool no_re_rank = getFlag_("force_no_re_rank");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    IdXMLFile x;
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    x.load(in_id, prot_ids, pep_ids);

    // maybe extra tool for this

    vector<ProteinIdentification> novo_prots;
    vector<PeptideIdentification> novo_peps;
    x.load(in_novo, novo_prots, novo_peps);

    MzMLFile m;
    PeakFileOptions op;
    op.addMSLevel(2);
    m.setOptions(op);
    PeakMap exp;
    m.load(in_spec, exp);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // db suitability


    double cut_off;
    if (!no_re_rank)
    {
      cut_off = getDecoyCutOff_(pep_ids, novo_fract);
      if (cut_off < 0)
      {
        OPENMS_LOG_ERROR << "Could not compute decoy cut off. Re-ranking impossible. If you want to ignore this, set the 'force_no_re_rank' flag." << endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    UInt64 count_db = 0;
    UInt64 count_novo = 0;
    UInt64 count_re_ranked = 0;

    for (const auto& pep_id : pep_ids)
    {
      vector<PeptideHit> hits = pep_id.getHits();

      if (hits.empty()) continue;

      PeptideHit top_hit = hits[0];

      // check if top hit is found in de novo protein
      set<String> accessions = top_hit.extractProteinAccessionsSet();
      bool is_novo = true;
      for (const String& acc : accessions)
      {
        if (acc.find(Constants::UserParam::CONCAT_PEPTIDE) == String::npos)
        {
          is_novo = false;
          break;
        }
      }

      if (is_novo) // top hit is de novo hit
      {
        if (hits.size() == 1)
        {
          ++count_novo;
          continue;
        }

        PeptideHit second_hit = hits[1];

        // check if second hit is db hit
        set<String> second_accessions = top_hit.extractProteinAccessionsSet();
        bool is_novo_too = true;
        for (const String& acc2 : second_accessions)
        {
          if (acc2.find(Constants::UserParam::CONCAT_PEPTIDE) == String::npos)
          {
            is_novo_too = false;
            break;
          }
        }

        if (is_novo_too) // second hit is also de novo hit
        {
          ++count_novo;
        }
        else // second hit is db hit
        {
          // check for re-ranking
          if (no_re_rank)
          {
            ++count_novo;
            continue;
          }

          if (double(top_hit.getMetaValue("MS:1002252")) - double(second_hit.getMetaValue("MS:1002252")) <= cut_off)
          {
            ++count_db;
            ++count_re_ranked;
          }
          else ++count_novo;
        }
      }
      else ++count_db; // top hit is db hit
    }

    // spectra quality

    UInt64 count_ms2_lvl = exp.size();
    UInt64 count_novo_seq = 0;
    set<AASequence> unique_novo;

    for (const auto& pep_id : novo_peps)
    {
      if (pep_id.getHits().empty()) continue;
      ++count_novo_seq;
      unique_novo.insert(pep_id.getHits()[0].getSequence());
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    OPENMS_LOG_INFO << count_db << " top hits that were found in the database." << endl;
    OPENMS_LOG_INFO << count_novo << " top hits that were only found in the concatenated de novo peptide." << endl;
    OPENMS_LOG_INFO << count_re_ranked << " top de novo hits where re-ranked using a decoy cut-off of " << cut_off << endl;
    OPENMS_LOG_INFO << "Database quality: " << double(count_db) / (count_db + count_novo) << endl << endl;
    OPENMS_LOG_INFO << count_novo_seq << " de novo sequences derived from a total of " << count_ms2_lvl << " ms2 spectra. Ratio: " << double(count_novo_seq)/count_ms2_lvl << endl << endl;

    return EXECUTION_OK;

  }

private:

  double getDecoyDiff_(const PeptideIdentification& pep_id)
  {
    double diff = -1;

    // get the score of the first two decoy hits
    double decoy_1 = -1;
    double decoy_2 = -1;
    UInt curr_hit = 1;

    for (const auto& hit : pep_id.getHits())
    {
      if (curr_hit > 10) break;
      ++curr_hit;

      if (!hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run before hand."));
      }

      if (!hit.metaValueExists("MS:1002252"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported right now."));
      }

      if (decoy_1 == -1 && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_1 = hit.getMetaValue("MS:1002252");
        continue;
      }
      if (decoy_1 > 0 && hit.getMetaValue("target_decoy") == "decoy")
      {
        decoy_2 = hit.getMetaValue("MS:1002252");
        break;
      }
    }

    if (decoy_2 > 0) // if there are two decoy hits
    {
      diff = abs(decoy_1 - decoy_2) / pep_id.getMZ(); // normalized by mw
    }

    // if there aren't two decoy hits -1 is returned
    return diff;
  }

  double getDecoyCutOff_(const vector<PeptideIdentification>& pep_ids, double novor_fract)
  {
    double cut_off = -1;

    // get all decoy diffs of peptide ids with at least two decoy hits
    vector<double> diffs;
    for (const auto& pep_id : pep_ids)
    {
      double diff = getDecoyDiff_(pep_id);
      if (diff > 0)
      {
        diffs.push_back(diff);
      }
    }

    if (double(diffs.size()) / pep_ids.size() < 0.2)
    {
      throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Under 20 % of peptide identifications have two decoy hits. This is not enough for re-ranking. Use the 'force_no_re_rank' flag to still compute a suitability score."));
    }

    // sort the diffs decreasing
    std::sort(diffs.begin(), diffs.end(), std::greater<double>());
    
    // create a vector of percentages according to the number of differences
    vector<double> percent;
    for (Int i = 1; i <= diffs.size(); ++i)
    {
      percent.push_back(i/double(diffs.size()));
    }

    // find the right cut_off for the wanted percent of novo peptides to capture
    double fract = 1 - novor_fract;
    for (Int i = 0; i < percent.size(); ++i)
    {
      if (percent[i] > fract) cut_off = diffs[i];
    }

    return cut_off;
  }
};


// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  DatabaseSuitability tool;
  return tool.main(argc, argv);
}
