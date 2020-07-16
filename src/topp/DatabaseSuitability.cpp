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

@brief Calculates the suitability of a database which was used a for peptide identification search. Also reports the quality of LC-MS spectra.

@dot
digraph sample_workflow {
  node [ style="solid,filled", color=black, fillcolor=grey90, width=1.5, fixedsize=true, shape=square, fontname=Helvetica, fontsize=10 ];
  edge [ arrowhead="open", style="solid" ];
  rankdir="LR";
  splines=ortho;
  mzml [ label="mzML file(s)" shape=oval fillcolor=white group=1];
  novor [ label="NovorAdapter" URL="\ref OpenMS::NovorAdapter" group=2];
  id_filter [ label="IDFilter" URL="\ref OpenMS::IDFilter" group=2];
  id_convert [ label="IDFileConverter" URL="\ref OpenMS::IDFileConverter" group=2];
  decoy_db [ label="DecoyDatabase" URL="\ref OpenMS::DecoyDatabase" group=2];
  comet [ label="CometAdapter" URL="\ref OpenMS::CometAdapter" group=1];
  pep_ind [ label="PeptideIndexer" URL="\ref OpenMS::PeptideIndexer" group=1];
  fdr [ label="FalseDiscoveryRate" URL="\ref OpenMS::FalseDiscoveryRate" group=1];
  db_suit [ label="DatabaseSuitability" fillcolor="#6F42C1" fontcolor=white group=3];
  tsv [ label="optional\ntsv output" shape=oval fillcolor=white group=3];
  {rank = same; db_suit; decoy_db;}
  mzml -> novor;
  mzml -> comet;
  comet -> pep_ind;
  pep_ind -> fdr;
  fdr -> db_suit [ xlabel="in_id" ];
  novor -> id_filter;
  id_filter -> id_convert;
  id_convert -> decoy_db;
  decoy_db -> comet;
  mzml -> db_suit [ xlabel="in_spec" ];
  novor -> db_suit [ xlabel="in_novor" ];
  db_suit -> tsv;
}
@enddot

The metric this tool uses to determine the suitability of a database is based on a de novo model. Therefore it is crucial that your workflow is set up the right way. Above you can see an example.@n
Most importantly the peptide identification search needs to be done with a combination of the database in question and a de novo "database".@n
To generate the de novo "database":
      - @ref UTILS_NovorAdapter calculates de novo sequences.
      - @ref TOPP_IDFilter can filter out unwanted ones.
      - @ref TOPP_IDFileConverter generates the de novo fasta file.

For re-ranking all cases where a peptide hit only found in the de novo "database" scores above a peptide hit found in the actual database are checked. In all these cases the cross-correlation scores of those peptide hits are compared. If they are similar enough, the database hit will be re-ranked to be on top of the de novo hit. You can control how much of cases with similar scores will be re-ranked by using the @p novor_fract option.@n
For this to work it is important that FDR filtering is done in this tool and not beforehand by @ref TOPP_FalseDiscoveryRate. You can provide the wanted FDR using the corresponding flag.

@note For identification search the only supported search engine for the time being is Comet because the Comet cross-correlation score is needed for re-ranking.@n
You can still uses other search engines and disable the re-ranking via the @p force_no_re_rank flag in this tool. This will probably result in an underestimated suitability though.@n


The results are written directly into the console. But you can provide an optional tsv output file where the most important results will be exported to.


This tool uses the metrics and algorithms first presented in:@n
<em>Assessing protein sequence database suitability using de novo sequencing. Molecular & Cellular Proteomics. January 1, 2020; 19, 1: 198-208. doi:10.1074/mcp.TIR119.001752.@n
Richard S. Johnson, Brian C. Searle, Brook L. Nunn, Jason M. Gilmore, Molly Phillips, Chris T. Amemiya, Michelle Heck, Michael J. MacCoss.</em>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_DatabaseSuitability.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_DatabaseSuitability.html
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
    registerDoubleOption_("FDR", "<double>", 0.01, "Filter peptide hits based on this q-value. (e.g., 0.05 = 5 % FDR)", false, true);
    setMinFloat_("FDR", 0);
    setMaxFloat_("FDR", 1);
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
    double FDR = getDoubleOption_("FDR");
    bool no_re_rank = getFlag_("force_no_re_rank");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    // load mzML file in scope because we only need the number of ms2 spectra and no data
    // this saves some memory
    Size count_ms2_lvl;
    {
      MzMLFile m;
      PeakFileOptions op;
      op.setMSLevels({ 2 }); // only ms2
      op.setFillData(false); // no data
      m.setOptions(op);
      PeakMap exp;
      m.load(in_spec, exp);
      count_ms2_lvl = exp.size();
    }
    
    IdXMLFile x;
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    x.load(in_id, prot_ids, pep_ids);

    vector<ProteinIdentification> novo_prots;
    vector<PeptideIdentification> novo_peps;
    x.load(in_novo, novo_prots, novo_peps);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // db suitability

    double cut_off{};
    if (!no_re_rank)
    {
      cut_off = getDecoyCutOff_(pep_ids, novo_fract);
      if (cut_off == DBL_MAX)
      {
        OPENMS_LOG_ERROR << "Could not compute decoy cut off. Re-ranking impossible. If you want to ignore this, set the 'force_no_re_rank' flag." << endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    Size count_db = 0;
    Size count_novo = 0;
    Size count_re_ranked = 0;
    Size count_interest = 0;

    for (auto& pep_id : pep_ids)
    {
      vector<PeptideHit>& hits = pep_id.getHits();
      bool q_value_score = (pep_id.getScoreType() == "q-value");

      if (hits.empty()) continue;

      // sort hits by q-value
      if (q_value_score)
      {
        sort(hits.begin(), hits.end(),
          [](const PeptideHit& a, const PeptideHit& b)
          {
            return a.getScore() < b.getScore();
          });
      }
      else
      {
        if (!hits[0].metaValueExists("q-value"))
        {
          throw(Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No q-value found at peptide identification nor at peptide hits. Make sure 'False Discovery Rate' is run beforehand."));
        }

        sort(hits.begin(), hits.end(),
          [](const PeptideHit& a, const PeptideHit& b)
          {
            return float(a.getMetaValue("q-value")) < float(b.getMetaValue("q-value"));
          });
      }


      const PeptideHit& top_hit = hits[0];

      // skip if the top hit is a decoy hit
      if (!top_hit.metaValueExists("target_decoy"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No target/decoy information found! Make sure 'PeptideIndexer' is run beforehand."));
      }
      if (top_hit.getMetaValue("target_decoy") == "decoy") continue;

      // skip if top hit is out ouf FDR
      if (scoreHigherThanFDR_(top_hit, FDR, q_value_score)) continue;

      // check if top hit is found in de novo protein
      if (!isNovoHit_(top_hit)) // top hit is db hit
      {
        ++count_db;
        continue;
      }

      // find the second target hit, skip all decoy or novo hits inbetween
      const PeptideHit* second_hit = nullptr;
      String target = "target";
      for (UInt i = 1; i < hits.size(); ++i)
        {
        // check for FDR
        if (scoreHigherThanFDR_(hits[i], FDR, q_value_score)) break;

        if (target.find(String(hits[i].getMetaValue("target_decoy"), 0)) == 0) // also check for "target+decoy" value
        {
          // check if hit is novo hit
          if (isNovoHit_(hits[i])) continue;
          
          second_hit = &hits[i];
          break;
        }
      }
      if (second_hit == nullptr) // no second target hit with given FDR found
      {
        ++count_novo;
        continue;
      }

      // second hit is db hit
      ++count_interest;

      // check for re-ranking
      if (no_re_rank)
      {
        ++count_novo;
        continue;
      }

      // check for xcorr score
      if (!top_hit.metaValueExists("MS:1002252") || !second_hit->metaValueExists("MS:1002252"))
      {
        throw(Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No cross correlation score found at peptide hit. Only Comet search engine is supported right now."));
      }
      
      double top_xscore_mw = double(top_hit.getMetaValue("MS:1002252")) / top_hit.getSequence().getMonoWeight();
      double second_xscore_mw = double(second_hit->getMetaValue("MS:1002252")) / second_hit->getSequence().getMonoWeight();
      if (top_xscore_mw - second_xscore_mw <= cut_off)
      {
        ++count_db;
        ++count_re_ranked;
      }
      else
      {
        ++count_novo;
      }
    }

    double suitability = double(count_db) / (count_db + count_novo); //db suitability

    // spectra quality

    Size count_novo_seq = 0;
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
    OPENMS_LOG_INFO << count_interest << " times scored a de novo hit above a database hit. Of those times " << count_re_ranked << " top de novo hits where re-ranked." << endl;
    OPENMS_LOG_INFO << "database suitability [0, 1]: " << suitability << endl << endl;
    OPENMS_LOG_INFO << unique_novo.size() << " / " << count_novo_seq << " de novo sequences are unique" << endl;
    OPENMS_LOG_INFO << count_ms2_lvl << " ms2 spectra found" << endl;
    OPENMS_LOG_INFO << "spectral quality (id rate of de novo sequences) [0, 1]: " << id_rate << endl << endl;

    if (!out.empty())
    {
      OPENMS_LOG_INFO << "Writing output to: " << out << endl << endl;

      std::ofstream os(out);
      os.precision(writtenDigits(double()));
      os << "key\tvalue\n";
      os << "#top_db_hits\t" << count_db << "\n";
      os << "#top_novo_hits\t" << count_novo << "\n";
      os << "db_suitability\t" << suitability << "\n";
      os << "#total_novo_seqs\t" << count_novo_seq << "\n";
      os << "#unique_novo_seqs\t" << unique_novo.size() << "\n";
      os << "#ms2_spectra\t" << count_ms2_lvl << "\n";
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
      diff = abs(decoy_1 - decoy_2) / pep_id.getHits()[0].getSequence().getMonoWeight(); // normalized by mw
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
    const set<String> accessions = hit.extractProteinAccessionsSet();
    for (const String& acc : accessions)
    {
      if (acc.find(Constants::UserParam::CONCAT_PEPTIDE) == String::npos)
      {
        return false;
      }
    }
    return true;
  }

  // Checks if the q-value of a peptide hit is higher than a given FDR.
  // Throws an error if no q-value is found.
  bool scoreHigherThanFDR_(const PeptideHit& hit, double FDR, bool q_value_score)
  {
    if (q_value_score) // score type is q-value
    {
      if (hit.getScore() > FDR) return true;
      return false;
    }
    
    if (hit.metaValueExists("q-value")) // look for q-value at metavalues
    {
      if (float(hit.getMetaValue("q-value")) > FDR) return true;
      return false;
    }
    
    // no q-value found
    throw(Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No q-value found at peptide identification nor at peptide hits. Make sure 'False Discovery Rate' is run beforehand."));
  }
};


// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  DatabaseSuitability tool;
  return tool.main(argc, argv);
}

/// @endcond
