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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <boost/math/special_functions/fpclassify.hpp> // for "isnan"
#include <vector>

using namespace OpenMS;
using namespace Math; //PosteriorErrorProbabilityModel
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_IDPosteriorErrorProbability IDPosteriorErrorProbability

    @brief  Tool to estimate the probability of peptide hits to be incorrectly assigned.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ IDPosteriorErrorProbability \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ConsensusID </td>
        </tr>
    </table>
    </CENTER>

    @experimental This tool has not been tested thoroughly and might behave not as expected!

    By default an estimation is performed using the (inverse) Gumbel distribution for incorrectly assigned sequences
    and a Gaussian distribution for correctly assigned sequences. The probabilities are calculated by using Bayes' law, similar to PeptideProphet.
    Alternatively, a second Gaussian distribution can be used for incorrectly assigned sequences.
    At the moment, IDPosteriorErrorProbability is able to handle X! Tandem, Mascot, MyriMatch and OMSSA scores.

    No target/decoy information needs to be provided, since the model fits are done on the mixed distribution.

    In order to validate the computed probabilities an optional plot output can be generated.
    There are two parameters for the plot:
    The scores are plotted in the form of bins. Each bin represents a set of scores in a range of '(highest_score - smallest_score) / number_of_bins' (if all scores have positive values).
    The midpoint of the bin is the mean of the scores it represents.
    The parameter 'out_plot' should be used to give the plot a unique name. Two files are created. One with the binned scores and one with all steps of the estimation.
    If parameter @p top_hits_only is set, only the top hits of each peptide identification are used for the estimation process.
    Additionally, if 'top_hits_only' is set, target/decoy information is available and a @ref TOPP_FalseDiscoveryRate run was performed previously, an additional plot will be generated with target and decoy bins ('out_plot' must not be empty).
    A peptide hit is assumed to be a target if its q-value is smaller than @p fdr_for_targets_smaller.
    The plots are saved as a Gnuplot file. An attempt is made to call Gnuplot, which will create a PDF file containing all steps of the estimation. If this fails, the user has to run Gnuplot manually - or adjust the PATH environment such that Gnuplot can be found and retry.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_IDPosteriorErrorProbability.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_IDPosteriorErrorProbability.html

    For the parameters of the algorithm section see the algorithms documentation: @n
        @ref OpenMS::Math::PosteriorErrorProbabilityModel "fit_algorithm" @n

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDPosteriorErrorProbability :
  public TOPPBase
{
public:
  TOPPIDPosteriorErrorProbability() :
    TOPPBase("IDPosteriorErrorProbability", "Estimates probabilities for incorrectly assigned peptide sequences and a set of search engine scores using a mixture model.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerOutputFile_("out_plot", "<file>", "", "txt file (if gnuplot is available, a corresponding PDF will be created as well.)", false);
    setValidFormats_("out_plot", ListUtils::create<String>("txt"));

    registerFlag_("split_charge", "The search engine scores are split by charge if this flag is set. Thus, for each charge state a new model will be computed.");
    registerFlag_("top_hits_only", "If set only the top hits of every PeptideIdentification will be used");
    registerDoubleOption_("fdr_for_targets_smaller", "<value>", 0.05, "Only used, when top_hits_only set. Additionally, target/decoy information should be available. The score_type must be q-value from an previous False Discovery Rate run.", false, true);
    registerFlag_("ignore_bad_data", "If set errors will be written but ignored. Useful for pipelines with many datasets where only a few are bad, but the pipeline should run through.");
    registerFlag_("prob_correct", "If set scores will be calculated as '1 - ErrorProbabilities' and can be interpreted as probabilities for correct identifications.");
    registerSubsection_("fit_algorithm", "Algorithm parameter subsection");
    addEmptyLine_();
  }

  //there is only one parameter at the moment
  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param p = PosteriorErrorProbabilityModel().getParameters();
    if (p.exists("out_plot"))
    { // hide from user -- we have a top-level param for that
      p.remove("out_plot");
    }
    else 
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "INTERNAL ERROR: Param 'out_plot' was removed from fit-algorithm. Please update param handling internally!");
    }
    return p;
  }

  double getScore_(String& engine, const PeptideHit& hit)
  {
    if (engine == "OMSSA")
    {
      return (-1) * log10(max(hit.getScore(), smallest_e_value_));
    }
    else if (engine == "MyriMatch")
    {
      //double e_val = exp(-hit.getScore());
      //double score_val = ((-1)* log10(max(e_val,smallest_e_value_)));
      //printf("myri score: %e ; e_val: %e ; score_val: %e\n",hit.getScore(),e_val,score_val);
      //return score_val;
      return hit.getScore();
    }
    else if (engine.compare("XTandem") == 0)
    {
      return (-1) * log10(max((double)hit.getMetaValue("E-Value"), smallest_e_value_));
    }
    else if (engine == "MASCOT")
    {
      // issue #740: unable to fit data with score 0
      if (hit.getScore() == 0.0) 
      {
        return numeric_limits<double>::quiet_NaN();
      }
      // end issue #740
      if (hit.metaValueExists("EValue"))
      {
        return (-1) * log10(max((double)hit.getMetaValue("EValue"), smallest_e_value_));
      }
      if (hit.metaValueExists("expect"))
      {
        return (-1) * log10(max((double)hit.getMetaValue("expect"), smallest_e_value_));
      }
    }
    else if (engine == "SpectraST")
    {
      return 100 * hit.getScore(); // SpectraST f-val
    }
    else if (engine == "SimTandem")
    {
      if (hit.metaValueExists("E-Value"))
      {
        return (-1) * log10(max((double)hit.getMetaValue("E-Value"), smallest_e_value_));
      }
    }
    else if ((engine == "MSGFPlus") || (engine == "MS-GF+"))
    {
      if (hit.metaValueExists("MS:1002053"))  // name: MS-GF:EValue
      {
        return (-1) * log10(max((double)hit.getMetaValue("MS:1002053"), smallest_e_value_));
      }
      else if (hit.metaValueExists("expect"))
      {
        return (-1) * log10(max((double)hit.getMetaValue("expect"), smallest_e_value_));
      }
    }
    else if (engine == "Comet")
    {
      if (hit.metaValueExists("MS:1002257")) // name: Comet:expectation value
      {
        return (-1) * log10(max((double)hit.getMetaValue("MS:1002257"), smallest_e_value_));
      }
      else if (hit.metaValueExists("expect"))
      {
        return (-1) * log10(max((double)hit.getMetaValue("expect"), smallest_e_value_));
      }
    }
    else
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No parameters for chosen search engine", "The chosen search engine is currently not supported");
    }

    // avoid compiler warning (every code path must return a value, even if there is a throw() somewhere)
    return std::numeric_limits<double>::max();
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");
    Param fit_algorithm = getParam_().copy("fit_algorithm:", true);
    fit_algorithm.setValue("out_plot", getStringOption_("out_plot")); // re-assemble full param (was moved to top-level)
    bool split_charge = getFlag_("split_charge");
    bool top_hits_only = getFlag_("top_hits_only");
    double fdr_for_targets_smaller = getDoubleOption_("fdr_for_targets_smaller");
    bool target_decoy_available = false;
    bool ignore_bad_data = getFlag_("ignore_bad_data");
    bool prob_correct = getFlag_("prob_correct");

    // Set fixed e-value threshold
    smallest_e_value_ = numeric_limits<double>::denorm_min();

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    IdXMLFile file;
    vector<ProteinIdentification> protein_ids;
    vector<PeptideIdentification> peptide_ids;
    file.load(inputfile_name, protein_ids, peptide_ids);
    vector<double> scores;
    vector<double> decoy;
    vector<double> target;
    set<Int> charges;
    PosteriorErrorProbabilityModel PEP_model;
    PEP_model.setParameters(fit_algorithm);
    StringList search_engines = ListUtils::create<String>("XTandem,OMSSA,MASCOT,SpectraST,MyriMatch,SimTandem,MSGFPlus,MS-GF+,Comet");
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    if (split_charge)
    {
      for (vector<PeptideIdentification>::iterator pep_it = peptide_ids.begin(); pep_it != peptide_ids.end(); ++pep_it)
      {
        vector<PeptideHit>& hits = pep_it->getHits();
        for (std::vector<PeptideHit>::iterator hit_it = hits.begin(); hit_it != hits.end(); ++hit_it)
        {
          charges.insert(hit_it->getCharge());
        }
      }
      if (charges.empty())
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "no charges found!");
      }
    }
    for (vector<PeptideIdentification>::iterator pep_it = peptide_ids.begin(); pep_it != peptide_ids.end(); ++pep_it)
    {
      if (!pep_it->getHits().empty())
      {
        target_decoy_available = ((pep_it->getScoreType() == "q-value") && pep_it->getHits()[0].metaValueExists("target_decoy"));
        break;
      }
    }

    set<Int>::iterator charge_it = charges.begin(); // charges can be empty, no problem if split_charge is not set
    if (split_charge && charges.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'split_charge' is set, but the list of charge states is empty");
    }
    map<String, vector<vector<double> > > all_scores;
    char splitter = ','; // to split the engine from the charge state later on
    do
    {
      for (StringList::iterator engine_it = search_engines.begin(); engine_it != search_engines.end(); ++engine_it)
      {
        for (vector<ProteinIdentification>::iterator prot_it = protein_ids.begin(); prot_it != protein_ids.end(); ++prot_it)
        {
          String searchengine = prot_it->getSearchEngine();
          if ((*engine_it == searchengine) || (*engine_it == searchengine.toUpper()))
          {
            for (vector<PeptideIdentification>::iterator pep_it = peptide_ids.begin(); pep_it != peptide_ids.end(); ++pep_it)
            {
              if (prot_it->getIdentifier() == pep_it->getIdentifier())
              {
                vector<PeptideHit>& hits = pep_it->getHits();
                if (top_hits_only)
                {
                  pep_it->sort();
                  if (!hits.empty() && (!split_charge || hits[0].getCharge() == *charge_it))
                  {
                    double score = getScore_(*engine_it, hits[0]);
                    if (!boost::math::isnan(score)) // issue #740: ignore scores with 0 values, otherwise you will get the error "unable to fit data"
                    {
                      scores.push_back(score);

                      if (target_decoy_available)
                      {
                        if (hits[0].getScore() < fdr_for_targets_smaller)
                        {
                          target.push_back(score);
                        }
                        else
                        {
                          decoy.push_back(score);
                        }
                      }
                    }
                  }
                }
                else
                {
                  for (std::vector<PeptideHit>::iterator hit_it = hits.begin(); hit_it != hits.end(); ++hit_it)
                  {
                    if (!split_charge || (hit_it->getCharge() == *charge_it))
                    {
                      double score = getScore_(*engine_it, *hit_it);
                      if (!boost::math::isnan(score)) // issue #740: ignore scores with 0 values, otherwise you will get the error "unable to fit data"
                      {
                        scores.push_back(score);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        if (scores.size() > 2)
        {
          vector<vector<double> > tmp;
          tmp.push_back(scores);
          tmp.push_back(target);
          tmp.push_back(decoy);
          if (split_charge)
          {
            String engine_with_charge_state = *engine_it + String(splitter) + String(*charge_it);
            all_scores.insert(make_pair(engine_with_charge_state, tmp));
          }
          else
          {
            all_scores.insert(make_pair(*engine_it, tmp));
          }
        }

        scores.clear();
        target.clear();
        decoy.clear();
      }

      if (split_charge) ++charge_it;
    }
    while (charge_it != charges.end());

    if (all_scores.empty())
    {
      writeLog_("No data collected. Check whether search engine is supported.");
      if (!ignore_bad_data) return INPUT_FILE_EMPTY;
    }

    String out_plot = fit_algorithm.getValue("out_plot").toString().trim();
    for (map<String, vector<vector<double> > >::iterator score_it = all_scores.begin(); score_it != all_scores.end(); ++score_it)
    {
      vector<String> engine_info;
      score_it->first.split(splitter, engine_info);
      String engine = engine_info[0];
      Int charge = -1;
      if (engine_info.size() == 2)
      {
        charge = engine_info[1].toInt();
      }
      if (split_charge)
      {
        // only adapt plot output if plot is requested (this badly violates the output rules and needs to change!)
        // one way to fix this: plot charges into a single file (no renaming of output file needed) - but this requires major code restructuring
        if (!out_plot.empty()) fit_algorithm.setValue("out_plot", out_plot + "_charge_" + String(charge));
        PEP_model.setParameters(fit_algorithm);
      }

      const bool return_value = PEP_model.fit(score_it->second[0]);
      if (!return_value) writeLog_("Unable to fit data. Algorithm did not run through for the following search engine: " + engine);
      if (!return_value && !ignore_bad_data) return UNEXPECTED_RESULT;

      if (return_value)
      {
        // plot target_decoy
        if (!out_plot.empty() && top_hits_only && target_decoy_available && (score_it->second[0].size() > 0))
        {
          PEP_model.plotTargetDecoyEstimation(score_it->second[1], score_it->second[2]); //target, decoy
        }

        bool unable_to_fit_data = true;
        bool data_might_not_be_well_fit = true;
        for (vector<ProteinIdentification>::iterator prot_it = protein_ids.begin(); prot_it != protein_ids.end(); ++prot_it)
        {
          String searchengine = prot_it->getSearchEngine();
          if ((engine == searchengine) || (engine == searchengine.toUpper()))
          {
            for (vector<PeptideIdentification>::iterator pep_it = peptide_ids.begin(); pep_it != peptide_ids.end(); ++pep_it)
            {
              if (prot_it->getIdentifier() == pep_it->getIdentifier())
              {
                String score_type = pep_it->getScoreType() + "_score";
                vector<PeptideHit> hits = pep_it->getHits();
                for (std::vector<PeptideHit>::iterator hit_it = hits.begin(); hit_it != hits.end(); ++hit_it)
                {
                  if (!split_charge || (hit_it->getCharge() == charge))
                  {
                    double score;
                    hit_it->setMetaValue(score_type, hit_it->getScore());

                    score = getScore_(engine, *hit_it);
                    if (boost::math::isnan(score)) // issue #740: ignore scores with 0 values, otherwise you will get the error "unable to fit data"
                    {
                      score = 1.0;
                    }
                    else 
                    { 
                      score = PEP_model.computeProbability(score);
                      if ((score > 0.0) && (score < 1.0)) unable_to_fit_data = false;  // only if all it->second[0] are 0 or 1 unable_to_fit_data stays true
                      if ((score > 0.2) && (score < 0.8)) data_might_not_be_well_fit = false;  //same as above
                    }
                    hit_it->setScore(score);
                    if (prob_correct)
                    {
                      hit_it->setScore(1.0 - score);
                    }
                    else
                    {
                      hit_it->setScore(score);
                    }
                  }
                }
                pep_it->setHits(hits);
              }
              if (prob_correct)
              {
                pep_it->setScoreType("Posterior Probability");
                pep_it->setHigherScoreBetter(true);
              }
              else
              {
                pep_it->setScoreType("Posterior Error Probability");
                pep_it->setHigherScoreBetter(false);
              }
            }
          }
        }
        if (unable_to_fit_data) writeLog_(String("Unable to fit data for search engine: ") + engine);
        if (unable_to_fit_data && !ignore_bad_data) return UNEXPECTED_RESULT;

        if (data_might_not_be_well_fit) writeLog_(String("Data might not be well fitted for search engine: ") + engine);
      }
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    file.store(outputfile_name, protein_ids, peptide_ids);
    return EXECUTION_OK;
  }

  //Used in several functions
  double smallest_e_value_;
};


int main(int argc, const char** argv)
{
  TOPPIDPosteriorErrorProbability tool;

  return tool.main(argc, argv);
}

/// @endcond
