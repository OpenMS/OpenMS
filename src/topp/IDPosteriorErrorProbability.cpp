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
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

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
    bool ignore_bad_data = getFlag_("ignore_bad_data");
    bool prob_correct = getFlag_("prob_correct");
    String outlier_handling = fit_algorithm.getValue("outlier_handling");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    IdXMLFile file;
    vector<ProteinIdentification> protein_ids;
    vector<PeptideIdentification> peptide_ids;
    file.load(inputfile_name, protein_ids, peptide_ids);
    PosteriorErrorProbabilityModel PEP_model;
    PEP_model.setParameters(fit_algorithm);
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // check if there is a q-value score and target_decoy information available
    bool target_decoy_available(false);
    for (PeptideIdentification const & pep_id : peptide_ids)
    {
      const vector<PeptideHit>& hits = pep_id.getHits();
      if (!hits.empty())
      {
        target_decoy_available = (pep_id.getScoreType() == "q-value"
                               && hits[0].metaValueExists("target_decoy"));
        break;
      }
    }
    
    // map identifier "engine,charge" (if split_charge==true) or "engine" 
    // to three extracted score vectors. The main score vector contains the PSM scores.
    // Second and third are optional and contain target and decoy scores.
    map<String, vector<vector<double> > > all_scores = PosteriorErrorProbabilityModel::extractAndTransformScores(
     protein_ids, 
     peptide_ids, 
     split_charge,
     top_hits_only,
     target_decoy_available,
     fdr_for_targets_smaller);

    if (all_scores.empty())
    {
      writeLog_("No data collected. Check whether search engine is supported.");
      if (!ignore_bad_data) { return INPUT_FILE_EMPTY; }
    }

    String out_plot = fit_algorithm.getValue("out_plot").toString().trim();

    for (auto & score : all_scores)
    {
      vector<String> engine_info;
      score.first.split(',', engine_info);
      String engine = engine_info[0];
      Int charge = (engine_info.size() == 2) ? engine_info[1].toInt() : -1;

      if (split_charge)
      {
        // only adapt plot output if plot is requested (this badly violates the output rules and needs to change!)
        // one way to fix this: plot charges into a single file (no renaming of output file needed) - but this requires major code restructuring
        if (!out_plot.empty()) fit_algorithm.setValue("out_plot", out_plot + "_charge_" + String(charge));
        PEP_model.setParameters(fit_algorithm);
      }

      // fit to score vector
      //TODO choose outlier handling based on search engine? If not set by user?
      //XTandem is prone to accumulation at min values/censoring
      //OMSSA is prone to outliers
      bool return_value = PEP_model.fit(score.second[0], outlier_handling);

      if (!return_value) 
      {
        writeLog_("Unable to fit data. Algorithm did not run through for the following search engine: " + engine);
        if (!ignore_bad_data) { return UNEXPECTED_RESULT; }
      }

      if (return_value)
      {
        // plot target_decoy
        if (!out_plot.empty() 
         && top_hits_only 
         && target_decoy_available 
         && (!score.second[0].empty()))
        {
          PEP_model.plotTargetDecoyEstimation(score.second[1], score.second[2]); //target, decoy
        }
        
        bool unable_to_fit_data(true), data_might_not_be_well_fit(true);
        PosteriorErrorProbabilityModel::updateScores(
         PEP_model,
         engine,
         charge,
         prob_correct,
         split_charge,
         protein_ids,
         peptide_ids,
         unable_to_fit_data,
         data_might_not_be_well_fit);

        if (unable_to_fit_data)
        {
          writeLog_(String("Unable to fit data for search engine: ") + engine);
          if (!ignore_bad_data) return UNEXPECTED_RESULT;
        }
        else if (data_might_not_be_well_fit) 
        {
          writeLog_(String("Data might not be well fitted for search engine: ") + engine);
        }
      }
    }
    // Unfortunately this cannot go into the algorithm since
    // you would overwrite some score types before they are extracted when you
    // do split_charge
    for (auto& pep : peptide_ids)
    {
      if (prob_correct)
      {
        pep.setScoreType("Posterior Probability");
        pep.setHigherScoreBetter(true);
      }
      else
      {
        pep.setScoreType("Posterior Error Probability");
        pep.setHigherScoreBetter(false);
      }
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    file.store(outputfile_name, protein_ids, peptide_ids);
    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPPIDPosteriorErrorProbability tool;
  return tool.main(argc, argv);
}

/// @endcond
