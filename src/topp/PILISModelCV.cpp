// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Johannes Junker $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISCrossValidation.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
  @page TOPP_PILISModelCV PILISModelCV

  @brief Perform a cross validation of the PILIS model parameters
    @experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ PILISModelCV \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PILISIdentification </td>
        </tr>
  <tr>
    <td ROWSPAN=1></td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PILISModelTrainer </td>
  </tr>
    </table>
    </CENTER>

  A cross validation is performed to find the best parameters.
    The ini file contains for each parameter that can be optimized a flag, whether it
    should be used, a min value, a max value and a step size. These parameters are used
    to perform a grid search on the parameter. The result is a model with best performing
    parameter set. More on the cross validation can be found at the docu of the
    PILISCrossValidation class.

    <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_PILISModelCV.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_PILISModelCV.html
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

// get a list of peptides and returns only those which are unique
void getUniquePeptides(vector<PILISCrossValidation::Peptide> & peptides)
{
  vector<PILISCrossValidation::Peptide> unique_peptides;
  Map<AASequence, Map<Size, vector<PILISCrossValidation::Peptide> > > sorted;
  for (vector<PILISCrossValidation::Peptide>::const_iterator it = peptides.begin(); it != peptides.end(); ++it)
  {
    sorted[it->sequence][it->charge].push_back(*it);
  }

  // TODO set tic_filter option
  TICFilter tic_filter;
  for (Map<AASequence, Map<Size, vector<PILISCrossValidation::Peptide> > >::ConstIterator it1 = sorted.begin(); it1 != sorted.end(); ++it1)
  {
    for (Map<Size, vector<PILISCrossValidation::Peptide> >::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
    {
      double max_tic(0);
      PILISCrossValidation::Peptide pep;
      for (vector<PILISCrossValidation::Peptide>::const_iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3)
      {
        RichPeakSpectrum spec = it3->spec;

        double tic(tic_filter.apply(spec));
        if (tic > max_tic)
        {
          max_tic = tic;
          pep = *it3;
        }
      }
      unique_peptides.push_back(pep);
    }
  }

  peptides = unique_peptides;
}

class TOPPPILISModelCV :
  public TOPPBase
{
public:
  TOPPPILISModelCV() :
    TOPPBase("PILISModelCV", "Perform a cross validation of the PILIS model parameters")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    // input
    registerInputFileList_("in", "<file>", StringList(), "Input files for the spectra in mzML or MSP format.", false);
    setValidFormats_("in", ListUtils::create<String>("mzML,msp"));
    registerInputFileList_("id_in", "<file>", StringList(), "Input files for the annotations in idXML format (if not given in MSP format).", false);
    setValidFormats_("id_in", ListUtils::create<String>("idXML"));
    registerInputFile_("model_file", "<file>", "", "Input model file, used for generation mode or as basis for training. If not given, a default parameters are used for training.", false);

    // output
    registerOutputFile_("trained_model_file", "<file>", "", "The output file of the trained model, used in training mode.", false);

    registerIntOption_("min_charge", "<charge>", 1, "The minimal charge state used for training (other peptides are ignored) and for 'generation' mode if peptides have charge 0.", false);
    setMinInt_("min_charge", 1);
    registerIntOption_("max_charge", "<charge>", 3, "The maximal charge state used for training (other peptides are ignored) and for 'generation' mode if peptides have charge 0.", false);
    setMinInt_("max_charge", 1);
    registerFlag_("score_filtering", "If this flag is enabled the used spectra for training or cross validation are filtered using the 'score_treshold' parameter.");
    registerDoubleOption_("score_threshold", "<score>", 0, "The score threshold that must be passed in order to be used for training if 'score_filtering' is enabled.", false);

    addEmptyLine_();

    // subsections
    registerSubsection_("PILIS_parameters", "PILIS model parameters");
    registerSubsection_("cross_validation_parameters", "Parameters for the PILIS cross validation.");
    registerSubsection_("grid_search_parameters", "Parameters for the PILIS grid search.");
  }

  Param getSubsectionDefaults_(const String & section) const
  {
    if (section == "PILIS_parameters")
    {
      return PILISModel().getParameters();
    }

    if (section == "cross_validation_parameters")
    {
      return PILISCrossValidation().getParameters();
    }

    if (section == "grid_search_parameters")
    {
      Param p;

      p.setValue("number_of_repeats", 2, "The grid search is performed 'number_of_repeats' times, to optimize the values.");
      p.setMinInt("number_of_repeats", 1);

      // lower_mz
      p.setValue("grid_search_lower_mz", "true", "Enables the grid search for the 'lower_mz' parameter", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_lower_mz", ListUtils::create<String>("true,false"));
      p.setValue("lower_mz_min", 0.0, "Minimal value of the 'lower_mz' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("lower_mz_max", 500.0, "Maximal value of the 'lower_mz' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("lower_mz_step_size", 20.0, "Step size for increasing the parameter 'lower_mz' during grid search", ListUtils::create<String>("advanced"));

      // charge_remote_threshold
      p.setValue("grid_search_charge_remote_threshold", "true", "Enables the grid search for the parameter 'charge_remote_threshold'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_charge_remote_threshold", ListUtils::create<String>("true,false"));
      p.setValue("charge_remote_threshold_min", 0.01, "Minimal value of the 'charge_remote_threshold' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("charge_remote_threshold_max", 0.8, "Maximal value of the 'charge_remote_threshold' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("charge_remote_threshold_step_size", 0.1, "Step size for increasing the parameter 'charge_remote_threshold' during the grid search.", ListUtils::create<String>("advanced"));

      // charge_directed_threshold
      p.setValue("grid_search_charge_directed_threshold", "true", "Enables the grid search for the parameter 'charge_directed_threshold'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_charge_directed_threshold", ListUtils::create<String>("true,false"));
      p.setValue("charge_directed_threshold_min", 0.0, "Minimal value of the 'charge_directed_threshold' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("charge_directed_threshold_max", 0.8, "Maximal value of the 'charge_directed_threshold' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("charge_directed_threshold_step_size", 0.1, "Step size for increasing the parameter 'charge_directed_threshold' during the grid search.", ListUtils::create<String>("advanced"));

      // min_enhancement_factor
      p.setValue("grid_search_min_enhancement_factor", "true", "Enables the grid search for the parameter 'min_enhancement_factor'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_min_enhancement_factor", ListUtils::create<String>("true,false"));
      p.setValue("min_enhancement_factor_min", 0.1, "Minimal value of the 'min_enhancement_factor' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_enhancement_factor_max", 2.0, "Maximal value of the 'min_enhancement_factor' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_enhancement_factor_step_size", 0.3, "Step size for increasing the parameter 'min_enhancement_factor' during the grid search.", ListUtils::create<String>("advanced"));

      // side_chain_activation
      p.setValue("grid_search_side_chain_activation", "true", "Enables the grid search for the parameter 'side_chain_activation'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_side_chain_activation", ListUtils::create<String>("true,false"));
      p.setValue("side_chain_activation_min", 0.0, "Minimal value of the 'side_chain_activation' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("side_chain_activation_max", 0.8, "Maximal value of the 'side_chain_activation' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("side_chain_activation_step_size", 0.05, "Step size for increasing the parameter 'side_chain_activation' during the grid search.", ListUtils::create<String>("advanced"));

      // model_depth
      p.setValue("grid_search_model_depth", "true", "Enables the grid search for the parameter 'model_depth'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_model_depth", ListUtils::create<String>("true,false"));
      p.setValue("model_depth_min", 4, "Minimal value of the 'model_depth' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("model_depth_max", 10, "Maximal value of the 'model_depth' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("model_depth_step_size", 1, "Step size for increasing the parameter 'model_depth' during the grid search.", ListUtils::create<String>("advanced"));

      // min_a_ion_intensity
      p.setValue("grid_search_min_a_ion_intensity", "true", "Enables the grid search for the parameter 'min_a_ion_intensity'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_min_a_ion_intensity", ListUtils::create<String>("true,false"));
      p.setValue("min_a_ion_intensity_min", 0.0, "Minimal value of the 'min_a_ion_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_a_ion_intensity_max", 0.5, "Maximal value of the 'min_a_ion_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_a_ion_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_a_ion_intensity' during the grid search.", ListUtils::create<String>("advanced"));

      // min_b_ion_intensity
      p.setValue("grid_search_min_b_ion_intensity", "true", "Enables the grid search for the parameter 'min_b_ion_intensity'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_min_b_ion_intensity", ListUtils::create<String>("true,false"));
      p.setValue("min_b_ion_intensity_min", 0.0, "Minimal value of the 'min_b_ion_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_b_ion_intensity_max", 0.8, "Maximal value of the 'min_b_ion_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_b_ion_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_b_ion_intensity' during the grid search.", ListUtils::create<String>("advanced"));

      // min_y_ion_intensity
      p.setValue("grid_search_min_y_ion_intensity", "true", "Enables the grid search for the parameter 'min_y_ion_intensity'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_min_y_ion_intensity", ListUtils::create<String>("true,false"));
      p.setValue("min_y_ion_intensity_min", 0.0, "Minimal value of the 'min_y_ion_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_y_ion_intensity_max", 0.8, "Maximal value of the 'min_y_ion_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_y_ion_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_y_ion_intensity' during the grid search.", ListUtils::create<String>("advanced"));

      // min_b_loss_intensity
      p.setValue("grid_search_min_b_loss_intensity", "true", "Enables the grid search for the parameter 'min_b_loss_intensity'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_min_b_loss_intensity", ListUtils::create<String>("true,false"));
      p.setValue("min_b_loss_intensity_min", 0.0, "Minimal value of the 'min_b_loss_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_b_loss_intensity_max", 0.5, "Maximal value of the 'min_b_loss_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_b_loss_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_b_loss_intensity' during the grid search.", ListUtils::create<String>("advanced"));

      // min_y_loss_intensity
      p.setValue("grid_search_min_y_loss_intensity", "true", "Enables the grid search for the parameter 'min_y_loss_intensity'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_min_y_loss_intensity", ListUtils::create<String>("true,false"));
      p.setValue("min_y_loss_intensity_min", 0.0, "Minimal value of the 'min_y_loss_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_y_loss_intensity_max", 0.5, "Maximal value of the 'min_y_loss_intensity' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("min_y_loss_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_y_loss_intensity' during the grid search.", ListUtils::create<String>("advanced"));

      // max_fragment_charge
      p.setValue("grid_search_max_fragment_charge", "true", "Enables the grid search for the parameter 'max_fragment_charge'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_max_fragment_charge", ListUtils::create<String>("true,false"));
      p.setValue("max_fragment_charge_min", 1, "Minimal value of the 'max_fragment_charge' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("max_fragment_charge_max", 3, "Maximal value of the 'max_fragment_charge' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("max_fragment_charge_step_size", 1, "Step size for increasing the parameter 'max_fragment_charge' during the grid search.", ListUtils::create<String>("advanced"));

      // max_isotope
      p.setValue("grid_search_max_isotope", "true", "Enables the grid search for the parameter 'max_isotope'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_max_isotope", ListUtils::create<String>("true,false"));
      p.setValue("max_isotope_min", 1, "Minimal value of the 'max_isotope' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("max_isotope_max", 4, "Maximal value of the 'max_isotope' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("max_isotope_step_size", 1, "Step size for increasing the parameter 'max_isotope' during the grid search.", ListUtils::create<String>("advanced"));

      // max_fragment_charge_training
      p.setValue("grid_search_max_fragment_charge_training", "true", "Enables the grid search for the parameter 'max_fragment_charge_training'.", ListUtils::create<String>("advanced"));
      p.setValidStrings("grid_search_max_fragment_charge_training", ListUtils::create<String>("true,false"));
      p.setValue("max_fragment_charge_training_min", 1, "Minimal value of the 'max_fragment_charge_training' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("max_fragment_charge_training_max", 3, "Maximal value of the 'max_fragment_charge_training' parameter.", ListUtils::create<String>("advanced"));
      p.setValue("max_fragment_charge_training_step_size", 1, "Step size for increasing the parameter 'max_fragment_charge_training' during the grid search.", ListUtils::create<String>("advanced"));

      return p;
    }

    return Param();
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input/output files
    StringList in(getStringList_("in"));
    StringList id_in(getStringList_("id_in"));
    String trained_model_file(getStringOption_("trained_model_file"));
    String model_file(getStringOption_("model_file"));
    bool score_filtering(getFlag_("score_filtering"));
    double score_threshold(getDoubleOption_("score_threshold"));
    Int min_charge(getIntOption_("min_charge"));
    Int max_charge(getIntOption_("max_charge"));

    if (in.empty())
    {
      writeLog_("Spectra and identification are needed.");
      return INCOMPATIBLE_INPUT_DATA;
    }

    //bool duplicates_by_tic(getFlag_("duplicates_by_tic"));
    //bool base_model_from_file(getFlag_("base_model_from_file"));

    // create model, either read from a model file, or initialize with default parameters
    PILISModel model;
    if (model_file != "")
    {
      writeDebug_("Reading model from file '" + model_file + "'", 1);
      model.readFromFile(model_file);
    }
    else
    {
      writeDebug_("Initializing model", 1);
      model.setParameters(getParam_().copy("PILIS_parameters:", true));
      model.init();
    }

    Param pilis_param(model.getParameters());
    ModificationDefinitionsSet mod_set(pilis_param.getValue("fixed_modifications"), pilis_param.getValue("variable_modifications"));

    // read spectra file (if available)
    vector<RichPeakMap> exp;
    vector<vector<ProteinIdentification> > prot_ids;
    vector<vector<PeptideIdentification> > pep_ids;

    if (!in.empty())
    {
      FileTypes::Type in_file_type = FileHandler().getType(in[0]);
      writeDebug_("File type of parameter 'in' estimated as '" + FileTypes::typeToName(in_file_type) + "'", 1);
      // TODO check all types
      if (in_file_type == FileTypes::MSP)
      {
        writeDebug_("Reading MSP file", 1);
        MSPFile f;
        exp.resize(in.size());
        pep_ids.resize(in.size());
        for (Size i = 0; i != in.size(); ++i)
        {
          f.load(in[i], pep_ids[i], exp[i]);
          for (Size j = 0; j != exp[i].size(); ++j)
          {
            exp[i][j].getPeptideIdentifications().push_back(pep_ids[i][j]);
          }
        }
      }

      if (in_file_type == FileTypes::MZML)
      {
        MzMLFile f;
        f.setLogType(log_type_);

        exp.resize(in.size());
        for (Size i = 0; i != in.size(); ++i)
        {
          f.load(in[i], exp[i]);
        }
      }
    }

    if (!id_in.empty())
    {
      prot_ids.resize(id_in.size());
      pep_ids.resize(id_in.size());
      IdXMLFile f;
      for (Size i = 0; i != id_in.size(); ++i)
      {
        f.load(id_in[i], prot_ids[i], pep_ids[i]);
      }
    }

    if (!id_in.empty() && !in.empty())
    {
      // map the
      if (id_in.size() != in.size())
      {
        writeLog_("If in parameter contains mzML files and id_in contains idXML files, the number should be equal to allow mapping of the identification to the spectra");
        return INCOMPATIBLE_INPUT_DATA;
      }

      // map the ids to the spectra
      IDMapper id_mapper;
      for (Size i = 0; i != exp.size(); ++i)
      {
        id_mapper.annotate(exp[i], pep_ids[i], prot_ids[i]);
      }
    }

    // get the peptides and spectra
    vector<PILISCrossValidation::Peptide> peptides;

    for (vector<RichPeakMap>::const_iterator it1 = exp.begin(); it1 != exp.end(); ++it1)
    {
      for (RichPeakMap::ConstIterator it2 = it1->begin(); it2 != it1->end(); ++it2)
      {
        if (it2->getPeptideIdentifications().empty())
        {
          continue;
        }

        PeptideHit hit;

        if (it2->getPeptideIdentifications().begin()->getHits().size() > 0)
        {
          hit = *it2->getPeptideIdentifications().begin()->getHits().begin();
        }
        else
        {
          continue;
        }

        // check whether the sequence contains a modification not modelled
        if (!mod_set.isCompatible(hit.getSequence()) || hit.getSequence().size() > (UInt)pilis_param.getValue("visible_model_depth"))
        {
          continue;
        }

        if (score_filtering &&
            ((hit.getScore() < score_threshold && it2->getPeptideIdentifications().begin()->isHigherScoreBetter()) ||
             (hit.getScore() > score_threshold && !it2->getPeptideIdentifications().begin()->isHigherScoreBetter())))
        {
          continue;
        }

        PILISCrossValidation::Peptide pep_struct;
        pep_struct.sequence = hit.getSequence();
        pep_struct.charge = hit.getCharge();
        pep_struct.spec = *it2;
        pep_struct.hits = it2->getPeptideIdentifications().begin()->getHits();

        // check charges
        if (pep_struct.charge < min_charge || pep_struct.charge > max_charge)
        {
          continue;
        }

        peptides.push_back(pep_struct);
      }
    }

    getUniquePeptides(peptides);
    writeDebug_("Number of (unique) peptides for training: " + String(peptides.size()), 1);

    //model.writeToFile("pilis_tmp.dat");

    PILISCrossValidation cv;
    Param cv_param = getParam_().copy("cross_validation_parameters:", true);
    cv.setParameters(cv_param);

    Param optimal_param = model.getParameters();

    Param grid_param = getParam_().copy("grid_search_parameters:", true);

    StringList double_parameters = ListUtils::create<String>("lower_mz,charge_remote_threshold,charge_directed_threshold,min_enhancement_factor,min_y_ion_intensity,min_b_ion_intensity,min_a_ion_intensity,min_b_loss_intensity,min_y_loss_intensity,side_chain_activation");
    StringList int_parameters = ListUtils::create<String>("max_isotope,max_fragment_charge,max_fragment_charge_training");   // todo add model_depth

    Size number_of_repeats = (UInt)grid_param.getValue("number_of_repeats");
    for (Size i = 0; i < number_of_repeats; ++i)
    {
      writeDebug_("Repeat " + String(i + 1) + " of " + String(number_of_repeats), 1);
      for (StringList::const_iterator it = double_parameters.begin(); it != double_parameters.end(); ++it)
      {
        // check whether this parameters should be used for optimization
        bool enabled = DataValue(grid_param.getValue("grid_search_" + *it)).toBool();
        if (!enabled)
        {
          continue;
        }

        writeDebug_("Optimizing parameter '" + *it + "'", 1);

        model.setParameters(optimal_param);
        cv.setOptions(Map<String, PILISCrossValidation::Option>());
        double min_value = (double)grid_param.getValue(*it + "_min");
        double max_value = (double)grid_param.getValue(*it + "_max");
        double step_size_value = (double)grid_param.getValue(*it + "_step_size");
        cv.setOption(*it, PILISCrossValidation::Option(PILISCrossValidation::Option::DOUBLE, min_value, max_value, step_size_value));
        cv.apply(optimal_param, model, peptides);
      }

      for (StringList::const_iterator it = int_parameters.begin(); it != int_parameters.end(); ++it)
      {
        bool enabled = DataValue(grid_param.getValue("grid_search_" + *it)).toBool();
        if (!enabled)
        {
          continue;
        }

        writeDebug_("Optimizing parameter '" + *it + "'", 1);

        model.setParameters(optimal_param);
        cv.setOptions(Map<String, PILISCrossValidation::Option>());
        Int min_value = (Int)grid_param.getValue(*it + "_min");
        Int max_value = (Int)grid_param.getValue(*it + "_max");
        Int step_size_value = (Int)grid_param.getValue(*it + "_step_size");
        cv.setOption(*it, PILISCrossValidation::Option(PILISCrossValidation::Option::INT, min_value, max_value, step_size_value));
        cv.apply(optimal_param, model, peptides);
      }
    }

    // finally set the optimal parameters
    model.setParameters(optimal_param);

    if (trained_model_file != "")
    {
      model.writeToFile(trained_model_file);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPPILISModelCV tool;
  return tool.main(argc, argv);
}

/// @endcond

