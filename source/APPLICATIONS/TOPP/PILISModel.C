// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
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
	@page TOPP_PILISModel PILISModel
	
	@brief Can be used to train the PILIS model with a given set of spectra an identifications

	This tool can be used in three different variants, 'training', 'cross_validation' and 
	'generation'.

	'training' mode:
	In training mode, the parameters for the fragmentation model needs to be set. Via 
	the -write_ini command line switch an ini file can be created, edited to the needs
	and used afterwards. Additionally, the spectra should be given as MSP file, which 
	already contains identifications or as mzML files. When using mzML files, idXML files
	must be used to get the peptide sequence information for the spectra.
	The tool trains then a model using the spectra and the peptides and writes it to
	the file given in the parameter 'trained_model_file'. Additionally, a model can
	be given as starting point via the parameter 'model_file'. With the min_charge and 
	max_charge parameters, the peptides can be restricted to the specified charge range.

	'cross_validation' mode:
	In cross validation mode a cross validation is performed to find the best parameters.
	The ini file contains for each parameter that can be optimized a flag, whether it
	should be used, a min value, a max value and a step size. These parameters are used
	to perform a grid search on the parameter. The result is a model with best performing
	parameter set. More on the cross validation can be found at the docu of the 
	PILISCrossValidation class.

	'generation' mode:
	This mode is used to generate spectra. A list of peptide must be given as
	idXML files. The peptides are used to generate spectra. Additionally a model file
	must be given, which contains the fragmentation model and its parameters. If a 
	peptide has charge 0, spectra for all charges from 'min_charge' to 'max_charge' 
	are generated.

	@experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PILISModel.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

// get a list of peptides and returns only those which are unique
void getUniquePeptides(vector<PILISCrossValidation::Peptide>& peptides)
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



class TOPPPILISModel
	: public TOPPBase
{
	public:
		TOPPPILISModel()
			: TOPPBase("PILISModel", "Used to trained the PILIS model with a given set of spectra an identifications")
		{
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			// input
			registerInputFileList_("in", "<file>", StringList(), "Input files for the spectra in MzML or MSP format.", false);
			setValidFormats_("in", StringList::create("mzML,MSP"));
			registerInputFileList_("id_in", "<file>", StringList(), "Input files for the annotations in IdXML format (if not given in MSP format).", false);
			setValidFormats_("id_in", StringList::create("idXML"));
			registerInputFile_("model_file", "<file>", "", "Input model file, used for generation mode or as basis for training. If not given, a default parameters are used for training.", false);

			// output
			registerOutputFile_("trained_model_file", "<file>", "", "The output file of the trained model, used in training mode.", false);
			registerOutputFile_("spectra_library_file", "<MSP-file>", "", "If this tool is used in generation mode, the spectral library is written into this MSP-file.", false);
			setValidFormats_("spectra_library_file", StringList::create("MSP"));
		
			// options
			registerStringOption_("type", "<usage-type>", "", "This parameter determines whether the model is used in 'training', 'cross_validation' or 'generation' mode.\n'training' is simply to train the model with the given spectra, using the parameters set in the ini file\n'cross_validation' performs a cross_validation using the identifications and the spectra, to find optimal parameters for the model\n'generation' generates a spectral library using a given model", true);
			setValidStrings_("type", StringList::create("training,cross_validation,generation"));

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

    Param getSubsectionDefaults_(const String& section) const
    {
			String type = getStringOption_("type");
      if (section == "PILIS_parameters")
      {
        return PILISModel().getParameters();
      }

			if (section == "cross_validation_parameters" && type == "cross_validation")
			{
				return PILISCrossValidation().getParameters();
			}

			if (section == "grid_search_parameters" && type == "cross_validation")
			{
				Param p;

				p.setValue("number_of_repeats", 2, "The grid search is performed 'number_of_repeats' times, to optimize the values.");
				p.setMinInt("number_of_repeats", 1);

				// lower_mz
				p.setValue("grid_search_lower_mz", "true", "Enables the grid search for the 'lower_mz' parameter", StringList::create("advanced"));
				p.setValidStrings("grid_search_lower_mz", StringList::create("true,false"));
				p.setValue("lower_mz_min", 0.0, "Minimal value of the 'lower_mz' parameter.", StringList::create("advanced"));
				p.setValue("lower_mz_max", 500.0, "Maximal value of the 'lower_mz' parameter.", StringList::create("advanced"));
				p.setValue("lower_mz_step_size", 20.0, "Step size for increasing the parameter 'lower_mz' during grid search", StringList::create("advanced"));

				// charge_remote_threshold
				p.setValue("grid_search_charge_remote_threshold", "true", "Enables the grid search for the parameter 'charge_remote_threshold'.", StringList::create("advanced"));
				p.setValidStrings("grid_search_charge_remote_threshold", StringList::create("true,false"));
      	p.setValue("charge_remote_threshold_min", 0.01, "Minimal value of the 'charge_remote_threshold' parameter.", StringList::create("advanced"));
      	p.setValue("charge_remote_threshold_max", 0.8, "Maximal value of the 'charge_remote_threshold' parameter.", StringList::create("advanced"));
      	p.setValue("charge_remote_threshold_step_size", 0.1, "Step size for increasing the parameter 'charge_remote_threshold' during the grid search.", StringList::create("advanced"));

				// charge_directed_threshold	
	      p.setValue("grid_search_charge_directed_threshold", "true", "Enables the grid search for the parameter 'charge_directed_threshold'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_charge_directed_threshold", StringList::create("true,false"));
        p.setValue("charge_directed_threshold_min", 0.0, "Minimal value of the 'charge_directed_threshold' parameter.", StringList::create("advanced"));
        p.setValue("charge_directed_threshold_max", 0.8, "Maximal value of the 'charge_directed_threshold' parameter.", StringList::create("advanced"));
        p.setValue("charge_directed_threshold_step_size", 0.1, "Step size for increasing the parameter 'charge_directed_threshold' during the grid search.", StringList::create("advanced"));

				// min_enhancement_factor
				p.setValue("grid_search_min_enhancement_factor", "true", "Enables the grid search for the parameter 'min_enhancement_factor'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_min_enhancement_factor", StringList::create("true,false"));
        p.setValue("min_enhancement_factor_min", 0.1, "Minimal value of the 'min_enhancement_factor' parameter.", StringList::create("advanced"));
        p.setValue("min_enhancement_factor_max", 2.0, "Maximal value of the 'min_enhancement_factor' parameter.", StringList::create("advanced"));
        p.setValue("min_enhancement_factor_step_size", 0.3, "Step size for increasing the parameter 'min_enhancement_factor' during the grid search.", StringList::create("advanced"));

				// side_chain_activation
        p.setValue("grid_search_side_chain_activation", "true", "Enables the grid search for the parameter 'side_chain_activation'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_side_chain_activation", StringList::create("true,false"));
        p.setValue("side_chain_activation_min", 0.0, "Minimal value of the 'side_chain_activation' parameter.", StringList::create("advanced"));
        p.setValue("side_chain_activation_max", 0.8, "Maximal value of the 'side_chain_activation' parameter.", StringList::create("advanced"));
        p.setValue("side_chain_activation_step_size", 0.05, "Step size for increasing the parameter 'side_chain_activation' during the grid search.", StringList::create("advanced"));

				// model_depth
        p.setValue("grid_search_model_depth", "true", "Enables the grid search for the parameter 'model_depth'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_model_depth", StringList::create("true,false"));
        p.setValue("model_depth_min", 4, "Minimal value of the 'model_depth' parameter.", StringList::create("advanced"));
        p.setValue("model_depth_max", 10, "Maximal value of the 'model_depth' parameter.", StringList::create("advanced"));
        p.setValue("model_depth_step_size", 1, "Step size for increasing the parameter 'model_depth' during the grid search.", StringList::create("advanced"));

				// min_a_ion_intensity
        p.setValue("grid_search_min_a_ion_intensity", "true", "Enables the grid search for the parameter 'min_a_ion_intensity'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_min_a_ion_intensity", StringList::create("true,false"));
        p.setValue("min_a_ion_intensity_min", 0.0, "Minimal value of the 'min_a_ion_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_a_ion_intensity_max", 0.5, "Maximal value of the 'min_a_ion_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_a_ion_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_a_ion_intensity' during the grid search.", StringList::create("advanced"));

				// min_b_ion_intensity
				p.setValue("grid_search_min_b_ion_intensity", "true", "Enables the grid search for the parameter 'min_b_ion_intensity'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_min_b_ion_intensity", StringList::create("true,false"));
        p.setValue("min_b_ion_intensity_min", 0.0, "Minimal value of the 'min_b_ion_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_b_ion_intensity_max", 0.8, "Maximal value of the 'min_b_ion_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_b_ion_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_b_ion_intensity' during the grid search.", StringList::create("advanced"));

				// min_y_ion_intensity
        p.setValue("grid_search_min_y_ion_intensity", "true", "Enables the grid search for the parameter 'min_y_ion_intensity'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_min_y_ion_intensity", StringList::create("true,false"));
        p.setValue("min_y_ion_intensity_min", 0.0, "Minimal value of the 'min_y_ion_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_y_ion_intensity_max", 0.8, "Maximal value of the 'min_y_ion_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_y_ion_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_y_ion_intensity' during the grid search.", StringList::create("advanced"));

				// min_b_loss_intensity
        p.setValue("grid_search_min_b_loss_intensity", "true", "Enables the grid search for the parameter 'min_b_loss_intensity'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_min_b_loss_intensity", StringList::create("true,false"));
        p.setValue("min_b_loss_intensity_min", 0.0, "Minimal value of the 'min_b_loss_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_b_loss_intensity_max", 0.5, "Maximal value of the 'min_b_loss_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_b_loss_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_b_loss_intensity' during the grid search.", StringList::create("advanced"));

				// min_y_loss_intensity
        p.setValue("grid_search_min_y_loss_intensity", "true", "Enables the grid search for the parameter 'min_y_loss_intensity'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_min_y_loss_intensity", StringList::create("true,false"));
        p.setValue("min_y_loss_intensity_min", 0.0, "Minimal value of the 'min_y_loss_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_y_loss_intensity_max", 0.5, "Maximal value of the 'min_y_loss_intensity' parameter.", StringList::create("advanced"));
        p.setValue("min_y_loss_intensity_step_size", 0.05, "Step size for increasing the parameter 'min_y_loss_intensity' during the grid search.", StringList::create("advanced"));

				// max_fragment_charge
        p.setValue("grid_search_max_fragment_charge", "true", "Enables the grid search for the parameter 'max_fragment_charge'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_max_fragment_charge", StringList::create("true,false"));
        p.setValue("max_fragment_charge_min", 1, "Minimal value of the 'max_fragment_charge' parameter.", StringList::create("advanced"));
        p.setValue("max_fragment_charge_max", 3, "Maximal value of the 'max_fragment_charge' parameter.", StringList::create("advanced"));
        p.setValue("max_fragment_charge_step_size", 1, "Step size for increasing the parameter 'max_fragment_charge' during the grid search.", StringList::create("advanced"));

				// max_isotope
        p.setValue("grid_search_max_isotope", "true", "Enables the grid search for the parameter 'max_isotope'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_max_isotope", StringList::create("true,false"));
        p.setValue("max_isotope_min", 1, "Minimal value of the 'max_isotope' parameter.", StringList::create("advanced"));
        p.setValue("max_isotope_max", 4, "Maximal value of the 'max_isotope' parameter.", StringList::create("advanced"));
        p.setValue("max_isotope_step_size", 1, "Step size for increasing the parameter 'max_isotope' during the grid search.", StringList::create("advanced"));

				// max_fragment_charge_training
        p.setValue("grid_search_max_fragment_charge_training", "true", "Enables the grid search for the parameter 'max_fragment_charge_training'.", StringList::create("advanced"));
        p.setValidStrings("grid_search_max_fragment_charge_training", StringList::create("true,false"));
        p.setValue("max_fragment_charge_training_min", 1, "Minimal value of the 'max_fragment_charge_training' parameter.", StringList::create("advanced"));
        p.setValue("max_fragment_charge_training_max", 3, "Maximal value of the 'max_fragment_charge_training' parameter.", StringList::create("advanced"));
        p.setValue("max_fragment_charge_training_step_size", 1, "Step size for increasing the parameter 'max_fragment_charge_training' during the grid search.", StringList::create("advanced"));

				return p;
			}

      return Param();
    }


		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------

			//input/output files
			StringList in(getStringList_("in"));
			StringList id_in(getStringList_("id_in"));
			String trained_model_file(getStringOption_("trained_model_file"));
			String model_file(getStringOption_("model_file"));
			String spectra_library_file(getStringOption_("spectra_library_file"));
			bool score_filtering(getFlag_("score_filtering"));
			DoubleReal score_threshold(getDoubleOption_("score_threshold"));
			Int min_charge(getIntOption_("min_charge"));
			Int max_charge(getIntOption_("max_charge"));

			String type = getStringOption_("type");


			if (type == "training")
			{
				if (in.size() == 0)
				{
					writeLog_("For 'training' mode spectra and identifications are needed.");
					return INCOMPATIBLE_INPUT_DATA;
				}
			} 
			else if (type == "cross_validation")
			{
				if (in.size() == 0)
				{
					writeLog_("For 'cross_validation' mode spectra and identification are needed.");
					return INCOMPATIBLE_INPUT_DATA;
				}
			} 
			else if (type == "generation")
			{
				// TODO
				if (spectra_library_file == "")
				{
					writeLog_("For 'generation' mode, the parameter 'spectra_library_file' must be given.");
					return MISSING_PARAMETERS;
				}

				if (model_file == "")
				{
					writeLog_("For 'generation' mode, the parameter 'model_file' must be given.");
					return MISSING_PARAMETERS;
				}
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
			ModificationDefinitionsSet mod_set((StringList)pilis_param.getValue("fixed_modifications"), (StringList)pilis_param.getValue("variable_modifications"));

			// read spectra file (if available)
			vector<RichPeakMap> exp;
			vector<vector<ProteinIdentification> > prot_ids;
			vector<vector<PeptideIdentification> > pep_ids;

			if (in.size() != 0)
			{
				FileTypes::Type in_file_type = FileHandler().getType(in[0]);
				writeDebug_("File type of parameter 'in' estimated as '" + FileHandler::typeToName(in_file_type) + "'", 1);
				// TODO check all types
				if (in_file_type == FileTypes::MSP)
				{
					writeDebug_("Reading MSP file" , 1);
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

			if (id_in.size() != 0)
			{
				prot_ids.resize(id_in.size());
				pep_ids.resize(id_in.size());
				IdXMLFile f;
				for (Size i = 0; i != id_in.size(); ++i)
				{
					f.load(id_in[i], prot_ids[i], pep_ids[i]);
				}
			}

			if (id_in.size() != 0 && in.size() != 0)
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
					if (it2->getPeptideIdentifications().size() == 0)
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

			if (type == "cross_validation")
			{
				PILISCrossValidation cv;
				Param cv_param = getParam_().copy("cross_validation_parameters:", true);
				cv.setParameters(cv_param);

				Param optimal_param = model.getParameters();

				Param grid_param = getParam_().copy("grid_search_parameters:", true);

				StringList double_parameters = StringList::create("lower_mz,charge_remote_threshold,charge_directed_threshold,min_enhancement_factor,min_y_ion_intensity,min_b_ion_intensity,min_a_ion_intensity,min_b_loss_intensity,min_y_loss_intensity,side_chain_activation");
				StringList int_parameters = StringList::create("max_isotope,max_fragment_charge,max_fragment_charge_training"); // todo add model_depth

				Size number_of_repeats = (UInt)grid_param.getValue("number_of_repeats");
				for (Size i = 0; i < number_of_repeats; ++i)
				{
					writeDebug_("Repeat " + String(i+1) + " of " + String(number_of_repeats), 1);
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
						DoubleReal min_value = (DoubleReal)grid_param.getValue(*it + "_min");
						DoubleReal max_value = (DoubleReal)grid_param.getValue(*it + "_max");
						DoubleReal step_size_value = (DoubleReal)grid_param.getValue(*it + "_step_size");
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
			}
			else if (type == "generation")
			{
				RichPeakMap exp;
				for (vector<vector<PeptideIdentification> >::const_iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
				{
					for (vector<PeptideIdentification>::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2)
					{
						if (it2->getHits().size() == 0)
						{
							continue;
						}
						PeptideHit hit = *it2->getHits().begin();
						Int charge = hit.getCharge();
						if (charge != 0)
						{
							RichPeakSpectrum spec;	
							model.getSpectrum(spec, hit.getSequence(), charge);
							spec.getPeptideIdentifications().push_back(*it2);
							exp.push_back(spec);
						}
						else
						{
							for (Int z = min_charge; z < max_charge; ++z)
							{
								RichPeakSpectrum spec;
								model.getSpectrum(spec, hit.getSequence(), z);
								
								PeptideIdentification id = *it2;
								vector<PeptideHit> hits = it2->getHits();
								hits.begin()->setCharge(z);
								id.setHits(hits);
								spec.getPeptideIdentifications().push_back(id);
								exp.push_back(spec);
							}
						}
					}
				}
			}
			else
			{
				model.setParameters(pilis_param);
				for (vector<PILISCrossValidation::Peptide>::const_iterator it = peptides.begin(); it != peptides.end(); ++it)
				{
					model.train(it->spec, it->sequence, it->charge);
				}
				model.evaluate();
			}

			if (trained_model_file != "")
			{
				model.writeToFile(trained_model_file);
			}


			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPPILISModel tool;
	return tool.main(argc,argv);
}

