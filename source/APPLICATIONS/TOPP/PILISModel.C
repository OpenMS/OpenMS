// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISModelGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_PILISModel PILISModel
	
	@brief Can be used to train the PILIS model with a given set of spectra an identifications

	The PILISModel can be trained with spectra and identifications to gain information about 
	the fragmentation of a specific instrument. Training spectra are read from the in parameter
	which should be in mzML format. The identifications of the spectra are read from file 
	given in the id_in parameter. The spectra and the identifications should be in the same order.
	
	The default parameters of the model are given in the model_file parameter, this can either be
	the default file stored in /data/PILIS/ or a custom modelfile with already trained 
	parameters. The trained parameters are written into the file trained_model_file. This parameter
	file can be used with @ref TOPP_PILISIdentification to generate identifications of MS/MS spectra.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PILISModel.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

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
			registerInputFile_("in", "<file>", "", "input file for the spectra in MzML format");
			registerInputFile_("id_in", "<file>", "", "input file for the annotations in IdXML format");
			registerOutputFile_("trained_model_file", "<file>", "", "the output file of the trained model");
		
			addEmptyLine_();		
			registerDoubleOption_("precursor_mass_tolerance", "<double>", 1.5, "precursor mass tolerance for the training", false);
			registerDoubleOption_("peak_mass_tolerance", "<double>", 0.3, "peak mass tolerance of the MS/MS spectra", false);
			

			registerFlag_("duplicates_by_tic", "duplicate sequence/charge combinations are filtered not by score but by TIC of the spectra");
			registerFlag_("use_tic_filtering", "Only use spectra filtered by the given TIC threshold");
			registerDoubleOption_("tic_threshold", "<double>", 10e10, "only spectra with TIC greater than the given threshold are used for training", false);
			registerFlag_("use_score_filtering", "Only use spectra filtered by the given score threshold");
			registerDoubleOption_("score_threshold", "<double>", 100, "only spectra with a score better than the given thresholde are used for training", false);
			registerStringOption_("fixed_modifications", "<mods>", "", "fixed modifications used for training", false);

			addEmptyLine_();
			registerFlag_("base_model_from_file", "if this flag is set, the model is not generated from scratch but read from the given 'model_file'");
			registerInputFile_("model_file", "<file>", "", "model file for training", false);
			registerIntOption_("model_depth", "<int>", 8, "model depth", false);
			registerIntOption_("visible_model_depth", "<int>", 40, "visible model depth", false);
			registerDoubleOption_("pseudo_counts", "<double>", 1e-15, "pseudo counts which are added when training the transition probabilties of the HMM", false);
			registerDoubleOption_("charge_remote_threshold", "<double>", 0.2, "charge remote threshold of the PILISModel", false);
			registerDoubleOption_("charge_directed_threshold", "<double>", 0.3, "charge directed threshold of the PILISModel", false);
			addEmptyLine_();
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input/output files
			String in(getStringOption_("in"));
			String id_in(getStringOption_("id_in"));
			//String model_file(getStringOption_("model_file"));
			String trained_model_file(getStringOption_("trained_model_file"));
			bool duplicates_by_tic(getFlag_("duplicates_by_tic"));
			bool base_model_from_file(getFlag_("base_model_from_file"));
			
			// fixed_modifications
			String fixed_modifications(getStringOption_("fixed_modifications"));
				

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      RichPeakMap exp;
      MzMLFile f;
      f.setLogType(log_type_);
      f.load(in, exp);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
			
			PILISModel* model = 0;
		
			writeDebug_("initializing the model", 5);
			
			if (base_model_from_file)
			{
				String model_file(getStringOption_("model_file"));
				model = new PILISModel();
				Param model_param(model->getParameters());
				model_param.setValue("model_depth", getIntOption_("model_depth"));
				model_param.setValue("visible_model_depth", getIntOption_("visible_model_depth"));
				model_param.setValue("fixed_modifications", fixed_modifications);
				model->setParameters(model_param);
				model->readFromFile(model_file);
			}
			else
			{
				PILISModelGenerator model_generator;
				Param p(model_generator.getParameters());
				p.setValue("model_depth", getIntOption_("model_depth"));
				p.setValue("visible_model_depth", getIntOption_("visible_model_depth"));
				p.setValue("fixed_modifications", fixed_modifications);
				model_generator.setParameters(p);
				model = new PILISModel(model_generator.getModel());
			}

			IdXMLFile id_in_file;
			vector<ProteinIdentification> prot_ids;
			vector<PeptideIdentification> peptide_ids;
			id_in_file.load(id_in, prot_ids, peptide_ids);

			Param model_param(model->getParameters());
			model_param.setValue("model_depth", getIntOption_("model_depth"));
			model_param.setValue("visible_model_depth", getIntOption_("visible_model_depth"));
			model_param.setValue("pseudo_counts", getDoubleOption_("pseudo_counts"));
			model_param.setValue("charge_remote_threshold", getDoubleOption_("charge_remote_threshold"));
			model_param.setValue("charge_directed_threshold", getDoubleOption_("charge_directed_threshold"));
			model_param.setValue("precursor_mass_tolerance", getDoubleOption_("precursor_mass_tolerance"));
			model_param.setValue("peak_mass_tolerance", getDoubleOption_("peak_mass_tolerance"));
			model_param.setValue("fixed_modifications", fixed_modifications);
			model->setParameters(model_param);


			Map<String, Map<UInt, Map<UInt, PeptideHit> > > ids; // [peptide][charge][index]
			for (Size i = 0; i != peptide_ids.size(); ++i)
			{
				if (peptide_ids[i].getHits().size() > 0)
				{
					PeptideHit hit = *peptide_ids[i].getHits().begin();
					if (hit.getCharge() <= 2) // TODO set option
					{
						String sequence = hit.getSequence().toString();
						replace(sequence.begin(), sequence.end(), 'L', 'I');
						ids[sequence][hit.getCharge()][i] = hit;
					}
				}
			}

			writeDebug_("There are " + String(exp.size()) + " spectra to train, with " + String(ids.size()) + " different sequences.", 1);

			if (peptide_ids.size() != exp.size())
			{
				throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "number of spectra and number of annotation from the file differ", "should be equal");
			}

			Map<UInt, PeptideHit> ids_to_train;

			if (!duplicates_by_tic)
			{
				// for each peptide sequence
				for (Map<String, Map<UInt, Map<UInt, PeptideHit> > >::ConstIterator it1 = ids.begin(); it1 != ids.end(); ++it1)
				{
					// for each charge
					for (Map<UInt, Map<UInt, PeptideHit> >::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
					{
						double score(-numeric_limits<double>::max());
						UInt max_idx(0);
						PeptideHit max_hit;
						bool has_max(false);
						// for each of the sequence charge spectra
						for (Map<UInt, PeptideHit>::ConstIterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3)
						{
							if (it3->second.getScore() >= score)
							{
								score = it3->second.getScore();
								max_idx = it3->first;
								max_hit = it3->second;
								has_max = true;
							}
						}
						if (has_max)
						{
							ids_to_train[max_idx] = max_hit;
						}
					}
				}
			}
			else
			{
				writeDebug_("Generating unique training spectra by TIC", 10);
				TICFilter tic_filter;
				// for each sequence
				for (Map<String, Map<UInt, Map<UInt, PeptideHit> > >::ConstIterator it1 = ids.begin(); it1 != ids.end(); ++it1)
				{
					// for each charge state
					for (Map<UInt, Map<UInt, PeptideHit> >::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
					{
						double tic(0);
						UInt max_idx(0);
						PeptideHit max_hit;
						// for each sequence charge combination spectra
						for (Map<UInt, PeptideHit>::ConstIterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3)
						{
							double actual_tic = tic_filter.apply(exp[it3->first]);
							if (actual_tic >= tic)
							{
								tic = actual_tic;
								max_idx = it3->first;
								max_hit = it3->second;
							}
						}
						ids_to_train[max_idx] = max_hit;
					}
				}
			}


			writeDebug_("Training with " + String(ids_to_train.size()) + " peptides", 1);

			double tic_threshold = getDoubleOption_("tic_threshold");
			double score_threshold = getDoubleOption_("score_threshold");
			bool use_tic_filtering = getFlag_("use_tic_filtering");
			bool use_score_filtering = getFlag_("use_score_filtering");
			UInt count(1);
			for (Map<UInt, PeptideHit>::ConstIterator it = ids_to_train.begin(); it != ids_to_train.end(); ++it, count++)
			{
				//if (it->second.getScore() >= threshold)
				//{
				TICFilter filter;
				double tic = filter.apply(exp[it->first]);
				if (!(use_tic_filtering && tic < tic_threshold) &&
						!(use_score_filtering && score_threshold < it->second.getScore()))
				{
					PeptideHit hit = it->second;
					String sequence = hit.getSequence().toString();
					replace(sequence.begin(), sequence.end(), 'L', 'I');
					try
					{
						AASequence seq(sequence);
						if (sequence.hasSubstring("("))
						{
							continue;
						}
					}
					catch(Exception::BaseException e)
					{
						writeDebug_(String("Error processing amino acid sequence: ")+e.what(), 1);
						continue;
					}
					writeDebug_(String(count) + String("/") + String(ids_to_train.size()) + String(", training model with peptide ") + sequence + String("z=") + String(hit.getCharge()), 3);
					model->train(exp[it->first], sequence, hit.getCharge());
				}
			}
	
			writeDebug_("Evaluating the model", 1);
			model->evaluate();
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
		
			writeDebug_("Writing the model file", 1);
			model->writeToFile(trained_model_file);
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPPILISModel tool;
	return tool.main(argc,argv);
}

