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
#include <OpenMS/ANALYSIS/ID/PILISIdentification.h>
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISScoring.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/DATASTRUCTURES/SuffixArrayPeptideFinder.h>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_PILISIdentification PILISIdentification
	
	@brief Performs an ProteinIdentification with PILIS

	The PILISIdentification TOPP tool performs a ProteinIdentification run with 
	the PILIS ProteinIdentification engine. As input the file given in the in 
	parameters is used. The identifications are written into an IdXML
	file given in the out parameter. Additionally the model_file must be 
	specified. To perform a search also a peptide database file should be
	used,given in the peptide_db_file parameter. This should contain a 
	peptide in a separate line, either only the sequence or additionally 
	with weight and charge in the second and third column.
	
	@todo Check for missing precursors (Andreas)
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PILISIdentification.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPILISIdentification
	: public TOPPBase
{
	public:
		TOPPPILISIdentification()
			: TOPPBase("PILISIdentification", "performs a peptide/protein identification with the PILIS engine")
		{
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "input file in MzML format", true);
			registerOutputFile_("out", "<file>", "", "output file in IdXML format", true);
			registerInputFile_("model_file", "<file", "", "the model file of the PILISModel", true);
			registerInputFile_("peptide_db_file", "<file>", "", "a file which should contain peptides in the format\n"
																														 "DFPIANGER 1019.09 1\n"
																														 "where the first column is the peptide, the second the m/z\n"
																														 "the third the charge. As a alternative the sequence file\n"
																														 "may contain only peptide sequences each in a separate line\n"
																														 "repectively", true);
			registerDoubleOption_("precursor_mass_tolerance", "<tol>", 2.0 , "the precursor mass tolerance", false);
			registerDoubleOption_("peak_mass_tolerance", "<tol>", 1.0, "the peak mass tolerance", false);
			registerIntOption_("max_pre_candidates", "<int>", 200, "number of candidates that are used for precise scoring", false);
			registerIntOption_("max_candidates", "<int>", 20, "number of candidates that are reported by PILIS", false);
      registerDoubleOption_("upper_mz", "<double>", 2000.0, "bla", false);
			registerDoubleOption_("lower_mz", "<double>", 200.0, "bla", false);
			registerStringOption_("fixed_modifications", "<mods>", "", "monoisotopic_mass@residues e.g.: 57.021464@C", false);
	
			addEmptyLine_();
			addText_("Parameters of PILISModel");
			registerDoubleOption_("charge_directed_threshold", "<double>", 0.3, "bla", false);
			registerDoubleOption_("charge_remote_threshold", "<double>", 0.2, "bla", false);
			registerDoubleOption_("charge_loss_factor", "<double>", 0.5, "bla", false);
			//registerDoubleOption_("min_main_ion_intensity", "<double>", 0.02, "bla", false);
			//registerDoubleOption_("min_loss_ion_intensity", "<double>", 0.005, "bla", false);
			registerDoubleOption_("min_y_ion_intensity", "<double>", 0.20, "", false);
	    registerDoubleOption_("min_b_ion_intensity", "<double>", 0.15, "", false);
			registerDoubleOption_("min_a_ion_intensity", "<double>", 0.05, "", false);
			registerDoubleOption_("min_y_loss_intensity", "<double>", 0.05, "", false);
			registerDoubleOption_("min_b_loss_intensity", "<double>", 0.02, "", false);

			registerIntOption_("visible_model_depth", "<int>", 30, "bla", false);
			registerIntOption_("model_depth", "<int>", 4, "bla", false);

			addEmptyLine_();
			addText_("Parameters of PILISScoring");
			registerFlag_("use_local_scoring", "bla");
			registerFlag_("do_not_use_evalue_scoring", "bla");
			registerIntOption_("survival_function_bin_size", "<int>", 20, "bla", false);
			registerDoubleOption_("global_linear_fitting_threshold", "<double>", 0.1, "bla", false);
			registerDoubleOption_("local_linear_fitting_threshold", "<double>", 0.5, "bla", false);

			addEmptyLine_();
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input/output files
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      RichPeakMap exp;
      MzMLFile f;
      f.setLogType(log_type_);
      f.load(in, exp);

			writeDebug_("Data set contains " + String(exp.size()) + " spectra", 1);
			
      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
		
			writeDebug_("Reading model file", 2);

			// create model an set the given options
			PILISModel* model = new PILISModel();
			model->readFromFile(getStringOption_("model_file"));
			Param model_param(model->getParameters());
			model_param.setValue("upper_mz", getDoubleOption_("upper_mz"));
			model_param.setValue("lower_mz", getDoubleOption_("lower_mz"));
			model_param.setValue("charge_directed_threshold", getDoubleOption_("charge_directed_threshold"));
			model_param.setValue("charge_remote_threshold", getDoubleOption_("charge_remote_threshold"));
			//model_param.setValue("min_main_ion_intensity", getDoubleOption_("min_main_ion_intensity"));
			//model_param.setValue("min_loss_ion_intensity", getDoubleOption_("min_loss_ion_intensity"));
			model_param.setValue("min_y_ion_intensity", getDoubleOption_("min_y_ion_intensity"));
			model_param.setValue("min_b_ion_intensity", getDoubleOption_("min_b_ion_intensity"));
			model_param.setValue("min_a_ion_intensity", getDoubleOption_("min_a_ion_intensity"));
			model_param.setValue("min_y_loss_intensity", getDoubleOption_("min_y_loss_intensity"));
			model_param.setValue("min_b_loss_intensity", getDoubleOption_("min_b_loss_intensity"));
			model_param.setValue("charge_loss_factor", getDoubleOption_("charge_loss_factor"));
			model_param.setValue("visible_model_depth", getIntOption_("visible_model_depth"));
			model_param.setValue("model_depth", getIntOption_("model_depth"));
			model_param.setValue("fixed_modifications", getStringOption_("fixed_modifications"));
			model->setParameters(model_param);

			writeDebug_("Reading sequence db", 2);

			// create sequence db 
			SuffixArrayPeptideFinder*  sapf = new SuffixArrayPeptideFinder(getStringOption_("peptide_db_file"), "trypticCompressed");
			sapf->setTolerance(getDoubleOption_("precursor_mass_tolerance"));
			sapf->setNumberOfModifications(0);
			sapf->setUseTags(false);

			//exp.resize(50); // TODO
			
			UInt max_charge(3), min_charge(1); // TODO
			vector<double> pre_weights;
			for (RichPeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
			{
				double pre_weight(it->getPrecursors()[0].getMZ());
				for (Size z = min_charge; z <= max_charge; ++z)
				{
					pre_weights.push_back((pre_weight * (double)z) - (double)z);
				}
			}

			sort(pre_weights.begin(), pre_weights.end());

			cerr << "Getting candidates from SA...";
			vector<vector<pair<pair<String, String>, String> > > candidates;
			sapf->getCandidates(candidates, pre_weights);
			cerr << "done" << endl;

			delete sapf;

			map<double, vector<pair<pair<String, String>, String> > > sorted_candidates;
			UInt count(0);
			for (Size count = 0; count != candidates.size(); ++count)
			{
				sorted_candidates[pre_weights[count]] = candidates[count];
			}
			candidates.clear();

			// create ProteinIdentification and set the options
			PILISIdentification PILIS_id;
			
			PILIS_id.setModel(model);

			Param id_param(PILIS_id.getParameters());
			id_param.setValue("precursor_mass_tolerance", getDoubleOption_("precursor_mass_tolerance"));
			id_param.setValue("max_candidates", getIntOption_("max_pre_candidates"));
			// disable evalue scoring, this is done separately to allow for a single id per spectrum
			id_param.setValue("use_evalue_scoring", 0);
			id_param.setValue("fixed_modifications", getStringOption_("fixed_modifications"));
			PILIS_id.setParameters(id_param);

			vector<PeptideIdentification> ids;

			// perform the ProteinIdentification of the given spectra
			UInt no(0);
			for (RichPeakMap::Iterator it = exp.begin(); it != exp.end(); ++it, ++no)
			{
				if (it->getMSLevel() == 0)
				{
					writeLog_("Warning: MSLevel is 0, assuming MSLevel 2");
					it->setMSLevel(2);
				}

				if (it->getMSLevel() == 2)
				{
					writeDebug_(String(no) + "/" + String(exp.size()), 1);
					PeptideIdentification id;

					map<String, UInt> cand;
				
					for (UInt z = min_charge; z <= max_charge; ++z)
					{
						double pre_weight = (it->getPrecursors()[0].getMZ()* (double)z) - (double)z;
						for (vector<pair<pair<String, String>, String> >::const_iterator cit = sorted_candidates[pre_weight].begin(); cit != sorted_candidates[pre_weight].end(); ++cit)
						{
							String seq = cit->first.second;
							if (seq.size() > 39)
							{
								continue;
							}
							UInt num_cleavages_sites(0);
							for (Size k = 0; k != seq.size(); ++k)
							{
								if (k != seq.size() - 1)
								{
									if ((seq[k] == 'K' || seq[k] == 'R') && seq[k+1] != 'P')
									{
										++num_cleavages_sites;
									}
								}
							}

							if (num_cleavages_sites > 1)
							{
								continue;
							}
									
							cand[seq] = z;
						}
					}

					cerr << "#cand=" << cand.size() << endl;
					PILIS_id.getIdentification(cand, id, *it);
		
					id.setMetaValue("RT", it->getRT());
					id.setMetaValue("MZ", it->getPrecursors()[0].getMZ());

					ids.push_back(id);

					if (id.getHits().size() != 0)
					{
						cerr << it->getPrecursors()[0].getMZ() << " " << AASequence(id.getHits().begin()->getSequence()).getAverageWeight() << endl;
						writeDebug_(id.getHits().begin()->getSequence().toString() + " (z=" + id.getHits().begin()->getCharge() + "), score=" + String(id.getHits().begin()->getScore()) , 10);
					}
				}
			}

			// perform the PILIS scoring to the spectra
			if (!getFlag_("do_not_use_evalue_scoring"))
			{
				PILISScoring scoring;
				Param scoring_param(scoring.getParameters());
				scoring_param.setValue("use_local_scoring", (int)getFlag_("use_local_scoring"));
				scoring_param.setValue("survival_function_bin_size", getIntOption_("survival_function_bin_size"));
				scoring_param.setValue("global_linear_fitting_threshold", getDoubleOption_("global_linear_fitting_threshold"));
				scoring_param.setValue("local_linear_fitting_threshold", getDoubleOption_("local_linear_fitting_threshold"));
				scoring.setParameters(scoring_param);
	
				scoring.getScores(ids);
			}

			// write the result to the IdentificationData structure for the storing
			UInt max_candidates = getIntOption_("max_candidates");
			for (Size i = 0; i != ids.size(); ++i)
			{
				if (ids[i].getHits().size() > max_candidates)
				{
					vector<PeptideHit> hits = ids[i].getHits();
					hits.resize(max_candidates);
					ids[i].setHits(hits);
				}
			}
			
			delete model;


			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------

			DateTime now;
			now.now();
			
			String date_string;
			//now.get(date_string); // @todo Fix it (Andreas)
			String identifier("PILIS_"+date_string);

			//UInt count(0);
			count = 0;
			for (RichPeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it)
			{
				if (it->getMSLevel() == 2)
				{
					ids[count].setMetaValue("RT", it->getRT());
					ids[count].setMetaValue("MZ", it->getPrecursors()[0].getMZ());

					ids[count].setIdentifier(identifier);
					ids[count++].setHigherScoreBetter(false);
				}
			}

			// search parameters
			ProteinIdentification::SearchParameters search_parameters;
			search_parameters.db = getStringOption_("peptide_db_file");
			search_parameters.db_version = "";
			search_parameters.taxonomy = "";
			//search_parameters.charges = getStringOption_("charges");
			search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
			vector<String> fixed_mods;
			getStringOption_("fixed_modifications").split(',', fixed_mods);
			search_parameters.fixed_modifications = fixed_mods;
			search_parameters.enzyme = ProteinIdentification::TRYPSIN;
			search_parameters.missed_cleavages = 1;
			search_parameters.peak_mass_tolerance = getDoubleOption_("peak_mass_tolerance");
			search_parameters.precursor_tolerance = getDoubleOption_("precursor_mass_tolerance");

			ProteinIdentification protein_identification;
			protein_identification.setDateTime(now);
			protein_identification.setSearchEngine("PILIS");
			protein_identification.setSearchEngineVersion("beta");
			protein_identification.setSearchParameters(search_parameters);
			protein_identification.setIdentifier(identifier);

			vector<ProteinIdentification> protein_identifications;
			protein_identifications.push_back(protein_identification);
			IdXMLFile().store(out, protein_identifications, ids);
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPPILISIdentification tool;
	return tool.main(argc,argv);
}

