// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/PILISIdentification.h>
#include <OpenMS/ANALYSIS/ID/PILISSequenceDB.h>
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISScoring.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page PILISIdentification PILISIdentification
	
	@brief Performs an identification with PILIS

	The PILISIdentification TOPP tool performs a identification run with 
	the PILIS identification engine. As input the file given in the in 
	parameters is used. The identifications are written into an IdXML
	file given in the out parameter. Additionally the model_file must be 
	specified. To perform a search also a peptide database file should be
	used,given in the peptide_db_file parameter. This should contain a 
	peptide in a separate line, either only the sequence or additionally 
	with weight and charge in the second and third column.
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPILISIdentification
	: public TOPPBase
{
	public:
		TOPPPILISIdentification()
			: TOPPBase("PILISIdentification", "performs an identification with the PILIS engine")
		{
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerStringOption_("in", "<file>", "", "input file in MzData format", true);
			registerStringOption_("out", "<file>", "", "output file in IdXML format", true);
			registerStringOption_("model_file", "<file", "", "the model file of the PILISModel", true);
			registerStringOption_("peptide_db_file", "<file>", "", "a file which should contain peptides in the format\n"
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
			registerStringOption_("fixed_modifications", "<mods>", "", "<monoisotopic_mass>@<residues> e.g.: 57.021464@C", false);
	
			addEmptyLine_();
			addText_("Parameters of PILISModel");
			registerDoubleOption_("charge_directed_threshold", "<double>", 0.3, "bla", false);
			registerDoubleOption_("charge_remote_threshold", "<double>", 0.2, "bla", false);
			registerDoubleOption_("charge_loss_factor", "<double>", 0.5, "bla", false);
			registerDoubleOption_("min_main_ion_intensity", "<double>", 0.02, "bla", false);
			registerDoubleOption_("min_loss_ion_intensity", "<double>", 0.005, "bla", false);
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
		
		ExitCodes main_(int , char**)
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

      PeakMap exp;
      MzDataFile f;
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
			model_param.setValue("min_main_ion_intensity", getDoubleOption_("min_main_ion_intensity"));
			model_param.setValue("min_loss_ion_intensity", getDoubleOption_("min_loss_ion_intensity"));
			model_param.setValue("charge_loss_factor", getDoubleOption_("charge_loss_factor"));
			model_param.setValue("visible_model_depth", getIntOption_("visible_model_depth"));
			model_param.setValue("model_depth", getIntOption_("model_depth"));
			model_param.setValue("fixed_modifications", getStringOption_("fixed_modifications"));
			model->setParameters(model_param);

			writeDebug_("Reading sequence db", 2);

			// create sequence db 
			PILISSequenceDB* db = new PILISSequenceDB();
			db->addPeptidesFromFile(getStringOption_("peptide_db_file"));
		

			// create identification and set the options
			PILISIdentification PILIS_id;
			
			PILIS_id.setSequenceDB(db);
			PILIS_id.setModel(model);

			Param id_param(PILIS_id.getParameters());
			id_param.setValue("precursor_mass_tolerance", getDoubleOption_("precursor_mass_tolerance"));
			id_param.setValue("max_candidates", getIntOption_("max_pre_candidates"));
			// disable evalue scoring, this is done separately to allow for a single id per spectrum
			id_param.setValue("use_evalue_scoring", 0);
			id_param.setValue("fixed_modifications", getStringOption_("fixed_modifications"));
			PILIS_id.setParameters(id_param);

			vector<Identification> ids;
			vector<IdentificationData> id_data;

			// perform the identification of the given spectra
			UInt no(1);
			for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it, ++no)
			{
				if (it->getMSLevel() == 0)
				{
					writeLog_("Warning: MSLevel is 0, assuming MSLevel 2");
					it->setMSLevel(2);
				}
							
				if (it->getMSLevel() == 2)
				{
					writeDebug_(String(no) + "/" + String(exp.size()), 1);
					Identification id;
					PILIS_id.getIdentification(id, *it);
		
					ids.push_back(id);

					IdentificationData id_data_tmp;
					id_data_tmp.rt = it->getRT();
					id_data_tmp.mz = it->getPrecursorPeak().getPosition()[0];
					id_data_tmp.id = id;
					id_data.push_back(id_data_tmp);
				}
			}

			// perform the PILIS scoring to the spectra
			PILISScoring scoring;
			Param scoring_param(scoring.getParameters());
			scoring_param.setValue("use_local_scoring", (int)getFlag_("use_local_scoring"));
			scoring_param.setValue("survival_function_bin_size", getIntOption_("survival_function_bin_size"));
			scoring_param.setValue("global_linear_fitting_threshold", getDoubleOption_("global_linear_fitting_threshold"));
			scoring_param.setValue("local_linear_fitting_threshold", getDoubleOption_("local_linear_fitting_threshold"));
			scoring.setParameters(scoring_param);

			scoring.getScores(ids);

			// write the result to the IdentificationData structure for the storing
			UInt max_candidates = getIntOption_("max_candidates");
			for (UInt i = 0; i != ids.size(); ++i)
			{
				id_data[i].id = ids[i];
				if (id_data[i].id.getPeptideHits().size() > max_candidates)
				{
					id_data[i].id.getPeptideHits().resize(max_candidates);
				}
			}
			
			delete model;
			delete db;


			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
		
			IdXMLFile().store(out, vector<ProteinIdentification>(), id_data);
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, char ** argv )
{
	TOPPPILISIdentification tool;
	return tool.main(argc,argv);
}

