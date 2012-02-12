// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
  @page TOPP_PILISModelTrainer PILISModelTrainer

  @brief Train the PILIS model with a given set of spectra and identifications
	@experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

	<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ PILISModelTrainer \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PILISIdentification </td>
		</tr>
  <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PILISModelCV </td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PILISSpectraGenerator </td>
  </tr>
	</table>
	</CENTER>

  In order to train the model, the parameters of the fragmentation model need to be set first. Via
	the -write_ini command line switch an ini file can be created, edited to the needs
	and used afterwards. Additionally, the spectra should be given as MSP file, which
	already contains identifications or as mzML files. When using mzML files, idXML files
	must be used to get the peptide sequence information for the spectra.
	The tool trains then a model using the spectra and the peptides and writes it to
	the file given in the parameter 'trained_model_file'. Additionally, a model can
	be given as starting point via the parameter 'model_file'. With the min_charge and
	max_charge parameters, the peptides can be restricted to the specified charge range.

	<B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_PILISModelTrainer.cli
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



class TOPPPILISModelTrainer
	: public TOPPBase
{
	public:
    TOPPPILISModelTrainer()
      : TOPPBase("PILISModelTrainer", "Train the PILIS model with a given set of spectra and identifications")
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

			registerIntOption_("min_charge", "<charge>", 1, "The minimal charge state used for training (other peptides are ignored) and for 'generation' mode if peptides have charge 0.", false);
			setMinInt_("min_charge", 1);
			registerIntOption_("max_charge", "<charge>", 3, "The maximal charge state used for training (other peptides are ignored) and for 'generation' mode if peptides have charge 0.", false);
			setMinInt_("max_charge", 1);
			registerFlag_("score_filtering", "If this flag is enabled the used spectra for training or cross validation are filtered using the 'score_treshold' parameter.");
			registerDoubleOption_("score_threshold", "<score>", 0, "The score threshold that must be passed in order to be used for training if 'score_filtering' is enabled.", false);

			addEmptyLine_();

			// subsections
			registerSubsection_("PILIS_parameters", "PILIS model parameters");
    }

    Param getSubsectionDefaults_(const String& section) const
    {
      if (section == "PILIS_parameters")
      {
        return PILISModel().getParameters();
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
			bool score_filtering(getFlag_("score_filtering"));
			DoubleReal score_threshold(getDoubleOption_("score_threshold"));
			Int min_charge(getIntOption_("min_charge"));
			Int max_charge(getIntOption_("max_charge"));

      if (in.empty())
      {
        writeLog_("For 'training' mode spectra and identifications are needed.");
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
			ModificationDefinitionsSet mod_set((StringList)pilis_param.getValue("fixed_modifications"), (StringList)pilis_param.getValue("variable_modifications"));

			// read spectra file (if available)
			vector<RichPeakMap> exp;
			vector<vector<ProteinIdentification> > prot_ids;
			vector<vector<PeptideIdentification> > pep_ids;

      if ( !in.empty())
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

      if ( !id_in.empty() )
			{
				prot_ids.resize(id_in.size());
				pep_ids.resize(id_in.size());
				IdXMLFile f;
				for (Size i = 0; i != id_in.size(); ++i)
				{
					f.load(id_in[i], prot_ids[i], pep_ids[i]);
				}
			}

      if ( !id_in.empty() && !in.empty() )
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

      model.setParameters(pilis_param);
      for (vector<PILISCrossValidation::Peptide>::const_iterator it = peptides.begin(); it != peptides.end(); ++it)
      {
        model.train(it->spec, it->sequence, it->charge);
      }
      model.evaluate();

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
  TOPPPILISModelTrainer tool;
	return tool.main(argc,argv);
}

