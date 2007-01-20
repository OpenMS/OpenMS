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
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page PILISModel PILISModel
	
	@brief Can be used to train the PILIS model with a given set of spectra an identifications

	The PILISModel can be trained with spectra and identifications to gain information about 
	the fragmentation of a specific instrument. Training spectra are read from the in parameter
	which should be in mzData format. The identifications of the spectra are read from file 
	given in the id_in parameter. The spectra and the identifications should be in the same order.
	
	The default parameters of the model are given in the model_file parameter, this can either be
	the default file stored in /data/PILIS/ or a custom modelfile with already trained 
	parameters. The trained parameters are written into the file trained_model_file. This parameter
	file can be used with PILISIdentification to generate identifications of MS/MS spectra.
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
			registerStringOption_("in", "<file>", "", "input file for the spectra in MzData format");
			registerStringOption_("id_in", "<file>", "", "input file for the annotations in AnalysisXML format");
			registerStringOption_("model_file", "<file>", "", "model file for training");
			registerStringOption_("trained_model_file", "<file>", "", "the output file of the trained model");
			addEmptyLine_();
		}
		
		ExitCodes main_(int , char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input/output files
			String in(getStringOption_("in"));
			String id_in(getStringOption_("id_in"));
			String model_file(getStringOption_("model_file"));
			String trained_model_file(getStringOption_("trained_model_file"));
						
      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      PeakMap exp;
      MzDataFile f;
      f.load(in, exp);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
			
			PILISModel* model = new PILISModel();
			model->readFromFile(model_file);

			AnalysisXMLFile id_in_file;
			vector<ProteinIdentification> prot_ids;
			vector<IdentificationData> peptide_ids;
			id_in_file.load(id_in, prot_ids, peptide_ids);

			if (peptide_ids.size() != exp.size())
			{
				throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "number of spectra and number of annotation from the file differ", "should be equal");
			}

			PeakMap::ConstIterator map_it = exp.begin();
			for (vector<IdentificationData>::const_iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it, ++map_it)
			{
				if (it->id.getPeptideHits().size() > 0)
				{
					PeptideHit hit = *it->id.getPeptideHits().begin();
					writeDebug_("Training with peptide: " + hit.getSequence() + " (z=" + String(hit.getCharge()) + ")", 1);
					model->train(*map_it, hit.getSequence(), hit.getCharge());
				}
			}
		
			model->evaluate();
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
		
			model->writeToFile(trained_model_file);
			
			// Write to model file
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, char ** argv )
{
	TOPPPILISModel tool;
	return tool.main(argc,argv);
}

