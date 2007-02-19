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
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page PILISIdentification PILISIdentification
	
	@brief Performs an identification with PILIS

	The PILISIdentification TOPP tool performs a identification run with 
	the PILIS identification engine. As input the file given in the in 
	parameters is used. The identifications are written into an AnalysisXML
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
			registerStringOption_("out", "<file>", "", "output file in AnalysisXML format", true);
			registerStringOption_("model_file", "<file", "", "the model file of the PILISModel", true);
			registerStringOption_("peptide_db_file", "<file>", "", "a file which should contain peptides in the format\n"
																														 "DFPIANGER 1019.09 1\n"
																														 "where the first column is the peptide, the second the m/z\n"
																														 "the third the charge. As a alternative the sequence file\n"
																														 "may contain only peptide sequences each in a separate line\n"
																														 "repectively", true);
			registerDoubleOption_("exponent", "<float>", 0.3, "exponent of the SpectrumAlignmentScore; see documentation of that class for more info", false);
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
      f.load(in, exp);

			writeDebug_("Data set contains " + String(exp.size()) + " spectra", 1);
			
      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
		
			writeDebug_("Reading model file", 2);

			PILISModel* model = new PILISModel();
			model->readFromFile(getStringOption_("model_file"));

			writeDebug_("Reading sequence db", 2);

			PILISSequenceDB* db = new PILISSequenceDB();
			db->addPeptidesFromFile(getStringOption_("peptide_db_file"));
			
			vector<Identification> ids;
			PILISIdentification PILIS_id;
			
			PILIS_id.setSequenceDB(db);
			PILIS_id.setModel(model);

			double exponent = getDoubleOption_("exponent");
			Param p(PILIS_id.getParameters());
			p.setValue("exponent", exponent);
			PILIS_id.setParameters(p);

			Size no(1);
			for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it, ++no)
			{
				if (it->getMSLevel() == 0)
				{
					writeLog_("Warning: MSLevel is 0, assuming MSLevel 2");
				}
							
				if (it->getMSLevel() == 2)
				{
					writeDebug_(String(no) + "/" + String(exp.size()), 1);
					//cerr << no << "/" << exp.size() << endl;
					Identification id;
					PILIS_id.getIdentification(id, *it);
					ids.push_back(id);
				}
			}
			
			delete model;
			delete db;


			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
		
			vector<IdentificationData> id_data;
			for (vector<Identification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
			{
				IdentificationData id_data_tmp;
				id_data_tmp.rt = 0;
				id_data_tmp.mz = 0;
				id_data_tmp.id = *it;
				id_data.push_back(id_data_tmp);
			}

			AnalysisXMLFile().store(out, vector<ProteinIdentification>(), id_data);
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, char ** argv )
{
	TOPPPILISIdentification tool;
	return tool.main(argc,argv);
}

