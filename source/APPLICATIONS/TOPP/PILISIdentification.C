// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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


#include <OpenMS/APPLICATIONS/TOPPBase2.h>
#include <OpenMS/ANALYSIS/ID/PILISIdentification.h>
#include <OpenMS/ANALYSIS/ID/PILISSequenceDB.h>
#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page PILISIdentification PILISIdentification
	
	@brief Performs an identification with PILIS
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPILISIdentification
	: public TOPPBase2
{
	public:
		TOPPPILISIdentification()
			: TOPPBase2("PILISIdentification", "can apply several spectra filters to the spectra")
		{
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerStringOption_("in", "<file>", "", "input file in MzData format");
			registerStringOption_("out", "<file>", "", "output file in MzData format");
			registerStringOption_("filters", "<filter1>[,<filter2>]", "", "filter to be applied");
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

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
			
			PILISModel* model = new PILISModel();
			model->initModel();

			PILISSequenceDB* db = new PILISSequenceDB();
			
			cerr << "beginning with ids" << endl;
			
			vector<Identification> ids;
			PILISIdentification PILIS_id;

			PILIS_id.setSequenceDB(db);
			PILIS_id.setModel(model);
			
			PILIS_id.getIdentifications(ids, exp);
		
			//debugging
			for (vector<Identification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
			{
				cerr << it->getPeptideHits().begin()->getSequence() << " " << it->getPeptideHits().begin()->getScore() << endl;
			}
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			f.store(out, exp);
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, char ** argv )
{
	TOPPPILISIdentification tool;
	return tool.main(argc,argv);
}

