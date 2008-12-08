// -*- mode: C++; tab-width: 2; -*-
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
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page DecoyDatabase DecoyDatabase
	
	@brief This small utility can create decoy databases used for decoy database searches.
		
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_DecoyDatabase.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDecoyDatabase
	: public TOPPBase
{
	public:
		TOPPDecoyDatabase()
			: TOPPBase("DecoyDatabase","Create decoy databases from normal ones", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input fasta file containing the database.");
			registerOutputFile_("out","<file>","","Output fasta file were the decoy database will be written to.");
			registerStringOption_("decoy_string", "<string>", "_ref", "String that is appended to the accession of the protein database to indicate a decoy protein.", false);
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));			
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			vector<FASTAFile::FASTAEntry> proteins;
			FASTAFile().load(in, proteins);
			

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------					

			String decoy_string(getStringOption_("decoy_string"));
			for (vector<FASTAFile::FASTAEntry>::iterator it = proteins.begin(); it != proteins.end(); ++it)
			{
				it->sequence.reverse();
				it->identifier += decoy_string;
			}
			
			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
		
			FASTAFile().store(out, proteins);
	
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPDecoyDatabase tool;
	return tool.main(argc,argv);
}
  
/// @endcond





