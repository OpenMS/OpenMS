// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Sven Nahnsen, Andreas Bertsch $
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
	@page UTILS_DecoyDatabase DecoyDatabase
	
	@brief Create decoy peptide databases from normal ones.
		
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
			: TOPPBase("DecoyDatabase","Create decoy peptide databases from normal ones.", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input fasta file containing the database.");
			registerOutputFile_("out","<file>","","Output fasta file where the decoy database will be written to.");
			registerStringOption_("decoy_string", "<string>", "_rev", "String that is appended to the accession of the protein database to indicate a decoy protein.", false);
			registerFlag_("append", "If this flag is used, the decoy database is appended to the target database, allowing combined target decoy searches.");
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			bool append = getFlag_("append");
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			vector<FASTAFile::FASTAEntry> proteins;
			FASTAFile().load(in, proteins);
			

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------					

			String decoy_string(getStringOption_("decoy_string"));
			Size num_proteins = proteins.size();
			set<String> identifiers;
			for (Size i = 0; i < num_proteins; ++i)
			{
				if (identifiers.find(proteins[i].identifier) != identifiers.end())
				{
					cerr << "DecoyDatabase: Warning, identifier is not unique to sequence file: '" << proteins[i].identifier << "'!" << endl;
				}
				identifiers.insert(proteins[i].identifier);

				if (append)
				{
					FASTAFile::FASTAEntry entry = proteins[i];
					entry.sequence.reverse();
					entry.identifier += decoy_string;
					proteins.push_back(entry);
				}
				else
				{
					proteins[i].sequence.reverse();
					proteins[i].identifier += decoy_string;
				}
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





