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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page IdXMLInfo IdXMLInfo
	
	@brief This application is used to retrieve information about IdXML files.
	
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIdXMLInfo
	: public TOPPBase
{
	public:
		TOPPIdXMLInfo()
			: TOPPBase("IdXMLInfo","prints information about IdXML files")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file");
		}

		ExitCodes main_(int , const char**)
		{
			IdXMLFile idXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<PeptideIdentification> identifications;
			vector<PeptideHit> temp_peptide_hits;
			UInt spectrum_count = 0;
			UInt peptide_hit_count = 0;
			UInt runs_count = 0;
			UInt protein_hit_count = 0;
			String inputfile_name = "";
			
			protein_identifications.push_back(ProteinIdentification());
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			inputfile_name = getStringOption_("in");			
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			idXML_file.load(inputfile_name, protein_identifications, identifications);
			
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			for(UInt i = 0; i < identifications.size(); ++i)
			{
				if (!identifications[i].empty())
				{
					++spectrum_count;
					peptide_hit_count += identifications[i].getHits().size();
				}
			} 
			for(UInt i = 0; i < protein_identifications.size(); ++i)
			{
				++runs_count;
				protein_hit_count += protein_identifications[i].getHits().size();
			} 

			cout << "Number of spectra: " << spectrum_count << endl;
			cout << "Number of peptide hits: " << peptide_hit_count << endl;
			cout << "Number of runs: " << runs_count << endl;
			cout << "Number of protein hits: " << protein_hit_count << endl;
			
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPIdXMLInfo tool;
	return tool.main(argc,argv);
}
  
/// @endcond





