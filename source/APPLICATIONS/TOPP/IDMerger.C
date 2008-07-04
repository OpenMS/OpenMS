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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page IDMerger IDMerger
	
	@brief Merges several IdXML files into one IdXML file.
	
	You can merge an unlimited number of files into one IdXML file. The file names
	that are to be merged are given at the '-in' parameter as a comma separated list.
	The output will be written to the file specified after the '-out' option.
	
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDMerger
	: public TOPPBase
{
 public:
	TOPPIDMerger()
		: TOPPBase("IDMerger","Merges several protein/peptide identification files into one file.")
	{
			
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<files>","","two or more IdXML files separated by comma (without blanks)");
		registerOutputFile_("out","<file>","","output file ");
		setValidFormats_("out",StringList::create("IdXML"));
	}
	
	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		//file list
		String file_list = getStringOption_("in");
		
		vector<String> file_names;
		file_list.split(',', file_names);
		if (file_names.size() < 2)
		{
			writeLog_("Less than two filenames given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}

		//output file names and types
		String out_file = getStringOption_("out");
				
		//-------------------------------------------------------------
		// testing whether input and output files are accessible
		//-------------------------------------------------------------

		for(UInt i = 0; i < file_names.size(); ++i)
		{
			inputFileReadable_(file_names[i]);
		}

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		IdXMLFile file;
		vector<ProteinIdentification> 	protein_identifications;
		vector<PeptideIdentification> identifications;
		vector<ProteinIdentification> 	additional_protein_identifications;
		vector<PeptideIdentification> additional_identifications;
		
		file.load(file_names[0], protein_identifications, identifications);

		vector<String> used_ids;
		for(UInt counter = 1; counter < file_names.size(); ++counter)
		{
			file.load(file_names[counter], additional_protein_identifications, additional_identifications);
			
			
			for (UInt i=0; i<additional_protein_identifications.size();++i)
			{
				if (find(used_ids.begin(), used_ids.end(), additional_protein_identifications[i].getIdentifier())!=used_ids.end())
				{
					writeLog_(String("Error: The idenitifier '") + additional_protein_identifications[i].getIdentifier() + "' was used before!");
					return INCOMPATIBLE_INPUT_DATA;
				}
				used_ids.push_back(additional_protein_identifications[i].getIdentifier());
			}
			
			protein_identifications.insert(protein_identifications.end(), additional_protein_identifications.begin(), additional_protein_identifications.end());
			identifications.insert(identifications.end(), additional_identifications.begin(), additional_identifications.end());
		}										
															
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		file.store(out_file, 
													protein_identifications, 
													identifications);
			
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPIDMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
