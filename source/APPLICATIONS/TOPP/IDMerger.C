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
		: TOPPBase("IDMerger","Merges several IdXML files into one IdXML file")
	{
			
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","two or more IdXML files separated by comma (without blanks)");
		registerStringOption_("out","<file>","","output file in IdXML format");
	}
	
	ExitCodes main_(int , char**)
	{
		vector<String> 									file_names;
		IdXMLFile 								IdXML_file;
		vector<Identification> 	protein_identifications;
		vector<PeptideIdentification> 					identifications;
		vector<Identification> 	additional_protein_identifications;
		vector<PeptideIdentification> 					additional_identifications;
		UInt											counter = 0;
		String 													out_file = "";
		String 													file_list	= "";


		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		//file list
		file_list = getStringOption_("in");
		file_list.split(',', file_names);
		if (file_names.size() < 2)
		{
			writeLog_("Less than two filenames given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}

		//output file names and types
		out_file = getStringOption_("out");
				
		//-------------------------------------------------------------
		// testing whether input and output files are accessible
		//-------------------------------------------------------------

		for(UInt i = 0; i < file_names.size(); ++i)
		{
			inputFileReadable_(file_names[i]);
		}
		
		outputFileWritable_(out_file);

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		IdXML_file.load(file_names[0],
													protein_identifications,
													identifications);

		for(counter = 1; counter < file_names.size(); ++counter)
		{
			IdXML_file.load(file_names[counter],
														additional_protein_identifications,
														additional_identifications);
			protein_identifications.insert(protein_identifications.end(), additional_protein_identifications.begin(), additional_protein_identifications.end());
			identifications.insert(identifications.end(), additional_identifications.begin(), additional_identifications.end());
		}										
															
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		IdXML_file.store(out_file, 
													protein_identifications, 
													identifications);
			
		return EXECUTION_OK;
	}
};


int main( int argc, char ** argv )
{
	TOPPIDMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
