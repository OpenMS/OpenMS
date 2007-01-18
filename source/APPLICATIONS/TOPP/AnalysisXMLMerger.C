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

#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/METADATA/Identification.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page AnalysisXMLMerger AnalysisXMLMerger
	
	@brief Merges several analysisXML files into one analysisXML file.
	
	You can merge an unlimited number of files into one analysisXML file. The file names
	that are to be merged are given at the '-in' parameter as a comma separated list.
	The output will be written to the file specified after the '-out' option.
	
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPAnalysisXMLMerger
	: public TOPPBase
{
 public:
	TOPPAnalysisXMLMerger()
		: TOPPBase("AnalysisXMLMerger","Merges several analysisXML files into one analysisXML file")
	{
			
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","two or more analysisXML files separated by comma (without blanks)");
		registerStringOption_("out","<file>","","output file in analysisXML format");
	}
	
	ExitCodes main_(int , char**)
	{
		vector<String> 									file_names;
		AnalysisXMLFile 								analysisXML_file;
		vector<ProteinIdentification> 	protein_identifications;
		vector<IdentificationData> 					identifications;
		vector<ProteinIdentification> 	additional_protein_identifications;
		vector<IdentificationData> 					additional_identifications;
		UnsignedInt											counter = 0;
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

		for(UnsignedInt i = 0; i < file_names.size(); ++i)
		{
			inputFileReadable_(file_names[i]);
		}
		
		outputFileWritable_(out_file);

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		analysisXML_file.load(file_names[0],
													protein_identifications,
													identifications);

		for(counter = 1; counter < file_names.size(); ++counter)
		{
			analysisXML_file.load(file_names[counter],
														additional_protein_identifications,
														additional_identifications);
			protein_identifications.insert(protein_identifications.end(), additional_protein_identifications.begin(), additional_protein_identifications.end());
			identifications.insert(identifications.end(), additional_identifications.begin(), additional_identifications.end());
		}										
															
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		analysisXML_file.store(out_file, 
													protein_identifications, 
													identifications);
			
		return EXECUTION_OK;
	}
};


int main( int argc, char ** argv )
{
	TOPPAnalysisXMLMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
