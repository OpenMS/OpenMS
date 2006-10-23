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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ContactPerson.h>

#include <qfileinfo.h>
#include <qfile.h>

#include "TOPPBase.h"

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
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPAnalysisXMLMerger
	: public TOPPBase
{
 public:
	TOPPAnalysisXMLMerger()
		: TOPPBase("AnalysisXMLMerger")
	{
			
	}
	
 protected:
	void printToolUsage_() const
	{
		cerr  << endl
					<< getToolName() << " -- Merges several analysisXML files into one analysisXML file." << endl
					<< "Version: " << VersionInfo::getVersion() << endl
					<< endl
					<< "Usage:" << endl
					<< " " << getToolName() << " [options]" << endl
					<< endl
					<< "Options are:" << endl
				 << "  -in               two ore more analysisXML files separated by comma (without blanks)" << endl
				 << "  -out              output analysisXML file name" << endl;
	}
	
	void printToolHelpOpt_() const
	{
		cerr << endl
				 << getToolName() << endl
				 << endl
				 << "INI options:" << endl
				 << "  in               two ore more analysisXML files separated by comma (without blanks)" << endl
				 << "  out              output analysisXML file name" << endl
				 << endl
				 << "INI File example section:" << endl
				 << "  <ITEM name=\"in\" value=\"file1.analysisXML,file2.analysisXML\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl;
	}
	
	void setOptionsAndFlags_()
	{
		options_["-in"] = "in";
		options_["-out"] = "out";
	}
	
	ExitCodes main_(int , char**)
	{
		vector<String> 									file_names;
		String 													actual_file_name;
		AnalysisXMLFile 								analysisXML_file;
		vector<ProteinIdentification> 	protein_identifications;
		vector<Identification> 					identifications;
		vector<Real> 										retention_times;		
		vector<Real> 										mz_values;		
		ContactPerson 									contact_person;
		vector<ProteinIdentification> 	additional_protein_identifications;
		vector<Identification> 					additional_identifications;
		vector<Real> 										additional_retention_times;		
		vector<Real> 										additional_mz_values;
		UnsignedInt											counter = 0;
		String 													out_file = "";
		String 													file_list	= "";	
		QFileInfo 											file_info;
		QFile 													file;


		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		//file list
		file_list = getParamAsString_("in");
		file_list.split(',', file_names);

		if (file_names.size() < 2)
		{
			writeLog_("Less than two filenames given. Aborting!");
			cout << "Less than two filenames given. Aborting!" << endl;
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
		
		writeDebug_(String("File list: ") + file_list, 1);

		//output file names and types
		out_file = getParamAsString_("out", "");
		if (out_file == "")
		{
			writeLog_("No output file specified. Aborting!");
			cout << "No output file specified. Aborting!" << endl;
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
		writeDebug_(String("Output file: ") + out_file, 1);
				
		//-------------------------------------------------------------
		// testing whether input and output files are accessible
		//-------------------------------------------------------------

		for(UnsignedInt i = 0; i < file_names.size(); ++i)
		{
			actual_file_name = file_names[i];
			file_info.setFile(actual_file_name.c_str());
			if (!file_info.exists())
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, actual_file_name);
			}
			if (!file_info.isReadable())
			{
				throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, actual_file_name);			
			}
	    if (file_info.size() == 0)
	    {
	      throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, actual_file_name);
	    }		
		}		
		file.setName(out_file.c_str());
		file.open( IO_WriteOnly );
		if (!file.isWritable())
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, out_file);
		}
		file.close();				

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		analysisXML_file.load(file_names[0],
													&protein_identifications,
													&identifications, 
													&retention_times,
													&mz_values,
													&contact_person);
		for(counter = 1; counter < file_names.size(); ++counter)
		{
			analysisXML_file.load(file_names[counter],
														&protein_identifications,
														&identifications, 
														&retention_times,
														&mz_values,
														&contact_person);
			protein_identifications.insert(protein_identifications.end(), additional_protein_identifications.begin(), additional_protein_identifications.end());
			identifications.insert(identifications.end(), additional_identifications.begin(), additional_identifications.end());
			retention_times.insert(retention_times.end(), additional_retention_times.begin(), additional_retention_times.end());
			mz_values.insert(mz_values.end(), additional_mz_values.begin(), additional_mz_values.end());						
		}										
															
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		analysisXML_file.store(out_file, 
													protein_identifications, 
													identifications, 
													retention_times, 
													mz_values);
			
		return OK;
	}
};


int main( int argc, char ** argv )
{
	TOPPAnalysisXMLMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
