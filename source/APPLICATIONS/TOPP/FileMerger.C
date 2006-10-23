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

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FileMerger FileMerger
	
	@brief Merges several single scan files into an mzData file.
	
	Input is a text file with a list of file names and (optional) retention times.
	Output is a mzData file that contains the merged scans.
	
	Only the first scan out of each file and its meta information are copied into the output file.
	The meta information that is valid for all scans is taken from the first file.
	
	If the retention time is not given in the file list, it is taken from the corresponding input file.
	Alternatively the FileMerger numbers the scans consecutively starting from 1 if the flag '-auto_number' is given.
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileMerger
	: public TOPPBase
{
 public:
	TOPPFileMerger()
		: TOPPBase("FileMerger")
	{
			
	}
	
 protected:
	void printToolUsage_() const
	{
		cerr  << endl
					<< getToolName() << " -- Merges several single scan files into a mzData file." << endl
					<< "Version: " << VersionInfo::getVersion() << endl
					<< endl
					<< "Usage:" << endl
					<< " " << getToolName() << " [options]" << endl
					<< endl
					<< "Options are:" << endl
					<< "  -file_list <file> a text file containing file names and retention times sparated by tab." << endl
					<< "                    If no retention time is given, it is taken from the input file." << endl
					<< "  -in_type <type>   input file type (default: determined from input file extension)." << endl
					<< "  -out <file>       output mzData file name" << endl
					<< "  -auto_number      automatically numbers the scans (starting at 1)." << endl
					<< endl
					<< "Valid input types are: 'mzData', 'mzXML', 'DTA2D', 'ANDIMS' (cdf), 'DTA'" << endl;
	}
	
	void printToolHelpOpt_() const
	{
		cerr << endl
				 << getToolName() << endl
				 << endl
				 << "INI options:" << endl
				 << "  file_list <file> a text file containing file names and retention times sparated by tab" << endl
				 << "  in_type <type>   input file type (default: determined from input file extension)" << endl
				 << "  out              output mzData file name" << endl
				 << "  auto_number      automatically numbers the scans (starting at 1)." << endl
				 << endl
				 << "INI File example section:" << endl
				 << "  <ITEM name=\"file_list\" value=\"input.txt\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"in_type\" value=\"DTA\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"auto_number\" value=\"\" type=\"string\"/>" << endl;
	}
	
	void setOptionsAndFlags_()
	{
		options_["-file_list"] = "file_list";
		options_["-in_type"] = "in_type";
		options_["-out"] = "out";
		flags_["-auto_number"] = "auto_number";
	}
	
	ExitCodes main_(int , char**)
	{

		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		//file list
		String file_list_name = getParamAsString_("file_list");
		writeDebug_(String("File list: ") + file_list_name, 1);

		//file type
		FileHandler fh;
		FileHandler::Type force_type = fh.nameToType(getParamAsString_("in_type",""));
	
		//output file names and types
		String out_file = getParamAsString_("out");
		writeDebug_(String("Output file: ") + out_file, 1);

		//auto numbering
		bool auto_number = getParamAsBool_("auto_number",false);
			
		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------
		TextFile file_list(file_list_name);
			
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
			
		//parse filename and rt
		float rt,auto_rt=1;
		String line, filename;
		vector<String> tmp;
		MSExperiment<DPeak<1> > out, in;
		out.reserve(file_list.size());
		bool first_scan = true;
			
		for (TextFile::Iterator it = file_list.begin(); it!=file_list.end(); ++it)
		{
			line = *it;
			line.trim();
			if (line == "")
			{
				continue;
			}
				
			if (line.find("	")!=string::npos)
			{
				line.split('	',tmp);
				filename = tmp[0];
				if (tmp[1].size()!=0)
				{
					rt = tmp[1].toFloat();
				}
				else
				{
					rt = -1;
				}
			}
			else
			{
				filename = line;
				rt = -1;
			}
				
			//load file 
			fh.loadExperiment(filename,in,force_type);
			if (in.size()==0)
			{
				writeLog_(String("Warning: Empty file '") + filename +"'!");
				continue;
			}
			else if (in.size()>1)
			{
				writeLog_(String("Warning: More than one scan in file '") + filename +"'!");
			}
				
			//handle rt
			if (auto_number)
			{
				rt = auto_rt++;
			}
			if (in[0].getRetentionTime()==-1 || rt!=-1)
			{
				in[0].setRetentionTime(rt);
			}
				
			if (in[0].getRetentionTime()==-1)
			{
				writeLog_(String("No retention time for file '") + filename + "' given. Aborting!");
				printUsage_();
				return PARSE_ERROR;	
			}
				
			out.push_back(in[0]);
			// copy experimental settings from first file
			if (first_scan)
			{
				out.ExperimentalSettings::operator=(in);
				first_scan = false;
			}

		}
			
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		MzDataFile().store(out_file,out);
			
		return EXECUTION_OK;
	}
};


int main( int argc, char ** argv )
{
	TOPPFileMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
