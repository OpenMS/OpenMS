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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FileMerger FileMerger
	
	@brief Merges several files into an mzData file.
	
	Input is a text file containing a list of file names and (optional) retention times.
	Output is a mzData file that contains the merged scans.
	
	The meta information that is valid for the whole experiment is taken from the first file.
	
	The retention times for the individual scans are taken from the input file(s) meta data,
	from the input file names or are generated.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileMerger
	: public TOPPBase
{
 public:
	TOPPFileMerger()
		: TOPPBase("FileMerger","Merges several MS files into one file.")
	{
			
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("file_list","<file>","","a text file containing one input file name per line");		
		registerStringOption_("in_type","<type>","","input file type (default: determined from file extension or content)\n", false);
		setValidStrings_("in_type",StringList::create("mzData,mzXML,DTA,DTA2D,cdf,mgf"));
		registerOutputFile_("out","<file>","","output file ");
		setValidFormats_("out",StringList::create("mzData"));
		
		registerFlag_("rt_auto","Assign retention times automatically (integers starting at 1)");
		registerFlag_("rt_file","Take retention times from file_list.\n"
														"If this flag is activated, the file list has to contain a filename and a\n"
														"retention time separated by tab in each line.");
		registerFlag_("rt_from_filename", "If this flag is set FileMerger tries to guess the rt of the spectrum.\n"
																			"This option is useful for merging DTA file, which should contain the string\n"
																			"'rt' directly followed by a floating point number:\n"
																			"i.e. my_spectrum_rt2795.15.dta"); 
		registerIntOption_("ms_level", "<num>", 2, "this option is useful for use with DTA files which does not \n"
																								"contain MS level information. The given level is assigned to the spectra.", false);
		registerFlag_("user_ms_level", "If this flag is set, the MS level given above is used");
		addEmptyLine_();
		addText_("Note: Meta data about the whole experiment is taken from the first file in the list!");
	}
	
	ExitCodes main_(int , const char**)
	{

		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		//file list
		String file_list_name = getStringOption_("file_list");

		//file type
		FileHandler fh;
		FileHandler::Type force_type = fh.nameToType(getStringOption_("in_type"));
	
		//output file names and types
		String out_file = getStringOption_("out");

		//auto numbering
		bool auto_number = getFlag_("rt_auto");
		bool rt_from_file = getFlag_("rt_file");
		bool user_ms_level = getFlag_("user_ms_level");
		bool rt_from_filename = getFlag_("rt_from_filename");
			
		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------
		TextFile file_list(file_list_name);
			
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		
		float rt_final,rt_file,rt_auto=0;
		String line, filename;
		vector<String> tmp;
		MSExperiment<Peak1D> out, in;
		out.reserve(file_list.size());
		bool first_file = true;
		
		for (TextFile::Iterator it = file_list.begin(); it!=file_list.end(); ++it)
		{
			//parsing of file list
			rt_file = -1;
			line = *it;
			line.trim();
			if (line == "")
			{
				continue;
			}
			
			//tab separator
			if (line.find("	")!=string::npos)
			{
				line.split('	',tmp);
				filename = tmp[0];
				if (tmp[1].size()!=0)
				{
					rt_file = tmp[1].toFloat();
				}
			}
			// space separator
			else if (line.find(" ")!=string::npos)
			{
				line.split(' ',tmp);
				filename = tmp[0];
				if (tmp[1].size()!=0)
				{
					rt_file = tmp[1].toFloat();
				}
			}
			else
			{
				filename = line;
			}
			
			//load file 
			fh.loadExperiment(filename,in,force_type,log_type_);
			if (in.size()==0)
			{
				writeLog_(String("Warning: Empty file '") + filename +"'!");
				continue;
			}
			else if (in.size()>1)
			{
				out.reserve(out.size()+in.size());
				if (rt_from_file)
				{
					writeLog_(String("Warning: More than one scan in file '") + filename +"'! All scans will have the same retention time!");
				}
			}
			
			for (MSExperiment<Peak1D>::const_iterator it2 = in.begin(); it2!=in.end(); ++it2)
			{ 
				//handle rt
				++rt_auto;
				rt_final = -1;
				if (auto_number)
				{
					rt_final = rt_auto;
				}
				else if (rt_from_file) 
				{
					rt_final = rt_file;
				}
				else
				{
					rt_final = it2->getRT();
				}
	
				// guess the retention time from filename
				if (rt_from_filename)
				{
					if (!filename.hasSubstring("rt"))
					{
						writeLog_(String("Warning: cannot guess retention time from filename as it does not contain 'rt'"));
					}
					for (UInt i = 0; i < filename.size(); ++i)
					{
						if (filename[i] == 'r' && ++i != filename.size() && filename[i] == 't' && ++i != filename.size() && isdigit(filename[i]))
						{
							String rt;
							while (i != filename.size() && (filename[i] == '.' || isdigit(filename[i])))
							{
								rt += filename[i++];
							}
							if (rt.size() > 0)
							{
								// remove dot from rt3892.98.dta
								//                          ^
								if (rt[rt.size() - 1] == '.')
								{
									// remove last character
									rt.erase(rt.end() - 1);
								}
							}
							try 
							{
								float tmp = rt.toFloat();
								rt_final = tmp;
							}
							catch (Exception::ConversionError)
							{
								 writeLog_(String("Warning: cannot convert the found retention time in a value '" + rt + "'."));
							}
						}
					}
				}

				// none of the rt methods were successful
        if(rt_final == -1)
				{
					writeLog_(String("Warning: No valid retention time for output scan '") + rt_auto +"' from file '" + filename + "'");
				}
				
				out.push_back(*it2);
				out.back().setRT(rt_final);
				if (user_ms_level)
				{
					out.back().setMSLevel((int)getIntOption_("ms_level"));
				}
			}

			// copy experimental settings from first file
			if (first_file)
			{
				out.ExperimentalSettings::operator=(in);
				first_file = false;
			}
		}
			
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		MzDataFile f;
		f.setLogType(log_type_);
		f.store(out_file,out);
			
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPFileMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
