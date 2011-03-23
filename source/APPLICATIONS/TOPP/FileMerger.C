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
// $Maintainer: Andreas Bertsch, Chris Bielow $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FileMerger FileMerger

	@brief Merges several files into an mzML file.
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FileMerger \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool/instrument producing mergeable files </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating merged files (e.g. @ref TOPP_XTandemAdapter) </td>
		</tr>
	</table>
</CENTER>

	The meta information that is valid for the whole experiment (e.g. MS instrument and sample)
	is taken from the first file.

	The retention times for the individual scans are taken from either:
  <ul>
  <li>the input file meta data (e.g. mzML)
  <li>from the input file names (name must contain 'rt' directly followed by a number, e.g. 'myscan_rt3892.98_MS2.dta')
  <li>as a list (one RT for each file)
  <li>or are auto-generated (starting at 1 with 1 second increment).
  </ul>

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FileMerger.cli
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
		registerInputFileList_("in","<files>",StringList(),"Input files separated by blank");
#ifdef USE_ANDIMS
		setValidFormats_("in",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,cdf,mgf,featureXML,fid"));
#else
		setValidFormats_("in",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,mgf,featureXML,fid"));
#endif
		registerStringOption_("in_type","<type>","","input file type (default: determined from file extension or content)\n", false);
#ifdef USE_ANDIMS
		setValidStrings_("in_type",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,cdf,mgf,featureXML,fid"));
#else
		setValidStrings_("in_type",StringList::create("mzData,mzXML,mzML,DTA,DTA2D,mgf,featureXML,fid"));
#endif
		registerOutputFile_("out","<file>","","output file");
		setValidFormats_("out",StringList::create("mzML,featureXML"));

    registerFlag_("annotate_file_origin","Store the original filename in each feature (MetaValue: file_origin).");

		addEmptyLine_();
		addText_ ("Flags for non-FeatureXML input/output:");
		registerFlag_("rt_auto","Assign retention times automatically (integers starting at 1)");
		registerDoubleList_("rt_custom","<rt>",DoubleList(),"List of custom retention times that are assigned to the files.\n"
		                                "The number of given retention times must be equal to the number of given input file.", false);
		registerFlag_("rt_filename", "If this flag is set FileMerger tries to guess the rt of the file name.\n"
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
		StringList file_list = getStringList_("in");

		//file type
		FileHandler fh;
		FileTypes::Type force_type;
		if (getStringOption_("in_type").size()>0)
		{ 
			force_type= fh.nameToType(getStringOption_("in_type"));
		}
		else
		{
			force_type= fh.getType(file_list[0]);
		}

		//output file names and types
		String out_file = getStringOption_("out");

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------

    bool annotate_file_origin =  getFlag_("annotate_file_origin");

		if (force_type == FileTypes::FEATUREXML)
		{
			FeatureMap<> out;
      for (Size i = 0; i < file_list.size(); ++i)
			{
				FeatureMap<> map;
				FeatureXMLFile fh;
				fh.load(file_list[i], map);

        if(annotate_file_origin)
        {
          for(FeatureMap<>::iterator it = map.begin(); it != map.end(); ++it)
          {
            it->setMetaValue("file_origin", DataValue(file_list[i]));
          }
        }
				out += map;
			}

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------

			//annotate output with data processing info
			addDataProcessing_(out, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

			FeatureXMLFile f;
			f.store(out_file,out);
						
		}
		else
		{
		  // we might want to combine different types, thus we only
		  // query in_type (which applies to all files)
		  // and not the suffix or content of a single file
			force_type = fh.nameToType(getStringOption_("in_type"));

			//rt
			bool rt_auto_number = getFlag_("rt_auto");
			bool rt_filename = getFlag_("rt_filename");
			bool rt_custom = false;
			DoubleList custom_rts = getDoubleList_("rt_custom");
			if (custom_rts.size()!=0)
			{
				rt_custom = true;
				if (custom_rts.size()!=file_list.size())
				{
					writeLog_("Custom retention time list must have as many elements as there are input files!");
					printUsage_();
					return ILLEGAL_PARAMETERS;
				}
			}

			//ms level
			bool user_ms_level = getFlag_("user_ms_level");

			
			MSExperiment<> out;
			out.reserve(file_list.size());
			UInt rt_auto = 0;
			UInt native_id = 0;
			for (Size i = 0; i < file_list.size();++i)
			{
				String filename = file_list[i];

				//load file
				MSExperiment<> in;
				fh.loadExperiment(filename,in,force_type,log_type_);
				if (in.size()==0)
				{
					writeLog_(String("Warning: Empty file '") + filename +"'!");
					continue;
				}
				out.reserve(out.size()+in.size());

				//warn if custom RT and more than one scan in input file
				if (rt_custom && in.size()>1)
				{
					writeLog_(String("Warning: More than one scan in file '") + filename +"'! All scans will have the same retention time!");
				}

				for (MSExperiment<>::const_iterator it2 = in.begin(); it2!=in.end(); ++it2)
				{
					//handle rt
					Real rt_final = it2->getRT();
					if (rt_auto_number)
					{
						rt_final = ++rt_auto;
					}
					else if (rt_custom)
					{
						rt_final = custom_rts[i];
					}
					else if (rt_filename)
					{
						if (!filename.hasSubstring("rt"))
						{
							writeLog_(String("Warning: cannot guess retention time from filename as it does not contain 'rt'"));
						}
						for (Size i = 0; i < filename.size(); ++i)
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
					out.back().setNativeID(native_id);
					if (user_ms_level)
					{
						out.back().setMSLevel((int)getIntOption_("ms_level"));
					}
					++native_id;
				}

				// copy experimental settings from first file
				if (i==0)
				{
					out.ExperimentalSettings::operator=(in);
				}
			}

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------

			//annotate output with data processing info
			addDataProcessing_(out, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

			MzMLFile f;
			f.setLogType(log_type_);
			f.store(out_file,out);

		}
		

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPFileMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
