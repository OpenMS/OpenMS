// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_TextImporter TextImporter
	
	@brief This application converts text files to %OpenMS XML formats.
	
	Currently only featureXML can we written.
	
	@todo Add import of msInspect and SpecArray feature data (Marc)
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_TextImporter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{

  class TOPPTextImporter 
  	: public TOPPBase
  {
    public:
      TOPPTextImporter() :
        TOPPBase("TextImporter", "Imports text files and converts them to XML.")
      {
      }

    protected:

      void registerOptionsAndFlags_()
      {
				registerInputFile_("in", "<file>", "", "(Excel readable) Text file (supported formats: see below)");
				registerInputFile_ ("template_ini", "<file>", "", "Template Ini file to augment", false);
        registerOutputFile_("out", "<file>", "", "Output XML file.");
        setValidFormats_("out",StringList::create( "featureXML,ini"));
				registerStringOption_("out_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
				setValidStrings_("out_type",StringList::create("featureXML,ini"));
        registerStringOption_( "separator", "<sep>", "", "The used separator characters in the input. If unset the 'tab' character is used.", false);
				addEmptyLine_();
				addText_("The following conversions are supported:");
				addText_("- CSV to featureXML");
				addText_("    Input text file containing the following columns: RT, m/z, intensity.");
				addText_("    Additionally meta data columns may follow.");
				addText_("    If meta data is used, meta data column names have to be specified in a header line.");
				addText_("- CSV to INI(ITRAQAnalyzer-settings)");
				addText_("    Input text file contains meta data as well as isotope correction matrix");
				addText_("    and channel assignments. The -template_ini option is mandatory and serves as template for the output ini file.");
				//addText_("    For the template CSV see share/OpenMS/");

      }

      ExitCodes main_( int, const char** )
      {
        //-------------------------------------------------------------
        // parameter handling
        //-------------------------------------------------------------
        String in = getStringOption_("in");
        String out = getStringOption_("out");

				FileHandler fh;
				FileTypes::Type out_type = fh.nameToType(getStringOption_("out_type"));
				if (out_type==FileTypes::UNKNOWN)
				{
					out_type = fh.getType(out);
					writeDebug_(String("Output file type: ") + fh.typeToName(out_type), 2);
				}

				if (out_type==FileTypes::UNKNOWN)
				{
					writeLog_("Error: Could not determine output file type!");
					return PARSE_ERROR;
				}

        String separator = getStringOption_("separator");
        if ( separator == "" ) separator = "\t";

        //-------------------------------------------------------------
        // load input
        //-------------------------------------------------------------
				TextFile text(in);



        //-------------------------------------------------------------
        // processing
        //-------------------------------------------------------------
				if (out_type==FileTypes::FEATUREXML)
				{
					//-------------------------------------------------------------
					// parsing header line
					//-------------------------------------------------------------
					vector<String> headers;
					text[0].split(separator[0], headers);
					int offset = 0;
					for (Size i=0; i<headers.size(); ++i)
					{
						headers[i].trim();
					}
					String header_trimmed = text[0];
					header_trimmed.trim();
					DoubleReal rt = 0.0;
					DoubleReal mz = 0.0;
					DoubleReal it = 0.0;
					// see if we have a header
					try
					{
						rt = headers[0].toDouble();
						mz = headers[1].toDouble();
						it = headers[2].toDouble();
					}
					catch (Exception::BaseException&)
					{
						offset=1;
						std::cout << "Detected a header line.\n";
					}
					//-------------------------------------------------------------
					// parsing features
					//-------------------------------------------------------------
					FeatureMap<> feature_map;
					feature_map.reserve(text.size());
					for (Size i=offset; i<text.size(); ++i)
					{
						//do nothing for empty lines
						String line_trimmed = text[i];
						line_trimmed.trim();
						if (line_trimmed=="")
						{
							if (i<text.size()-1) writeLog_(String("Notice: Empty line ignored (line ") + (i+1) + ").");
							continue;
						}
						
						//split line to tokens
						vector<String> parts;
						text[i].split(separator[0], parts);
						
						//abort if line does not contain enough fields
						if (parts.size()<3)
						{
							writeLog_("Error: Invalid input line: At least three columns are needed!");
							writeLog_(String("Offending line: '") + line_trimmed + "'  (line " + (i+1) + ")");
							return INPUT_FILE_CORRUPT;
						}
						
						//convert coordinate columns to doubles
						try
						{
							rt = parts[0].toDouble();
							mz = parts[1].toDouble();
							it = parts[2].toDouble();
						}
						catch (Exception::BaseException&)
						{
							writeLog_("Error: Invalid input line: Could not convert the first three columns to float!");
							writeLog_("       Is the correct separator specified?");
							writeLog_(String("Offending line: '") + line_trimmed + "'  (line " + (i+1) + ")");
							return INPUT_FILE_CORRUPT;
						}
						Feature f;
						f.setMZ(mz);
						f.setRT(rt);
						f.setIntensity(it);

						//parse meta data
						for (Size j=3; j<parts.size(); ++j)
						{
							String part_trimmed = parts[j];
							part_trimmed.trim();
							if (part_trimmed!="")
							{
								//check if column name is ok
								if (headers.size()<=j || headers[j]=="")
								{
									writeLog_(String("Error: Missing meta data header for column ") + (j+i) + "!");
									writeLog_(String("Offending header line: '") + header_trimmed + "'  (line 1)");
									return INPUT_FILE_CORRUPT;
								}
								//add meta value
								f.setMetaValue(headers[j],part_trimmed);
							}

						}
						
						//insert feature to map
						feature_map.push_back(f);
					}
				
					//-------------------------------------------------------------
					// write output
					//-------------------------------------------------------------
					
					//annotate output with data processing info
					addDataProcessing_(feature_map, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));
					
					FeatureXMLFile().store(out, feature_map);
				}
				else // PARAM
				{
					Param p;
					String ini_file("");
					ini_file = getStringOption_("template_ini");
					if (File::exists(ini_file))	p.load(ini_file);
					else
					{
						std::cerr << "For INI file output this tool requires a template ini file to augment. Please use the -template_ini argument!\n";
						exit(MISSING_PARAMETERS);
					}

					enum mode {ITRAQ_METADATA, ITRAQ_CHANNELALLOC, ITRAQ_MATRIX};
					int imode = ITRAQ_METADATA;
					StringList channel_alloc;
					StringList isotope_matrix;
					
					// get current instance and assume its the one we want to change
					StringList subs;
					getIniLocation_().split(':',subs,false);
					String instance = subs[1];

					for (Size i=0; i<text.size(); ++i)
					{
						//do nothing for empty lines
						String line_trimmed = text[i];
						line_trimmed.trim();

						//split line to tokens
						vector<String> parts;
						text[i].split(separator[0], parts, true);

						if (line_trimmed=="" || parts[0].hasPrefix("**COMMENT") || parts[0].trim()=="")
						{
							if (i<text.size()-1) writeLog_(String("Notice: Empty/Comment line ignored (line ") + (i+1) + ").");
							continue;
						}

						if (parts[0].has(':') || parts[1].has(':'))
						{
							writeLog_(String("Invalid character ':' found in line ") + (i+1) + String(". Aborting."));
							exit(INPUT_FILE_CORRUPT);
						}
						
						if (parts[0].hasPrefix("**METADATA")) imode = ITRAQ_METADATA;
						else if (parts[0].hasPrefix("**ITRAQ [CHANNELALLOC]")) imode = ITRAQ_CHANNELALLOC;
						else if (parts[0].hasPrefix("**ITRAQ [ISOTOPE_4PLEX_CORRECTION]")) imode = ITRAQ_MATRIX;
						else
						{ // actual content

							switch (imode)
							{
								case ITRAQ_METADATA:
									// add to meta section
									p.setValue("ITRAQAnalyzer:"+instance+":algorithm:MetaInformation:"+parts[0],parts[1].trim(), "MetaValue" ,StringList::create("advanced"));
									break;
								case ITRAQ_CHANNELALLOC:
									if (parts[1].trim()=="") break;
									parts[0].split(' ', subs);
									try{subs[1].toInt();}
									catch (...)
									{
										writeLog_(String("Channel allocation entry in column 1 in line ") + String(i+1) + String(" does not have the format <String> <Number> <String> in CSV file! Terminating..."));
										exit(INCOMPATIBLE_INPUT_DATA);
									}
									channel_alloc.push_back(subs[1] + ":" + parts[1].trim());
									break;

								case ITRAQ_MATRIX:
									// determine channel
									Int channel=parts[0].substr(7,3).toInt();
									// is something filled in?
									if (parts.size()<5)
									{
										writeLog_(String("CSV file does not have enough matrix correction entries for channel ")+String(channel)+String("! Terminating..."));
										exit(INCOMPATIBLE_INPUT_DATA);
									}
									for (Int i=1;i<5;++i)
									{
										try{parts[i].toDouble();}
										catch (...)
										{
											writeLog_(String("Correction matrix entry #") + String(i) + String(" for channel ")+String(channel)+String(" in CSV file is not a number or missing! Terminating..."));
											exit(INCOMPATIBLE_INPUT_DATA);
										}
									}
									// create string
									isotope_matrix.push_back(String(channel)+":"+parts[1]+"/"+parts[2]+"/"+parts[3]+"/"+parts[4]);
									// result: 114:0/1/5.9/0.2
									break;
							}
						}

					}
					
					if (channel_alloc.size()==0)
					{
						writeLog_(String("CSV file does not contain compulsory channel allocation information!"));
						exit(INCOMPATIBLE_INPUT_DATA);
					}

					p.setValue("ITRAQAnalyzer:"+instance+":algorithm:Extraction:channel_active", 
										 channel_alloc,
										 p.getDescription("ITRAQAnalyzer:"+instance+":algorithm:Extraction:channel_active"),
										 p.getTags("ITRAQAnalyzer:"+instance+":algorithm:Extraction:channel_active"));

					if (isotope_matrix.size()!=4)
					{
						writeLog_(String("CSV file does not contain complete isotope correction matrix! Terminating..."));
						exit(INCOMPATIBLE_INPUT_DATA);
					}
					p.setValue("ITRAQAnalyzer:"+instance+":algorithm:Quantification:isotope_correction_values", 
										 isotope_matrix,
										 p.getDescription("ITRAQAnalyzer:"+instance+":algorithm:Quantification:isotope_correction_values"),
										 p.getTags("ITRAQAnalyzer:"+instance+":algorithm:Quantification:isotope_correction_values"));

					// store result
					p.store(out);
				}

        return EXECUTION_OK;
      }
  };
}

int main( int argc, const char** argv )
{
  TOPPTextImporter t;
  return t.main(argc, argv);
}

/// @endcond
