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
// $Maintainer: $
// $Authors: Marc Sturm $
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
        registerOutputFile_("out", "<file>", "", "Output XML file.");
        setValidFormats_("out",StringList::create( "featureXML"));
        registerStringOption_( "separator", "<sep>", "", "The used separator characters in the input. If unset the 'tab' character is used.", false);
				registerStringOption_( "mode", "<mode>", "default", "Conversion mode (see below).", false);
				setValidStrings_("mode",StringList::create("default,msInspect,SpecArray"));
				addEmptyLine_();
				addText_("The following conversion modes are supported:");
				addText_("- default");
				addText_("    Input text file containing the following columns: RT, m/z, intensity.");
				addText_("    Additionally meta data columns may follow.");
				addText_("    If meta data is used, meta data column names have to be specified in a header line.");
				addText_("- msInspect");
				addText_("    Imports an msInspect feature file.");
				addText_("- SpecArray");
				addText_("    Imports a SpecArray feature file.");
      }

      ExitCodes main_( int, const char** )
      {
        // parameter handling
        String in = getStringOption_("in");
        String out = getStringOption_("out");
        String mode = getStringOption_("mode");

        String separator = getStringOption_("separator");
        if ( separator == "" ) separator = "\t";

        // load input
				TextFile input(in);
				
				// init output
				FeatureMap<> feature_map;

        //-------------------------------------------------------------
        // default
        //-------------------------------------------------------------
				if (mode=="default")
				{
					// parsing header line
					vector<String> headers;
					input[0].split(separator[0], headers);
					int offset = 0;
					for (Size i=0; i<headers.size(); ++i)
					{
						headers[i].trim();
					}
					String header_trimmed = input[0];
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

					// parsing features
					feature_map.reserve(input.size());
					for (Size i=offset; i<input.size(); ++i)
					{
						//do nothing for empty lines
						String line_trimmed = input[i];
						line_trimmed.trim();
						if (line_trimmed=="")
						{
							if (i<input.size()-1) writeLog_(String("Notice: Empty line ignored (line ") + (i+1) + ").");
							continue;
						}
						
						//split line to tokens
						vector<String> parts;
						input[i].split(separator[0], parts);
						
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
				}
        //-------------------------------------------------------------
        // msInspect
        //-------------------------------------------------------------
				else if (mode=="msInspect")
				{
					bool first_line = true;
					for (Size i=1; i<input.size(); ++i)
					{
						String line = input[i];
						
						//ignore comment lines
						if (line.empty() || line[0]=='#') continue;
				
						//skip leader line
						if (first_line)
						{
							first_line = false;
							continue;
						}
						
						//split lines: scan	time	mz	accurateMZ	mass	intensity	charge	chargeStates	kl	background	median	peaks	scanFirst	scanLast	scanCount	totalIntensity	sumSquaresDist	description
						std::vector< String > parts;
						line.split('\t', parts);
						
						//create feature
						Feature f;
						f.setMZ(parts[2].toDouble());
						f.setCharge(parts[6].toInt());
						f.setRT(parts[1].toDouble());
						f.setOverallQuality(parts[8].toDouble());
						f.setIntensity(parts[5].toDouble());
						f.setMetaValue("accurateMZ",parts[3]);
						f.setMetaValue("mass",parts[4].toDouble());
						f.setMetaValue("chargeStates",parts[7].toInt());
						f.setMetaValue("background",parts[9].toDouble());
						f.setMetaValue("median",parts[10].toDouble());
						f.setMetaValue("peaks",parts[11].toInt());
						f.setMetaValue("scanFirst",parts[12].toInt());
						f.setMetaValue("scanLast",parts[13].toInt());
						f.setMetaValue("scanCount",parts[14].toInt());
						f.setMetaValue("totalIntensity",parts[15].toDouble());
						f.setMetaValue("sumSquaresDist",parts[16].toDouble());
						f.setMetaValue("description",parts[17]);
						feature_map.push_back(f);
					}
				}
        //-------------------------------------------------------------
        // SpecArray
        //-------------------------------------------------------------
				else if (mode=="SpecArray")
				{
					for (Size i=1; i<input.size(); ++i)
					{
						String line = input[i];
						
						Feature f;
						f.setMZ(line.substr(0,12).toDouble());
						f.setCharge(line.substr(36,12).toInt());
						f.setRT(line.substr(12,12).toDouble() *60.0);
						f.setIntensity(line.substr(48,12).toDouble());
						f.setMetaValue("s/n",line.substr(24,12).toDouble());
						feature_map.push_back(f);
					}
				}
				
				// assign unique ids
				feature_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

				//annotate output with data processing info
				addDataProcessing_(feature_map, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

				// write output
				FeatureXMLFile().store(out, feature_map);

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
