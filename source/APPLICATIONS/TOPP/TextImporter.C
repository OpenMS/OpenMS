// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/CONCEPT/Constants.h>

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
				setValidStrings_("mode",StringList::create("default,msInspect,SpecArray,Kroenik"));
				addEmptyLine_();
				addText_("The following conversion modes are supported:");
				addText_("- default");
				addText_("    Input text file containing the following columns: RT, m/z, intensity.");
				addText_("    Additionally meta data columns may follow.");
				addText_("    If meta data is used, meta data column names have to be specified in a header line.");
				addText_("    If a meta column named 'charge' with numeric data exists, the charge of the features will be set accordingly.");
				addText_("- msInspect");
				addText_("    Imports an msInspect feature file.");
				addText_("- SpecArray");
				addText_("    Imports a SpecArray feature file.");
				addText_("- Kroenik");
				addText_("    Imports a Kroenik (Hardkloer sibling) feature file.");
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
            if (headers.size()>3) throw Exception::BaseException(); // there is meta-data, so these must be their names
            else if (headers.size()<3) throw Exception::BaseException(); // not enough data columns in first line...
            // try to convert... if not: thats a header
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
								if (headers[j] == "charge")
								{
									try
									{
										f.setCharge(part_trimmed.toInt());
									}
									catch (...)
									{
										writeLog_(String("Failed to convert metavalue 'charge' into integer (line '") + (i+1) + ")");
									}
								}
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
						Size column_to_convert=0;
						try
						{
							column_to_convert = 1;
							f.setRT(parts[1].toDouble());
							column_to_convert = 2;
							f.setMZ(parts[2].toDouble());
							column_to_convert = 5;
							f.setIntensity(parts[5].toDouble());
							column_to_convert = 6;
							f.setCharge(parts[6].toInt());
							column_to_convert = 8;
							f.setOverallQuality(parts[8].toDouble());

							column_to_convert = 3;
							f.setMetaValue("accurateMZ",parts[3]);
							column_to_convert = 4;
							f.setMetaValue("mass",parts[4].toDouble());
							column_to_convert = 7;
							f.setMetaValue("chargeStates",parts[7].toInt());
							column_to_convert = 9;
							f.setMetaValue("background",parts[9].toDouble());
							column_to_convert = 10;
							f.setMetaValue("median",parts[10].toDouble());
							column_to_convert = 11;
							f.setMetaValue("peaks",parts[11].toInt());
							column_to_convert = 12;
							f.setMetaValue("scanFirst",parts[12].toInt());
							column_to_convert = 13;
							f.setMetaValue("scanLast",parts[13].toInt());
							column_to_convert = 14;
							f.setMetaValue("scanCount",parts[14].toInt());
							column_to_convert = 15;
							f.setMetaValue("totalIntensity",parts[15].toDouble());
							column_to_convert = 16;
							f.setMetaValue("sumSquaresDist",parts[16].toDouble());
						}
						catch (Exception::BaseException /*&e*/)
						{
							writeLog_(String("Failed to convert value in column ") + String(column_to_convert+1) + "into a number (line '" + (i+1) + ")");
						}
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
						try
						{						
							f.setMZ(line.substr(0,12).toDouble());
							f.setCharge(line.substr(36,12).toInt());
							f.setRT(line.substr(12,12).toDouble() *60.0);
							f.setIntensity(line.substr(48,12).toDouble());
							f.setMetaValue("s/n",line.substr(24,12).toDouble());
						}
						catch (Exception::BaseException /*&e*/)
						{
							writeLog_(String("Failed to convert value into a number (line '") + (i+1) + ")");
						}
						feature_map.push_back(f);
					}
				}
        //-------------------------------------------------------------
        // HardKlör
        //-------------------------------------------------------------
				else if (mode=="Kroenik")
				{
					for (Size i=1; i<input.size(); ++i)
					{
						String line = input[i];
					
						//split lines: File,	First Scan,	Last Scan,	Num of Scans,	Charge,	Monoisotopic Mass, Base Isotope Peak,	Best Intensity,	Summed Intensity,	First RTime,	Last RTime,	Best RTime,	Best Correlation,	Modifications
						std::vector< String > parts;
						line.split('\t', parts);
						
						if (parts.size() != 14)
						{
							std::cerr << "Line #" << (i+1) << " does not have the expected 14 tab-separated entries. Skipping this line!\n";
							continue;
						}
						//create feature
						Feature f;
						f.setCharge(parts[4].toInt());
						f.setMZ(parts[5].toDouble()/f.getCharge() + Constants::PROTON_MASS_U);
						f.setRT(parts[11].toDouble());
						f.setOverallQuality(parts[12].toDouble());
						f.setIntensity(parts[8].toDouble());
						ConvexHull2D hull;
						ConvexHull2D::PointType point;
						
						point.setX(parts[9].toDouble());
						point.setY(f.getMZ());
						hull.addPoint(point);

						point.setX(parts[9].toDouble());
						point.setY(f.getMZ()+3.0/(DoubleReal)f.getCharge());
						hull.addPoint(point);

						point.setX(parts[10].toDouble());
						point.setY(f.getMZ()+3.0/(DoubleReal)f.getCharge());
						hull.addPoint(point);

						point.setX(parts[10].toDouble());
						point.setY(f.getMZ());
						hull.addPoint(point);
						
						point.setX(parts[9].toDouble());
						point.setY(f.getMZ());
						hull.addPoint(point);
						
						std::vector< ConvexHull2D > hulls;
						hulls.push_back(hull);
						f.setConvexHulls(hulls);
						f.setMetaValue("Mass",parts[5].toDouble());
						f.setMetaValue("FirstScan",parts[1].toDouble());
						f.setMetaValue("LastScan",parts[2].toInt());
						f.setMetaValue("NumOfScans",parts[3].toDouble());
						f.setMetaValue("AveragineModifications",parts[13]);
						feature_map.push_back(f);
					}

					std::cout << "Hint: The convex hulls are approximated in m/z dimension (Kroenik lacks this information)!\n";
				}
				
				std::cout << "Converted " << feature_map.size() << " features!\n";
				
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
