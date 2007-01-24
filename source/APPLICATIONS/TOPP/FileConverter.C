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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>


#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FileConverter FileConverter
	
	@brief Converts between different MS file formats.
	
	Supported input file types are: 'mzData', 'mzXML', 'DTA2D', 'ANDIMS' (cdf).<BR>
	'FeatureFile' (OpenMS features) is also supported but will lose feature specific information.
	
	Supported output file types are: 'mzData', 'mzXML', 'DTA2D'
	'FeatureFile' can be generated using defaults for feature specific information.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileConverter
	: public TOPPBase
{
 public:
	TOPPFileConverter()
		: TOPPBase("FileConverter","converts between different MS file formats")
	{
			
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","input file");
		registerStringOption_("in_type", "<type>", "",
													"input file type (default: determined from file extension or content)\n"
													"Valid input types are: 'mzData', 'mzXML', 'DTA2D', 'ANDIMS'.\n"
													"'FeatureFile' can be converted, but will lose feature specific information", false);
		registerStringOption_("out","<file>","","output file");
		registerStringOption_("out_type", "<type>", "",
													"output file type (default: determined from output file extension)\n"
													"Valid output types are: 'mzData', 'mzXML', 'DTA2D'.\n"
													"'FeatureFile' can be generated using defaults for feature specific information", false);
	}
	
	ExitCodes main_(int , char**)
	{
	
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		//input file names
		String in = getStringOption_("in");
		inputFileReadable_(in);	
			
		//input file type
		FileHandler fh;
		FileHandler::Type in_type = fh.nameToType(getStringOption_("in_type"));
					
		if (in_type==FileHandler::UNKNOWN)
		{
			in_type = fh.getTypeByFileName(in);
			writeDebug_(String("Input file type (from file extention): ") + fh.typeToName(in_type), 2);
		}

		if (in_type==FileHandler::UNKNOWN)
		{
			in_type = fh.getTypeByContent(in);
			writeDebug_(String("Input file type (from content): ") + fh.typeToName(in_type), 2);
		}

		if (in_type==FileHandler::UNKNOWN)
		{
			writeLog_("Error: Could not determine input file type!");
			return PARSE_ERROR;
		}

	
		//output file names and types
		String out = getStringOption_("out");
		FileHandler::Type out_type = fh.nameToType(getStringOption_("out_type"));
			
		if (out_type==FileHandler::UNKNOWN)
		{
			out_type = fh.getTypeByFileName(out);
		}

		if (out_type==FileHandler::UNKNOWN)
		{
			writeLog_("Error: Could not determine output file type!");
			return PARSE_ERROR;
		}

		writeDebug_(String("Output file type: ") + fh.typeToName(out_type), 1);
			
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
		typedef MSExperiment< DPeak<1> > MSExperimentType;
		MSExperimentType exp;
		
		typedef MSExperimentType::SpectrumType SpectrumType;

		typedef DFeatureMap<2> FeatureMapType;

		writeDebug_(String("Loading input file"), 1);
			
		if (in_type == FileHandler::FEATURE)
		{
			// This works because DFeature<DIM> is derived from DPeak<DIM>.
			// However you will lose information and waste memory.
			// Enough reasons to issue a warning!
			writeLog_("Warning: Converting features to peaks. You will lose information!");	
			FeatureMapType fm;
			DFeatureMapFile().load(in,fm);
			fm.sortByPosition();
			exp.set2DData(fm);
		}
		else
		{
			fh.loadExperiment(in,exp,in_type);
		}
	
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		writeDebug_(String("Writing output file"), 1);
			
		if (out_type == FileHandler::MZDATA)
		{
			MzDataFile().store(out,exp);			
		}
		else if (out_type == FileHandler::MZXML)
		{
			MzXMLFile().store(out,exp);				
		}
		else if (out_type == FileHandler::DTA2D)
		{
			DTA2DFile().store(out,exp);			
		}
		else if (out_type == FileHandler::FEATURE)
		{
			// This works because DFeature<DIM> is derived from DPeak<DIM>.
			// However the feature specific information is only defaulted.
			// Enough reasons to issue a warning!
			writeLog_("Warning: Converting peaks into features.  This is only a hack - use at your own risk!");	
			FeatureMapType feature_map;
			static_cast<ExperimentalSettings>(feature_map) = exp;
			feature_map.reserve(exp.getSize());
			typedef FeatureMapType::FeatureType FeatureType;
			FeatureType feature;
			feature.setQuality(0,1); // override default
			feature.setQuality(1,1); // override default
			feature.setOverallQuality(1); // override default
			for ( MSExperimentType::ConstIterator spec_iter = exp.begin();
						spec_iter != exp.end();
						++spec_iter
					)
			{
				feature.setPos( DimensionDescription<LCMS_Tag>::RT,
												spec_iter->getRetentionTime() );
				for ( SpectrumType::ConstIterator peak1_iter = spec_iter->begin();
							peak1_iter != spec_iter->end();
							++peak1_iter
						)
				{
					feature.setPos( DimensionDescription<LCMS_Tag>::MZ,
												  peak1_iter->getPos() );
					feature.setIntensity(peak1_iter->getIntensity());
					feature_map.push_back(feature);
				}
			}
			feature_map.updateRanges();
			DFeatureMapFile().store(out,feature_map);
		}
		else
		{
			writeLog_("Unknown output file type given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;					
		}
			
		return EXECUTION_OK;
	}
};

int main( int argc, char ** argv )
{
	TOPPFileConverter tool;
	return tool.main(argc,argv);
}

/// @endcond
