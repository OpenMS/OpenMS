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

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>


#include "TOPPBase.h"

using namespace OpenMS;
using namespace std;


/**
	@page FileConverter FileConverter
	
	@brief Converts between different MS file formats.
	
	Supported input file types are: 'mzData', 'mzXML', 'DTA2D', 'cdf' (ANDI/MS).
	'feat' (features) is also supported but will lose feature specific information.
	
	Supported output file types are: 'mzData', 'mzXML', 'DTA2D'
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu -> @cond
/// @cond 

class TOPPFileConverter
	: public TOPPBase
{
	public:
		TOPPFileConverter()
			: TOPPBase("FileConverter")
		{
			
		}
	
	protected:
		void printToolUsage_()
		{
			cerr << endl
		       << tool_name_ << " -- converts between different MS file formats." << endl
		       << endl
		       << "Usage:" << endl
					 << "  " << tool_name_ << " [options]" << endl
					 << endl
					 << "Options are:" << endl
					 << "  -in <file>        input file" << endl
					 << "  -out <file>       output file" << endl
					 << "  -in_type <type>   input file type (default: determined from input file extension)" << endl
					 << "  -out_type <type>  output file type (default: determined from output file extension)" << endl
					 << endl
					 << "Valid input types are: 'mzData', 'mzXML', 'DTA2D', 'ANDIMS'" << endl
					 << "                       'feat' (features) can be converted, but will lose feature specific information" << endl
					 << "Valid output types are: 'mzData', 'mzXML', 'DTA2D'" << endl;	
		}
	
		void printToolHelpOpt_()
		{
			cerr << endl
		       << tool_name_ << endl
		       << endl
		       << "INI options:" << endl
					 << "  in         input file" << endl
					 << "  out        output file" << endl
					 << "  in_type    input file type (default: determined from input file extension)" << endl
					 << "  out_type   output file type (default: determined from output file extension)" << endl
					 << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"in_type\" value=\"MZDATA\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"out\" value=\"output.mzXML\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"out_type\" value=\"MZXML\" type=\"string\"/>" << endl;
		}
	
		void setOptionsAndFlags_()
		{
			options_["-out"] = "out";
			options_["-in"] = "in";
			options_["-out_type"] = "out_type";
			options_["-in_type"] = "in_type";
		}
	
		ExitCodes main_(int , char**)
		{
	
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input file names
			String in = getParamAsString_("in");
			writeDebug_(String("Input file: ") + in, 1);
			
			//input file type
			FileHandler fh;
			FileHandler::Type in_type = fh.nameToType(getParamAsString_("in_type",""));
			
			writeDebug_(String("Input file type (from command line): ") + fh.typeToName(in_type), 1);
			
			if (in_type==FileHandler::UNKNOWN)
			{
				in_type = fh.getTypeByFileName(in);
				writeDebug_(String("Input file type (from file extention): ") + fh.typeToName(in_type), 1);
			}	

			if (in_type==FileHandler::UNKNOWN)
			{
				in_type = fh.getTypeByContent(in);
				writeDebug_(String("Input file type (from file content): ") + fh.typeToName(in_type), 1);
			}
	
			//output file names and types
			String out = getParamAsString_("out");
			FileHandler::Type out_type = fh.nameToType(getParamAsString_("out_type",""));
			
			if (out_type==FileHandler::UNKNOWN)
			{
				out_type = fh.getTypeByFileName(out);
			}
			
			writeDebug_(String("Output file: ") + out, 1);
			writeDebug_(String("Output file type: ") + fh.typeToName(out_type), 1);
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			MSExperiment< DPeak<1> > exp;
			
			writeDebug_(String("Loading input file"), 1);
			
			if (in_type == FileHandler::FEATURE)
			{
				// This works because DFeature<DIM> is derived from DPeak<DIM>.
				// However you will lose information and waste memory.
				// Enough reasons to issue a warning!
				writeLog_("Warning:  Converting features to peaks.  You will lose information!");	
				DFeatureMap<2> fm;
				DFeatureMapFile().load(in,fm);
				fm.sortByPosition();
				exp.set2DData(fm);
			}
			else if (in_type != FileHandler::UNKNOWN)
			{
				fh.loadExperiment(in,exp,in_type);
			}
			else
			{
				writeLog_("Unknown input file type given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;			
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
			else
			{
				writeLog_("Unknown output file type given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;					
			}
			
			return OK;
		}
};

///@endcond

int main( int argc, char ** argv )
{
	TOPPFileConverter tool;
	return tool.main(argc,argv);
}
