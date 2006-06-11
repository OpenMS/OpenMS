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
// $Id: FileConverter.C,v 1.18 2006/06/09 14:46:55 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#ifdef ANDIMS_DEF
#include <OpenMS/FORMAT/ANDIFile.h>		
#endif

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
					 << "Valid input types are: 'mzData', 'mzXML', 'DTA2D', 'cdf' (ANDI/MS)" << endl
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
	
			//input file names and types
			String in = getParamAsString_("in");
			String in_type = getParamAsString_("in_type","");
		
			if (in_type=="")
			{
				in_type = in.suffix('.');
			}	
			in_type.toUpper();
			
			writeDebug_(String("Input file: ") + in, 1);
			writeDebug_(String("Input file type: ") + in_type, 1);
	
			//output file names and types
			String out = getParamAsString_("out");
			String out_type = getParamAsString_("out_type","");
		
			if (out_type=="")
			{
				out_type = out.suffix('.');
			}	
			out_type.toUpper();

			writeDebug_(String("Output file: ") + out, 1);
			writeDebug_(String("Output file type: ") + out_type, 1);
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
		
			//load input file data
			MSExperiment< DPeak<1> > exp;
			if (in_type == "MZDATA")
			{
				MzDataFile f;
				f.load(in,exp);			
			}
			else if (in_type == "MZXML")
			{
				MzXMLFile f;
				f.load(in,exp);				
			}
			else if (in_type == "CDF")
			{
	#ifdef ANDIMS_DEF
				ANDIFile f;
				f.load(in,exp);			
	#else
				writeLog_( String(" Unsupported file type '") + in_type + "' given. Aborting!");
				return INPUT_FILE_NOT_READABLE;			
	#endif	
			}
			else if (in_type == "DTA2D")
			{
				DTA2DFile f;
				f.load(in,exp);			
			}
			else if (in_type == "FEAT")
			{
				// This works because DFeature<DIM> is derived from DPeak<DIM>.
				// However you will lose information and waste memory.
				// Enough reasons to issue a warning!
				writeLog_("Warning:  Converting features to peaks.  You will lose information!");	
				DFeatureMapFile f;
				DFeatureMap<2> fm;
				f.load(in,fm);
				fm.sortByPosition();
				exp.set2DData(fm);
			}
			else
			{
				writeLog_( String("Unknown input file type '") + in_type + "' given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;			
			}
	
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
		
			if (out_type == "MZDATA")
			{
				MzDataFile f;
				f.store(out,exp);			
			}
			else if (out_type == "MZXML")
			{
				MzXMLFile f;
				f.store(out,exp);				
			}
			else if (out_type == "DTA2D")
			{
				DTA2DFile f;
				f.store(out,exp);			
			}
			else
			{
				writeLog_( String("Unknown output file type '") + out_type + "' given. Aborting!");
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
