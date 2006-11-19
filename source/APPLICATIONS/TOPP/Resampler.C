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
// $Maintainer: Clemens Gröpl $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/MATH/MISC/BilinearInterpolation.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>


using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page Resampler Resampler
	
	@brief doc missing
	
	longer doc missing
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPResampler
	: public TOPPBase
{
 public:
	TOPPResampler()
		: TOPPBase("Resampler")
	{
			
	}
	
 protected:
	void printToolUsage_() const
	{
		cerr << endl
				 << getToolName() << " -- transform a LC-MS map into a resampled pgm image." << endl
				 << "Version: " << VersionInfo::getVersion() << endl
				 << endl
				 << "Usage:" << endl
				 << "  " << getToolName() << " [options]" << endl
				 << endl <<
			"Options are:\n"
			"   -in <file>        input file\n"
			"   -in_type <type>   input file type (default: determined from input file extension)\n"
			"                    (Valid input types are: 'mzData', 'mzXML', 'DTA2D', 'ANDIMS' (cdf).)\n"
			"   -out <file>       output file (PGM format)\n"
			" Parameters affecting the resampling\n"
			"   -rt min:max       retention time range to be resampled for output\n"
			"   -rows <number>    number of rows in output\n"
			"   -mz min:max       mass-to-charge range to be resampled for output\n"
			"   -cols <number>    number of columns in output\n"
			"\n"
			" Parameters affecting the conversion from intensity to brightness:\n"
			"   -maxval <number>  maximum brightness\n"
			"   -scale <number>   scaling factor for brightness\n"
			"   -reverse          flag to switch on reverse video\n"
				 << endl;
	}
	
	void printToolHelpOpt_() const
	{
		cerr << endl
				 << getToolName() << endl
				 << endl
				 << "INI options:" << endl
				 << "" << endl
				 << "... to be documented ..." << endl
				 << "" << endl
				 << "  in        input file name" << endl
				 << "  in_type   input file type (default: determined from input file name extension)" << endl
				 << endl
				 << "INI File example section:" << endl
				 << "  <ITEM name=\"in\" value=\"example.mzData\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"in_type\" value=\"MZDATA\" type=\"string\"/>" << endl
				 << "  ... and so on ..." << endl;
	}
	
	void setOptionsAndFlags_()
	{
		options_["-in"] = "in";
		options_["-in_type"] = "in_type";
		options_["-out"] = "out";
		options_["-rows"] = "rows";
		options_["-cols"] = "cols";
		options_["-mz"] = "mz";
		options_["-rt"] = "rt";
		options_["-maxval"] = "maxval";
		options_["-scale"] = "scale";
		flags_["-reverse"] = "reverse";
	}
	
	ExitCodes main_(int , char**)
	{
	
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		//file names
		String in = getParamAsString_("in");
		writeDebug_(String("Input file: ") + in, 1);
			
		//file names
		String out = getParamAsString_("out");
		writeDebug_(String("Output file: ") + out, 1);
			
		//file type
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

		cout << endl
				 << "file name: " << in << endl
				 << "file type: " <<  fh.typeToName(in_type) << endl
				 << endl;
			
		typedef MSExperiment< DPeak<1> > MSExperimentType;

		MSExperimentType exp;
	
		std::stringstream comments;
		
		if (! fh.loadExperiment(in,exp,in_type) )
		{
			writeLog_("Unsupported or corrupt input file. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;			
		}
		
		//basic info
		exp.updateRanges();

		comments <<
			"number of peaks: " << exp.getSize() << "\n"
			"ranges before resampling:\n"
			"  RT: " << exp.getMinRT()  << ":" << exp.getMaxRT()  << "\n"
			"  MZ: " << exp.getMinMZ()  << ":" << exp.getMaxMZ()  << "\n"
			"  IT: " << exp.getMinInt() << ":" << exp.getMaxInt() << "\n"
			;

		BilinearInterpolation<double,double> bilip;

		int rows = getParamAsInt_("rows",100);
		int cols = getParamAsInt_("cols",100);
		bilip.getData().resize(rows,cols);

		String rt = getParamAsString_("rt",":");
		String mz = getParamAsString_("mz",":");
		String tmp;
		double rt_l, rt_u, mz_l, mz_u;

		//convert bounds to numbers
		try
		{
			//rt
			tmp = rt.prefix(':');
			if (tmp!="")
			{
				rt_l = tmp.toDouble();
			}
			else
			{
				rt_l = exp.getMinRT();
			}
			tmp = rt.suffix(':');
			if (tmp!="")
			{
				rt_u = tmp.toDouble();
			}
			else
			{
				rt_u = exp.getMaxRT();
			}
			writeDebug_("rt lower:upper bound: " + String(rt_l) + " : " + String(rt_u),1);	
				
			//mz
			tmp = mz.prefix(':');
			if (tmp!="")
			{
				mz_l = tmp.toDouble();
			}
 			else
			{
				mz_l = exp.getMinMZ();
			}
			tmp = mz.suffix(':');
			if (tmp!="")
			{
				mz_u = tmp.toDouble();
			}
 			else
			{
				mz_u = exp.getMaxMZ();
			}
			writeDebug_("mz lower:upper bound: " + String(mz_l) + " : " + String(mz_u),1);	
				
		}
		catch(Exception::ConversionError& e)
		{
			writeLog_(String("Invalid boundary '") + tmp + "' given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;			
		}

		bilip.setMapping_0( 0, rt_l, rows-1, rt_u );
		bilip.setMapping_1( 0, mz_l, cols-1, mz_u );

		comments <<
			"ranges after resampling:\n"
			"  RT: " << rt_l << ":" << rt_u << "\n"
			"  MZ: " << mz_l << ":" << mz_u << "\n"
			;

		enum DimensionId
			{
				RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
				MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
			};
		for ( MSExperimentType::PeakIterator<DPeak<1> > iter(exp.peakBegin()); iter != exp.peakEnd(); ++iter )
		{
			bilip.addValue(iter.getRt(),iter->getPos(),iter->getIntensity());
		}
				
		int maxval = getParamAsInt_("maxval",255);
		double scale = getParamAsDouble_("scale",0);
		bool reverse = getParamAsBool_("reverse",false);

		std::ofstream out_file(out.c_str());
		bilip.getData().writePGM(out_file,maxval,scale,reverse,"generated by TOPP Resampler on "+Date::now()+'\n'+comments.str() );
		out_file.close();
			
	return EXECUTION_OK;
	}	

};


int main( int argc, char ** argv )
{
	TOPPResampler tool;
	return tool.main(argc,argv);
}

/// @endcond


// TODO
/*
	- output IT range
	- apply gamma transformation (optional)
	- copy points to a MSExperiment
	- write resampled image out as mzData etc.
*/
