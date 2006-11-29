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
#include <OpenMS/FORMAT/MzDataFile.h>
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
			"   -out <file>       output file (mzData format)\n"
			"   -pgm <file>       output file (plain PGM image format)\n"
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
			"   -transpose        flag to transpose the resampled matrix (RT vs. m/z)\n"
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
				 << endl
				 << "INI File example section:" << endl
				 << "  <ITEM name=\"in\" value=\"example.mzData\" type=\"string\"/>" << endl
				 << "  ... and so on ..." << endl;
	}
	
	void setOptionsAndFlags_()
	{
		options_["-in"] = "in";
		options_["-out"] = "out";
		options_["-pgm"] = "pgm";
		options_["-rows"] = "rows";
		options_["-cols"] = "cols";
		options_["-mz"] = "mz";
		options_["-rt"] = "rt";
		options_["-maxval"] = "maxval";
		options_["-scale"] = "scale";
		flags_["-reverse"] = "reverse";
		flags_["-transpose"] = "transpose";
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
		writeDebug_(String("Output file (mzData format): ") + out, 1);

		//file names
		String pgm = getParamAsString_("pgm");
		writeDebug_(String("Output file (plain PGM image format): ") + pgm, 1);
			
		typedef MSExperiment< DPeak<1> > MSExperimentType;

		MSExperimentType exp;
		MzDataFile().load(in,exp);			
	
		std::stringstream comments;
		
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

		comments <<
			"ranges after resampling:\n"
			"  RT: " << rt_l << ":" << rt_u << "\n"
			"  MZ: " << mz_l << ":" << mz_u << "\n"
			;

		bool transpose = getParamAsBool_("transpose",false);
		if ( !transpose )
		{ // normal ranges, no transposition

			bilip.setMapping_0( 0, rt_l, rows-1, rt_u );
			bilip.setMapping_1( 0, mz_l, cols-1, mz_u );

			for ( MSExperimentType::PIterator iter = exp.peakBegin(); iter != exp.peakEnd(); ++iter )
			{
				bilip.addValue(iter.getRt(),iter->getPos(),iter->getIntensity());
			}

		}
		else
		{ // flipped ranges, transposed matrix

			bilip.setMapping_0( 0, mz_l, cols-1, mz_u );
			bilip.setMapping_1( 0, rt_l, rows-1, rt_u );

			for ( MSExperimentType::PIterator iter = exp.peakBegin(); iter != exp.peakEnd(); ++iter )
			{
				bilip.addValue(iter->getPos(),iter.getRt(),iter->getIntensity());
			}

		} // if transpose
				
		int maxval = getParamAsInt_("maxval",255);
		double scale = getParamAsDouble_("scale",0);
		bool reverse = getParamAsBool_("reverse",false);


		if ( !pgm.empty() )
		{
			// all data in the matrix is directly written to file in pgm format

			std::ofstream pgm_file(pgm.c_str());
			bilip .getData()
				.writePGM ( pgm_file,
										maxval,
										scale,
										reverse,
										"generated by TOPP Resampler on "+Date::now()+'\n'+comments.str()
									);
			pgm_file.close();
		}

		if ( !out.empty() )
		{
			// all data in the matrix is copied to an MSExperiment,
			// which is then written to an mzData file.

			MSExperimentType exp_resampled;
			exp_resampled.resize(rows);

			for ( int row_index = 0; row_index < rows; ++row_index )
			{
				typedef MSExperimentType::SpectrumType SpectrumType;
				SpectrumType & spectrum = exp_resampled[row_index];

				spectrum.setRetentionTime( bilip.index2key_0( row_index ) );
				spectrum.setMSLevel(1);
				spectrum.resize(cols);

				for ( int col_index = 0; col_index < cols; ++col_index )
				{
					typedef SpectrumType::PeakType PeakType;
					PeakType & peak = spectrum[col_index];

					peak.setIntensity( bilip.getData()(row_index,col_index) );
					peak.setPos( bilip.index2key_1( col_index ) );

				} // col_index

			} // row_index

			MzDataFile().store(out,exp_resampled);			

		} // !out.empty()
		
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
