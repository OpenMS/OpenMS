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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/MATH/MISC/BilinearInterpolation.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>


using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page Resampler Resampler
	
	@brief Resampler can be used to transform an LC/MS map into a resampled map or a pgm image.

	The input is first resampled into a matrix using bilinear interpolation.
	Then the content of the matrix is written into a mzData File or a pgm image
 	(pgm = portable network graphics, a very simple image file format).
	The output has a uniform spacing in both dimensions regardless of the input.
	You can output the data in transposed order, reverse video, with gamma correction, etc.

	@todo fix test (Clemens)
	@todo support for a better graphics format like png - use Qt (Clemens)
	@todo output IT range (Clemens)
	@todo maybe we could include support for one-dimensional resampling ("-cols auto -rows auto") for mzData output (Clemens)

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPResampler
	: public TOPPBase
{
 public:
	TOPPResampler()
	: TOPPBase("Resampler", "transform an LC/MS map into a resampled map or a pgm image")
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","input file in MzData format");

		// Note that we can have two output files.  At least one of them should be
		// specified.
		registerStringOption_("out","<file>","","output file in MzData format", false);
		registerStringOption_("pgm","<file>","","output file in plain PGM format", false);
		addText_("(Either -out or -pgm must be specified.)");

		addEmptyLine_();
		addText_("Parameters affecting the resampling:");
		registerStringOption_("mz","[min]:[max]",":","mass-to-charge range in input to be resampled", false);
		registerStringOption_("rt","[min]:[max]",":","retention time range in input to be resampled", false);
		registerIntOption_("cols_mz","<number>",101,"peaks per spectrum in output (image width)", false);
		registerIntOption_("rows_rt","<number>",101,"number of spectra in output (image height)", false);
		registerFlag_("transpose","flag to transpose the resampled matrix (RT vs. m/z)");

		addEmptyLine_();
		addText_("Parameters affecting the conversion from intensity to brightness:");
		registerIntOption_("maxval","<number>",255,"maximum brightness",false);
		registerDoubleOption_("scale","<factor>",0,"scaling factor for brightness",false);
		registerDoubleOption_("gamma","<value>",1.,"apply gamma correction",false);
		registerFlag_("reverse","flag to switch on reverse video");

		addEmptyLine_();
		addText_("In mzData output, peaks are ordered ascending in RT and m/z.");
		addText_("In pgm output, dimensions run bottom-up in RT and left-right in m/z.");
	}

	ExitCodes main_(int , char**)
	{
	
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		//file names

		String in = getStringOption_("in");
		inputFileReadable_(in);
					
		String out = getStringOption_("out");
		TOPPBase::ParameterInformation const & pi_out = findEntry_("out");
		bool has_out = (out != pi_out.default_value);

		String pgm = getStringOption_("pgm");
		TOPPBase::ParameterInformation const & pi_pgm = findEntry_("pgm");
		bool has_pgm = (pgm != pi_pgm.default_value);

		if ( !has_out && !has_pgm )
		{
			writeLog_("You need to specify an output destination using parameters \"out\" or \"pgm\".");
			return MISSING_PARAMETERS;
		}

		if ( has_out )
		{
			outputFileWritable_(out);
		}
		
		if ( has_pgm )
		{
			outputFileWritable_(pgm);
		}
		
		//parse RT and m/z range
		String rt = getStringOption_("rt");
		String mz = getStringOption_("mz");
		double rt_l, rt_u, mz_l, mz_u;
		//initialize ranges
		mz_l = rt_l = -1 * numeric_limits<double>::max();
		mz_u = rt_u = numeric_limits<double>::max();
				
		//rt
		parseRange_(rt,rt_l,rt_u);
		writeDebug_("rt lower/upper bound: " + String(rt_l) + " / " + String(rt_u),1);	
		//mz
		parseRange_(mz,mz_l,mz_u);
		writeDebug_("mz lower/upper bound: " + String(mz_l) + " / " + String(mz_u),1);	

		//load needed data
		typedef MSExperiment< DPeak<1> > MSExperimentType;
		typedef MSExperimentType::SpectrumType SpectrumType;
		MSExperimentType exp;
		MzDataFile f;
		f.getOptions().setRTRange(DRange<1>(rt_l,rt_u));
		f.getOptions().setMZRange(DRange<1>(mz_l,mz_u));
		f.load(in,exp);			
	
		std::stringstream comments;
		
		//basic info
		exp.updateRanges();

		//update RT and m/z to the real data if no boundary was given
		if (rt_l == -1 * numeric_limits<double>::max()) rt_l = exp.getMinRT();
		if (rt_u == numeric_limits<double>::max()) rt_u = exp.getMaxRT();
		if (mz_l == -1 * numeric_limits<double>::max()) mz_l = exp.getMinMZ();
		if (mz_u == numeric_limits<double>::max()) mz_u = exp.getMaxMZ();

		comments <<
			"generated by TOPP Resampler on " << Date::now() << "\n"
			"number of peaks: " << exp.getSize() << "\n"
 			"ranges before resampling:\n"
			"  RT: " << exp.getMinRT()  << ":" << exp.getMaxRT()  << "\n"
			"  MZ: " << exp.getMinMZ()  << ":" << exp.getMaxMZ()  << "\n"
			"  IT: " << exp.getMinInt() << ":" << exp.getMaxInt() << "\n"
			;

		int rows = getIntOption_("rows_rt");
		if ( rows < 1 )
		{
			writeLog_("Error: must have at least 1 row.");
			return ILLEGAL_PARAMETERS;
		}

		int cols = getIntOption_("cols_mz");
		if ( cols < 1 )
		{
			writeLog_("Error: must have at least 1 column.");
			return ILLEGAL_PARAMETERS;
		}


		BilinearInterpolation<double,double> bilip;
		bilip.getData().resize(rows,cols);

		bool transpose = getFlag_("transpose");
		if ( !transpose )
		{ // not transposed
			
			bilip.setMapping_0( 0, rt_u, rows-1, rt_l ); // scans run bottom-up
			bilip.setMapping_1( 0, mz_l, cols-1, mz_u ); // peaks run left-right

			for ( MSExperimentType::ConstIterator spec_iter = exp.begin();
						spec_iter != exp.end();
						++spec_iter
					)
			{
				double const rt = spec_iter->getRetentionTime();
				for ( SpectrumType::ConstIterator peak1_iter = spec_iter->begin();
							peak1_iter != spec_iter->end();
							++peak1_iter
						)
				{
					bilip.addValue(rt,peak1_iter->getPos(),peak1_iter->getIntensity());
				}
			}
		}
		else
		{ // transposed

			bilip.setMapping_0( 0, mz_u, rows-1, mz_l ); // spectra run bottom-up
			bilip.setMapping_1( 0, rt_l, cols-1, rt_u ); // scans run left-right

			for ( MSExperimentType::ConstIterator spec_iter = exp.begin();
						spec_iter != exp.end();
						++spec_iter
					)
			{
				double const rt = spec_iter->getRetentionTime();
				for ( SpectrumType::ConstIterator peak1_iter = spec_iter->begin();
							peak1_iter != spec_iter->end();
							++peak1_iter
						)
				{
					bilip.addValue(peak1_iter->getPos(),rt,peak1_iter->getIntensity());
				}
			}

		} // if
				
		int maxval = getIntOption_("maxval");
		double scale = getDoubleOption_("scale");
		double gamma = getDoubleOption_("gamma");
		bool reverse = getFlag_("reverse");


		if ( !pgm.empty() )
		{
			// all data in the matrix is directly written to file in pgm format

			std::ofstream pgm_file(pgm.c_str());
			bilip .getData()
				.writePGM ( pgm_file,
										maxval,
										scale,
										gamma,
										reverse,
										comments.str()
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
				SpectrumType & spectrum = exp_resampled[rows-row_index-1]; // reversed order so that retention times are increasing again

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
