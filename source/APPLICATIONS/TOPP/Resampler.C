// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/MATH/MISC/BilinearInterpolation.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/MultiGradient.h>

#include <QtGui/QImage>
#include <QtGui/QColor>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page Resampler Resampler

 @brief Resampler can be used to transform an LC/MS map into a resampled map or a png image.

 The input is first resampled into a matrix using bilinear interpolation.
 Then the content of the matrix is written into a mzData File or a png image.
 The output has a uniform spacing in both dimensions regardless of the input.

 @improvement maybe we could include support for one-dimensional resampling e.g. "-cols auto -rows auto" for mzData output (Clemens)
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPResampler
	: public TOPPBase
{
 public:
	TOPPResampler()
		: TOPPBase("Resampler", "Transforms an LC/MS map into a resampled map or a png image.")
	{
	}

 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "input file ");
		setValidFormats_("in",StringList::create("mzData"));

		// Note that we can have two output files.  At least one of them should be specified.
		registerOutputFile_("out", "<file>", "", "output file ", false);
		setValidFormats_("out",StringList::create("mzData"));
		registerOutputFile_("png", "<file>", "", "output file in PNG format", false);
		addText_("(Either -out or -png must be specified.)");

		addEmptyLine_();
		addText_("Parameters affecting the resampling:");
		registerStringOption_("mz", "[min]:[max]", ":", "mass-to-charge range in input to be resampled", false);
		registerStringOption_("rt", "[min]:[max]", ":", "retention time range in input to be resampled", false);
		registerIntOption_("cols_mz", "<number>", 101, "peaks per spectrum in output (image width); use 0 for one col per Th", false);
		registerIntOption_("rows_rt", "<number>", 101, "number of spectra in output (image height); use 0 for one row per scan", false);
		registerFlag_("transpose", "flag to transpose the resampled matrix (RT vs. m/z)");

		addEmptyLine_();
		addText_("Parameters affecting the image:");
		registerStringOption_("gradient", "<gradient>", "", "Intensity gradient that defines colors for the range "
													"between 0 and 100. Example: '0,#FFFFFF;50,#FF0000;100,#000000'", false);
		registerDoubleOption_("maxintensity", "<maxintensity>", 0,
													"Maximum peak intensity used to determine range for colors.  "
													"If 0, this is determined from data.", false);
    registerFlag_("log_intensity", "apply logarithm to intensity values");
		addEmptyLine_();
		addText_("In mzData output, peaks are ordered ascending in RT and m/z.");
		addText_("In png output, dimensions run bottom-up in RT and left-right in m/z.");
	}

	ExitCodes main_(int , const char**)
	{

		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------

		String in = getStringOption_("in");

		bool output_defined=false;
		String out = getStringOption_("out");
		if (out!="")
		{
			output_defined = true;
		}
		String png = getStringOption_("png");
		if (png!="")
		{
			output_defined = true;
		}

		if (!output_defined)
		{
			writeLog_("You need to specify an output destination using parameters \"out\" or \"png\".");
			return MISSING_PARAMETERS;
		}

		//parse RT and m/z range
		String rt = getStringOption_("rt");
		String mz = getStringOption_("mz");
		double rt_l, rt_u, mz_l, mz_u;
		//initialize ranges
		mz_l = rt_l = -1 * numeric_limits<double>::max();
		mz_u = rt_u = numeric_limits<double>::max();

		//rt
		parseRange_(rt, rt_l, rt_u);
		writeDebug_("rt lower/upper bound: " + String(rt_l) + " / " + String(rt_u), 1);
		//mz
		parseRange_(mz, mz_l, mz_u);
		writeDebug_("mz lower/upper bound: " + String(mz_l) + " / " + String(mz_u), 1);

		//load needed data
		typedef MSExperiment<Peak1D> MSExperimentType;
		typedef MSExperimentType::SpectrumType SpectrumType;
		MSExperimentType exp;
		MzDataFile f;
		f.setLogType(log_type_);
		f.getOptions().setRTRange(DRange<1>(rt_l, rt_u));
		f.getOptions().setMZRange(DRange<1>(mz_l, mz_u));
		f.load(in, exp);

		//basic info
		exp.updateRanges(1);

		//update RT and m/z to the real data if no boundary was given
		if (rt_l == -1 * numeric_limits<double>::max()) rt_l = exp.getMinRT();
		if (rt_u == numeric_limits<double>::max()) rt_u = exp.getMaxRT();
		if (mz_l == -1 * numeric_limits<double>::max()) mz_l = exp.getMinMZ();
		if (mz_u == numeric_limits<double>::max()) mz_u = exp.getMaxMZ();

		int rows = getIntOption_("rows_rt");
    // one row for each scan
    if (rows == 0)
    {
      rows = exp.size();
      writeDebug_("row count: "+ String(rows) + " [" + String(rt_l) + " - " + String(rt_u) + "]", 1);  		
    }
    if ( rows < 1 )
		{
			writeLog_("Error: must have at least 1 row.");
			return ILLEGAL_PARAMETERS;
		}

		int cols = getIntOption_("cols_mz");
    // one row for each Thomson
    if (cols == 0)
    {
      cols = int(mz_u - mz_l); 
      writeDebug_("column count: "+ String(cols) + " [" + String(mz_l) + " - " + String(mz_u) + "]", 1);     
    }
		if ( cols < 1 )
		{
			writeLog_("Error: must have at least 1 column.");
			return ILLEGAL_PARAMETERS;
		}


		BilinearInterpolation<double, double> bilip;
		bilip.getData().resize(rows, cols);

		bool transpose = getFlag_("transpose");
		if ( !transpose )
		{ // not transposed

			bilip.setMapping_0( 0, rt_u, rows-1, rt_l ); // scans run bottom-up
			bilip.setMapping_1( 0, mz_l, cols-1, mz_u ); // peaks run left-right

			for ( MSExperimentType::Iterator spec_iter = exp.begin();
						spec_iter != exp.end();
						++spec_iter
					)
			{
				if (spec_iter->getMSLevel()!=1)
				{
					continue;
				}
				double const rt = spec_iter->getRT();
				for ( SpectrumType::ConstIterator peak1_iter = spec_iter->begin();
							peak1_iter != spec_iter->end();
							++peak1_iter )
				{
					bilip.addValue(rt, peak1_iter->getMZ(), peak1_iter->getIntensity());
				}
			}
		}
		else
		{ // transposed

			bilip.setMapping_0( 0, mz_u, rows-1, mz_l ); // spectra run bottom-up
			bilip.setMapping_1( 0, rt_l, cols-1, rt_u ); // scans run left-right

			for ( MSExperimentType::Iterator spec_iter = exp.begin();
						spec_iter != exp.end();
						++spec_iter
					)
			{
				if (spec_iter->getMSLevel()!=1)
				{
					continue;
				}
				double const rt = spec_iter->getRT();
				for ( SpectrumType::ConstIterator peak1_iter = spec_iter->begin();
							peak1_iter != spec_iter->end();
							++peak1_iter
						)
				{
					bilip.addValue(peak1_iter->getMZ(), rt, peak1_iter->getIntensity());
				}
			}

		}
  
		if(!png.empty())
		{
			UInt scans = bilip.getData().sizePair().first;
			UInt peaks = bilip.getData().sizePair().second;

			MultiGradient gradient;
			String gradient_str = getStringOption_("gradient");
			if (gradient_str!="")
			{
				gradient.fromString(String("Linear|") + gradient_str);
			}
			else
			{
				gradient.fromString("Linear|0,#FFFFFF;2,#FFFF00;11,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000");
			}
      
      bool use_log = getFlag_("log_intensity");
      writeDebug_("log_intensity: " + String(use_log), 1);
			
      QImage image(peaks, scans, QImage::Format_RGB32);
			DoubleReal factor = getDoubleOption_("maxintensity");
			if ( factor == 0 )
			{
				factor = (*std::max_element(bilip.getData().begin(), bilip.getData().end()));
			}
      // logarithmize maxintensity as well
      if (use_log)
      {
        factor = std::log(factor); 
      }
			factor /= 100.0;
      // apply logarithm to intensities
      if (use_log)
      {
  			for (UInt i=0; i<scans; ++i)
  			{
  				for (UInt j=0; j<peaks; ++j)
  				{
  					image.setPixel(j, i, gradient.interpolatedColorAt(std::log(bilip.getData().getValue(i, j))/factor).rgb());
  				}
  			}
      }
      else
      {
        for (UInt i=0; i<scans; ++i)
        {
          for (UInt j=0; j<peaks; ++j)
          {
            image.setPixel(j, i, gradient.interpolatedColorAt(bilip.getData().getValue(i, j)/factor).rgb());
          }
        }
      }
			image.save(png.c_str(), "PNG");
		}

		if ( !out.empty() )
		{
			// all data in the matrix is copied to an MSExperiment,
			// which is then written to an mzData file.
			MSExperiment< Peak1D > exp_resampled;
			exp_resampled.resize(rows);
			for ( int row_index = 0; row_index < rows; ++row_index )
			{
				SpectrumType & spectrum = exp_resampled[rows-row_index-1]; // reversed order so that retention times are increasing again

				spectrum.setRT( bilip.index2key_0( row_index ) );
				spectrum.setMSLevel(1);
				spectrum.resize(cols);

				for ( int col_index = 0; col_index < cols; ++col_index )
				{
					typedef SpectrumType::PeakType PeakType;
					PeakType & peak = spectrum[col_index];

					peak.setIntensity( bilip.getData()(row_index, col_index) );
					peak.setMZ( bilip.index2key_1( col_index ) );

				} // col_index

			} // row_index

			MzDataFile f;
			f.setLogType(log_type_);
			f.store(out, exp_resampled);

		} // !out.empty()

		return EXECUTION_OK;
	}

};


int main( int argc, const char** argv )
{
	TOPPResampler tool;
	return tool.main(argc, argv);
}

/// @endcond
