// -*- mode: C++; tab-width: 2; -*-
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
#include <OpenMS/MATH/MISC/BilinearInterpolation.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include <QtGui/QImage>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_Resampler Resampler
	
	@brief Resampler can be used to transform an LC/MS map into a resampled map or a png image.
	
	When writing an mzData file, all spectra are resampled with a new sampling
	rate.  The number of spectra does not change.
	
	When writing an image, the input is first resampled into a matrix using
	bilinear forward resampling.  Then the content of the matrix is written to
	a PNG file.  The output has a uniform spacing in both dimensions regardless
	of the input.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_Resampler.cli
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
		registerOutputFile_("out", "<file>", "", "output file in mzData format or png format");
		registerFlag_("image", "Activates image mode (a png is written instead of a mzData file.");

		addEmptyLine_();
		addText_("Parameters affecting the MzData file:");
		registerDoubleOption_("sampling_rate", "<rate>", 0.1, "New sampling rate in m/z dimension", false);
		setMinFloat_("sampling_rate",0.0);
		
		addEmptyLine_();
		addText_("Parameters affecting the PNG file:");
		registerIntOption_("width", "<number>", 1000, "Number of pixels in m/z dimension.\nIf 0, for one pixel per Th.", false);
		setMinInt_("width",0);
		registerIntOption_("height", "<number>", 1000, "Number of pixels in RT dimension.\nIf 0, for one pixel per spectrum.", false);
		setMinInt_("height",0);
		registerStringOption_("gradient", "<gradient>", "", "Intensity gradient that defines colors for the range between 0 and 100.\n"
													"Example: '0,#FFFFFF;50,#FF0000;100,#000000'", false);
		registerDoubleOption_("maxintensity", "<int>", 0,
													"Maximum peak intensity used to determine range for colors.\n"
													"If 0, this is determined from data.", false);
		registerFlag_("log_intensity", "Apply logarithm to intensity values");
		registerFlag_("transpose", "flag to transpose the resampled matrix (RT vs. m/z).\n"
															 "Per default, dimensions run bottom-up in RT and left-right in m/z.");
	}

	ExitCodes main_(int , const char**)
	{
		//----------------------------------------------------------------
		// load data
		//----------------------------------------------------------------
		String in = getStringOption_("in");
		String out = getStringOption_("out");
		MSExperiment<> exp;
		MzDataFile f;
		f.setLogType(log_type_);
		f.load(in, exp);
		
		//----------------------------------------------------------------
		// PNG image
		//----------------------------------------------------------------
		if (getFlag_("image"))
		{
			exp.updateRanges(1);
	    
			Int rows = getIntOption_("height");
	    if (rows == 0)
	    {
	      rows = exp.size();
	    }
	    if ( rows <= 0 )
			{
				writeLog_("Error: Zero rows is not possible.");
				return ILLEGAL_PARAMETERS;
			}
	
			Int cols = getIntOption_("width");
	    if (cols == 0)
	    {
	      cols = UInt(ceil(exp.getMaxMZ() - exp.getMinMZ()));
	    }
			if ( cols <= 0 )
			{
				writeLog_("Error: Zero columns is not possible.");
				return ILLEGAL_PARAMETERS;
			}
			
			//----------------------------------------------------------------
			//Do the actual resampling
			BilinearInterpolation<DoubleReal, DoubleReal> bilip;
			bilip.getData().resize(rows, cols);
			if (!getFlag_("transpose"))
			{
				bilip.setMapping_0( 0, exp.getMaxRT(), rows-1, exp.getMinRT() ); // scans run bottom-up
				bilip.setMapping_1( 0, exp.getMinMZ(), cols-1, exp.getMaxMZ() ); // peaks run left-right
		
				for ( MSExperiment<>::Iterator spec_iter = exp.begin(); spec_iter != exp.end(); ++spec_iter)
				{
					if (spec_iter->getMSLevel()!=1) continue;
					for ( MSExperiment<>::SpectrumType::ConstIterator peak1_iter = spec_iter->begin(); peak1_iter != spec_iter->end(); ++peak1_iter )
					{
						bilip.addValue(spec_iter->getRT(), peak1_iter->getMZ(), peak1_iter->getIntensity());
					}
				}
			}
			else // transpose
			{ 
				bilip.setMapping_0( 0, exp.getMaxMZ(), rows-1, exp.getMinMZ() ); // spectra run bottom-up
				bilip.setMapping_1( 0, exp.getMinRT(), cols-1, exp.getMaxRT() ); // scans run left-right
	
				for ( MSExperiment<>::Iterator spec_iter = exp.begin(); spec_iter != exp.end(); ++spec_iter )
				{
					if (spec_iter->getMSLevel()!=1) continue;
					for ( MSExperiment<>::SpectrumType::ConstIterator peak1_iter = spec_iter->begin(); peak1_iter != spec_iter->end(); ++peak1_iter)
					{
						bilip.addValue(peak1_iter->getMZ(),  spec_iter->getRT(), peak1_iter->getIntensity());
					}
				}
			}
			
			//----------------------------------------------------------------
			//create and store image
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
			image.save(out.toQString(), "PNG");
		}
		//----------------------------------------------------------------
		// MzData file
		//----------------------------------------------------------------
		else
		{
			DoubleReal sampling_rate = getDoubleOption_("sampling_rate");
			
			LinearResampler lin_resampler;
			Param resampler_param;
			resampler_param.setValue("spacing",sampling_rate);
			lin_resampler.setParameters(resampler_param);
	
      // resample and filter every scan
      for (UInt i = 0; i < exp.size(); ++i)
      {
      	MSExperiment<>::SpectrumType resampled_spectrum;
        lin_resampler.raster(exp[i],resampled_spectrum);
        exp[i].swap(resampled_spectrum);
        exp[i].getMetaDataArrays().clear();
      }
			MzDataFile f;
			f.setLogType(log_type_);
			f.store(out, exp);
		}
		
		return EXECUTION_OK;
	}

};


int main( int argc, const char** argv )
{
	TOPPResampler tool;
	return tool.main(argc, argv);
}

/// @endcond
