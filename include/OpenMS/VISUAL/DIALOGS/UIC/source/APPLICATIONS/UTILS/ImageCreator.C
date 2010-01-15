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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
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
	@page UTILS_ImageCreator ImageCreator

	@brief Transforms an LC/MS map into a png image.

	The input is first resampled into a matrix using
	bilinear forward resampling.  Then the content of the matrix is written to
	a PNG file.  The output has a uniform spacing in both dimensions regardless
	of the input.

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_ImageCreator.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPImageCreator
	: public TOPPBase
{
 public:
	TOPPImageCreator()
		: TOPPBase("ImageCreator",
							 "Transforms an LC/MS map into a PNG image.",false)
	{
	}

 protected:

	void addMS2Point_(int x, int y, QImage& image, QColor color = Qt::black,
									 Size size = 2)
		{		
			int h = image.height(), w = image.width();
			vector<int> xs(1, x), ys(1, y);
			if (size == 2)
			{
				int xtemp[] = {x - 1, x, x, x + 1};
				int ytemp[] = {y, y - 1, y + 1, y};
				xs = vector<int>(xtemp, xtemp + 4);
				ys = vector<int>(ytemp, ytemp + 4);
			}
			else if (size == 3)
			{
				int xtemp[] = {x - 2, x - 1, x - 1, x, x, x + 1, x + 1, x + 2};
				int ytemp[] = {y, y + 1, y - 1, y + 2, y - 2, y + 1, y - 1, y};
				xs = vector<int>(xtemp, xtemp + 8);
				ys = vector<int>(ytemp, ytemp + 8);
			}			
			for (Size i = 0; i < xs.size(); ++i)
			{
				int xi = xs[i], yi = ys[i];
				if ((xi > 0) && (xi < w) && (yi > 0) && (yi < h))
				{
					image.setPixel(xi, yi, color.rgb());
				}
			}
		}

	
	void markMS2Locations_(MSExperiment<>& exp, QImage& image, bool transpose,
												QColor color, Size size)
		{
			double xcoef = image.width(), ycoef = image.height();
			if (transpose)
			{
				xcoef /= exp.getMaxRT() - exp.getMinRT();
				ycoef /= exp.getMaxMZ() - exp.getMinMZ();
			}
			else
			{
				xcoef /= exp.getMaxMZ() - exp.getMinMZ();
				ycoef /= exp.getMaxRT() - exp.getMinRT();
			}
			for (MSExperiment<>::Iterator spec_iter = exp.begin();
					 spec_iter != exp.end(); ++spec_iter)
			{
				if (spec_iter->getMSLevel() == 2)
				{
					double mz = spec_iter->getPrecursors()[0].getMZ();
					double rt = exp.getPrecursorSpectrum(spec_iter)->getRT();
					int x, y;
					if (transpose)
					{
						x = int(xcoef * (rt - exp.getMinRT()));
						y = int(ycoef * (exp.getMaxMZ() - mz));
					}
					else
					{
						x = int(xcoef * (mz - exp.getMinMZ()));
						y = int(ycoef * (exp.getMaxRT() - rt));
					}
					addMS2Point_(x, y, image, color, size);
				}
			}
		}
	

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "input file ");
		setValidFormats_("in", StringList::create("mzML"));
		registerOutputFile_("out", "<file>", "",
												"output file in PNG format");
		setValidFormats_("out", StringList::create("PNG"));

		registerIntOption_("width", "<number>", 1024, "Number of pixels in m/z dimension.\nIf 0, one pixel per Th.", false);
		setMinInt_("width",0);
		registerIntOption_("height", "<number>", 1024, "Number of pixels in RT dimension.\nIf 0, one pixel per spectrum.", false);
		setMinInt_("height",0);
		registerStringOption_("gradient", "<gradient>", "", "Intensity gradient that defines colors for the range between 0 and 100.\n"
													"Example: '0,#FFFFFF;50,#FF0000;100,#000000'", false);
		registerDoubleOption_("maxintensity", "<int>", 0, "Maximum peak intensity used to determine range for colors.\n"
													"If 0, this is determined from the data.", false);
		registerFlag_("log_intensity", "Apply logarithm to intensity values");
		registerFlag_("transpose", "flag to transpose the resampled matrix (RT vs. m/z).\n"
															 "Per default, dimensions run bottom-up in RT and left-right in m/z.");
		registerFlag_("precursors", "Mark locations of MS2 precursors.\n"
									"Implied if 'precursor_color' or 'precursor_size' are set.");
		registerStringOption_("precursor_color", "<color>", "#000000", "Color for precursor marks (color code or word, e.g. 'black')", false);
		registerIntOption_("precursor_size", "<number>", 2,
											 "Size of the precursor marks", false);
		setMinInt_("precursor_size", 1);
		setMaxInt_("precursor_size", 3);
	}

	ExitCodes main_(int , const char**)
	{
		//----------------------------------------------------------------
		// load data
		//----------------------------------------------------------------
		String in = getStringOption_("in");
		String out = getStringOption_("out");
		MSExperiment<> exp;
		MzMLFile f;
		f.setLogType(log_type_);
		f.load(in, exp);

		exp.updateRanges(1);

		SignedSize rows = getIntOption_("height");
    if (rows == 0)
    {
      rows = exp.size();
    }
    if ( rows <= 0 )
		{
			writeLog_("Error: Zero rows is not possible.");
			return ILLEGAL_PARAMETERS;
		}

    SignedSize cols = getIntOption_("width");
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
			// scans run bottom-up:
			bilip.setMapping_0(0, exp.getMaxRT(), rows-1, exp.getMinRT());
			// peaks run left-right:
			bilip.setMapping_1(0, exp.getMinMZ(), cols-1, exp.getMaxMZ());

			for (MSExperiment<>::Iterator spec_iter = exp.begin();
					 spec_iter != exp.end(); ++spec_iter)
			{
				if (spec_iter->getMSLevel() != 1) continue;
				for (MSExperiment<>::SpectrumType::ConstIterator peak1_iter =
							 spec_iter->begin(); peak1_iter != spec_iter->end();
						 ++peak1_iter)
				{
					bilip.addValue(spec_iter->getRT(), peak1_iter->getMZ(),
												 peak1_iter->getIntensity());
				}
			}
		}
		else // transpose
		{
			// spectra run bottom-up:
			bilip.setMapping_0(0, exp.getMaxMZ(), rows-1, exp.getMinMZ());
			// scans run left-right:
			bilip.setMapping_1(0, exp.getMinRT(), cols-1, exp.getMaxRT());

			for (MSExperiment<>::Iterator spec_iter = exp.begin();
						spec_iter != exp.end(); ++spec_iter)
			{
				if (spec_iter->getMSLevel()!=1) continue;
				for (MSExperiment<>::SpectrumType::ConstIterator peak1_iter =
							 spec_iter->begin(); peak1_iter != spec_iter->end();
						 ++peak1_iter)
				{
					bilip.addValue(peak1_iter->getMZ(), spec_iter->getRT(),
												 peak1_iter->getIntensity());
				}
			}
		}

		//----------------------------------------------------------------
		//create and store image
		int scans = (int) bilip.getData().sizePair().first;
		int peaks = (int) bilip.getData().sizePair().second;

		MultiGradient gradient;
		String gradient_str = getStringOption_("gradient");
		if (gradient_str!="")
		{
			gradient.fromString(String("Linear|") + gradient_str);
		}
		else
		{
			gradient.fromString("Linear|0,#FFFFFF;2,#FFFF00;11,#FFAA00;32,#FF0000;55,#AA00FF;78,#5500FF;100,#000000");
		}

    bool use_log = getFlag_("log_intensity");
    writeDebug_("log_intensity: " + String(use_log), 1);

    QImage image(peaks, scans, QImage::Format_RGB32);
		DoubleReal factor = getDoubleOption_("maxintensity");
		if ( factor == 0 )
		{
			factor = (*std::max_element(bilip.getData().begin(),
																	bilip.getData().end()));
		}
    // logarithmize max. intensity as well:
    if (use_log) factor = std::log(factor);
		
		factor /= 100.0;
		for (int i = 0; i < scans; ++i)
		{
			for (int j = 0; j < peaks; ++j)
			{
				double value = bilip.getData().getValue(i, j);
				if (use_log) value = std::log(value);
				image.setPixel(j, i, gradient.interpolatedColorAt(value/factor).
											 rgb());
			}
		}

		if (getFlag_("precursors") || setByUser_("precursor_color") ||
				setByUser_("precursor_size"))
		{
			markMS2Locations_(exp, image, getFlag_("transpose"),
											  getStringOption_("precursor_color").toQString(),
											  Size(getIntOption_("precursor_size")));
		}
		
		image.save(out.toQString(), "PNG");
		return EXECUTION_OK;
	}

};


int main( int argc, const char** argv )
{
	TOPPImageCreator tool;
	return tool.main(argc, argv);
}

/// @endcond
