// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
	@page TOPP_Resampler Resampler

	@brief Resampler can be used to transform an LC/MS map into a resampled map.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ Resampler \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterSGolay  </td>
		</tr>
	</table>
</CENTER>

	When writing an peak file, all spectra are resampled with a new sampling
	rate. The number of spectra does not change.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_Resampler.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_Resampler.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPResampler
	: public TOPPBase
{
 public:
	TOPPResampler()
		: TOPPBase("Resampler",
							 "Transforms an LC/MS map into a resampled map or a PNG image.")
	{
	}

 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "input file ");
		setValidFormats_("in", StringList::create("mzML"));
		registerOutputFile_("out", "<file>", "",
												"output file in mzML format");
		setValidFormats_("out", StringList::create("mzML"));

		registerDoubleOption_("sampling_rate", "<rate>", 0.1,
													"New sampling rate in m/z dimension", false);
		setMinFloat_("sampling_rate",0.0);

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

		DoubleReal sampling_rate = getDoubleOption_("sampling_rate");

		LinearResampler lin_resampler;
		Param resampler_param;
		resampler_param.setValue("spacing",sampling_rate);
		lin_resampler.setParameters(resampler_param);

    // resample every scan
    for (Size i = 0; i < exp.size(); ++i)
    {
      lin_resampler.raster(exp[i]);
    }

    //clear meta data because they are no longer meaningful
    exp.clearMetaDataArrays();

    //annotate output with data processing info
		addDataProcessing_(exp,
											 getProcessingInfo_(DataProcessing::DATA_PROCESSING));

    //store output
		f.store(out, exp);

		return EXECUTION_OK;
	}

};


int main( int argc, const char** argv )
{
	TOPPResampler tool;
	return tool.main(argc, argv);
}

/// @endcond
