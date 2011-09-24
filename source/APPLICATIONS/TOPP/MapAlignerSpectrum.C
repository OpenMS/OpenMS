// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>

#include "MapAlignerBase.C"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_MapAlignerSpectrum MapAlignerSpectrum

		@brief Corrects retention time distortions between maps by aligning spectra.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MapAlignerSpectrumAlignment \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided @n (or another feature finding algorithm) </td>
		</tr>
	</table>
</CENTER>

	This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them.

	Here, an experimental algorithm based on spectrum alignment is implemented. It is only applicable to peak maps (mzML format).
	For more details and algorithm-specific parameters (set in the ini file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmSpectrumAlignment "algorithm documentation".

	Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm - see the @p model section of the parameters. This algorithm has been tested mostly with the "interpolated" model. The different available models are:
	- @ref OpenMS::TransformationModelLinear "linear": Linear model.
	- @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
	- @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

	@see @ref TOPP_MapAlignerIdentification @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapAlignerApplyTransformation

	<B>The command line parameters of this tool are:</B> @n
	@verbinclude TOPP_MapAlignerSpectrum.cli

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerSpectrum
  : public TOPPMapAlignerBase
{

public:
	TOPPMapAlignerSpectrum()
		: TOPPMapAlignerBase("MapAlignerSpectrum", "Corrects retention time distortions between maps by spectrum alignment.")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		String formats = "mzML";
		TOPPMapAlignerBase::registerOptionsAndFlags_(formats);
		// no support for a reference file yet
		registerSubsection_("algorithm", "Algorithm parameters section");
		registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
	}

	Param getSubsectionDefaults_(const String& section) const
	{
		if (section == "algorithm")
		{
			MapAlignmentAlgorithmSpectrumAlignment algo;
			return algo.getParameters();
		}
		Param params;
		if (section == "model")
		{
			getModelDefaults_(params, "interpolated");
		}
		return params;
	}

	ExitCodes main_(int, const char**)
	{
		MapAlignmentAlgorithmSpectrumAlignment algorithm;
		return TOPPMapAlignerBase::commonMain_(&algorithm);
	}
};


int main(int argc, const char** argv)
{
  TOPPMapAlignerSpectrum tool;
  return tool.main(argc, argv);
}

/// @endcond
