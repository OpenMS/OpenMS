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

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>

#include "MapAlignerBase.C"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
	 @page TOPP_MapAlignerPoseClustering MapAlignerPoseClustering

	 @brief Corrects retention time distortions between maps, using a pose clustering approach.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MapAlignerPoseClustering \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided @n (or another feature finding algorithm) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinker </td>
		</tr>
	</table>
</CENTER>

  This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them.

	The alignment algorithm implemented here is the pose clustering algorithm as described in doi:10.1093/bioinformatics/btm209. It is used to find an affine transformation, which is further refined by a feature grouping step. This algorithm can be applied to features (featureXML) and peaks (mzML), but it has mostly been developed and tested on features.
	For more details and algorithm-specific parameters (set in the ini file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmPoseClustering "algorithm documentation".

	Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm - see the @p model section of the parameters. This algorithm has been tested mostly with the "linear" model. The different available models are:
	- @ref OpenMS::TransformationModelLinear "linear": Linear model.
	- @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
	- @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

	@see @ref TOPP_MapAlignerIdentification @ref TOPP_MapAlignerSpectrum @ref TOPP_MapAlignerSpectrum

	<B>The command line parameters of this tool are:</B> @n
	@verbinclude TOPP_MapAlignerIdentification.cli

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerPoseClustering
  : public TOPPMapAlignerBase
{

public:
	TOPPMapAlignerPoseClustering()
		: TOPPMapAlignerBase("MapAlignerPoseClustering", "Corrects retention time distortions between maps using a pose clustering approach.")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		String formats = "mzML,featureXML";
		TOPPMapAlignerBase::registerOptionsAndFlags_(formats);
		registerTOPPSubsection_("reference", "Options to define a reference file");
		registerInputFile_("reference:file", "<file>", "", "File to use as reference (same file format as input files required)", false);
		setValidFormats_("reference:file", StringList::create(formats));
		registerIntOption_("reference:index", "<number>", 0, "Use one of the input files as reference ('1' for the first file, etc.).\nIf '0', no explicit reference is set - the algorithm will select a reference.", false);
		setMinInt_("reference:index", 0);
		registerSubsection_("algorithm", "Algorithm parameters section");
		registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
	}

	Param getSubsectionDefaults_(const String& section) const
	{
		if (section == "algorithm")
		{
			MapAlignmentAlgorithmPoseClustering algo;
			return algo.getParameters();
		}
		Param params;
		if (section == "model")
		{
			getModelDefaults_(params, "linear");
		}
		return params;
	}

	ExitCodes main_(int, const char**)
	{
		MapAlignmentAlgorithmPoseClustering algorithm;
		handleReference_(&algorithm);
		return TOPPMapAlignerBase::commonMain_(&algorithm);
	}
};


int main(int argc, const char** argv)
{
  TOPPMapAlignerPoseClustering tool;
  return tool.main(argc, argv);
}

/// @endcond
