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

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

#include "MapAlignerBase.C"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_MapAlignerIdentification MapAlignerIdentification

		@brief Corrects retention time distortions between maps, using information from peptides identified in diferent maps.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ MapAlignerIdentification \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter @n (or another search engine adapter) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMerger </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureLinker </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
		</tr>
	</table>
</CENTER>

	This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them.

	The alignment algorithm implemented here is based on peptide identifications, and thus applicable to files containing peptide IDs (idXML, annotated featureXML/consensusXML). It finds peptide sequences that different input files have in common and uses them as points of correspondence between the inputs. For more details and algorithm-specific parameters (set in the ini file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmIdentification "algorithm documentation".

	Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm - see the @p model section of the parameters. This algorithm has been tested mostly with the "b_spline" model. The different available models are:
	- @ref OpenMS::TransformationModelLinear "linear": Linear model.
	- @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
	- @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

	@see @ref TOPP_MapAlignerApplyTransformation @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapAlignerSpectrum

	<B>The command line parameters of this tool are:</B> @n
	@verbinclude TOPP_MapAlignerIdentification.cli

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerIdentification
  : public TOPPMapAlignerBase
{

public:
	TOPPMapAlignerIdentification()
		: TOPPMapAlignerBase("MapAlignerIdentification", "Corrects retention time distortions between maps based on common peptide identifications.")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		String formats = "featureXML,consensusXML,idXML";
		TOPPMapAlignerBase::registerOptionsAndFlags_(formats);
		registerTOPPSubsection_("reference", "Options to define a reference file");
		registerInputFile_("reference:file", "<file>", "", "File to use as reference", false);
		setValidFormats_("reference:file", StringList::create(formats));
		registerIntOption_("reference:index", "<number>", 0, "Use one of the input files as reference ('1' for the first file, etc.).\nIf '0', no explicit reference is set - the algorithm will use an average of all inputs as reference.", false);
		setMinInt_("reference:index", 0);
		registerSubsection_("algorithm", "Algorithm parameters section");
		registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
	}

	Param getSubsectionDefaults_(const String& section) const
	{
		if (section == "algorithm")
		{
			MapAlignmentAlgorithmIdentification algo;
			return algo.getParameters();
		}
		Param params;
		if (section == "model")
		{
			getModelDefaults_(params, "b_spline");
		}
		return params;
	}

	ExitCodes main_(int, const char**)
	{
		MapAlignmentAlgorithmIdentification algorithm;
		handleReference_(&algorithm);
		return TOPPMapAlignerBase::commonMain_(&algorithm);
	}
};


int main(int argc, const char** argv)
{
  TOPPMapAlignerIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
