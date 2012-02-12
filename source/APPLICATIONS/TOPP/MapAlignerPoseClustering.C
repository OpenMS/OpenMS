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
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>

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
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
		</tr>
	</table>
</CENTER>

  This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them. Retention time adjustment may be necessary to correct for chromatography differences e.g. before data from multiple LC-MS runs can be combined (feature grouping), or when one run should be annotated with peptide identifications obtained in a different run.

	All map alignment tools (MapAligner...) collect retention time data from the input files and - by fitting a model to this data - compute transformations that map all runs to a common retention time scale. They can apply the transformations right away and return output files with aligned time scales (parameter @p out), and/or return descriptions of the transformations in trafoXML format (parameter @p trafo_out). Transformations stored as trafoXML can be applied to arbitrary files with the @ref TOPP_MapRTTransformer tool.

	The map alignment tools differ in how they obtain retention time data for the modeling of transformations, and consequently what types of data they can be applied to. The alignment algorithm implemented here is the pose clustering algorithm as described in doi:10.1093/bioinformatics/btm209. It is used to find an affine transformation, which is further refined by a feature grouping step. This algorithm can be applied to features (featureXML) and peaks (mzML), but it has mostly been developed and tested on features.
	For more details and algorithm-specific parameters (set in the INI file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmPoseClustering "algorithm documentation".

	@see @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapAlignerSpectrum @ref TOPP_MapRTTransformer

	Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm. This algorithm has been tested mostly with the "linear" model. The different available models are:
	- @ref OpenMS::TransformationModelLinear "linear": Linear model.
	- @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
	- @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

	The following parameters control the modeling of RT transformations (they can be set in the "model" section of the INI file):
	@htmlinclude OpenMS_MapAlignerPoseClusteringModel.parameters @n

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
		if (section == "model")
		{
			return getModelDefaults("linear");
		}
		return Param(); // shouldn't happen
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
