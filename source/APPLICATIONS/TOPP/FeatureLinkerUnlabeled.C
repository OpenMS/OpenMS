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
// $Authors: Marc Sturm, Clemens Groepl, Steffen Sass $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>

#include "FeatureLinkerBase.C"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FeatureLinkerUnlabeled FeatureLinkerUnlabeled

	@brief Groups corresponding features from multiple maps.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ FeatureLinkerUnlabeled \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided @n (or another feature detection algorithm) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_MapAlignerPoseClustering @n (or another map alignment algorithm) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_SeedListGenerator </td>
		</tr>
	</table>
</CENTER>


	This tool provides an algorithm for grouping corresponding features in multiple runs of label-free experiments. For more details and algorithm-specific parameters (set in the ini file) see "Detailed Description" in the @ref OpenMS::FeatureGroupingAlgorithmUnlabeled "algorithm documentation".

	FeatureLinkerUnlabeled takes several feature maps (featureXML files) and stores the corresponding features in a consensus map (consensusXML file). Feature maps can be created from MS experiments (peak data) using one of the FeatureFinder TOPP tools.

	@see @ref TOPP_FeatureLinkerUnlabeledQT @ref TOPP_FeatureLinkerLabeled

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FeatureLinkerUnlabeled.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerUnlabeled
  : public TOPPFeatureLinkerBase
{

public:
	TOPPFeatureLinkerUnlabeled()
		: TOPPFeatureLinkerBase("FeatureLinkerUnlabeled", "Groups corresponding features from multiple maps.")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		TOPPFeatureLinkerBase::registerOptionsAndFlags_();
		registerSubsection_("algorithm", "Algorithm parameters section");
	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		FeatureGroupingAlgorithmUnlabeled algo;
		Param p = algo.getParameters();
		return p;
	}

	ExitCodes main_(int , const char**)
	{
		FeatureGroupingAlgorithmUnlabeled algo;
		return TOPPFeatureLinkerBase::common_main_(&algo);
	}
};


int main(int argc, const char** argv)
{
  TOPPFeatureLinkerUnlabeled tool;
  return tool.main(argc, argv);
}

/// @endcond
