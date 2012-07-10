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
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>

#include "FeatureLinkerBase.C"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FeatureLinkerLabeled FeatureLinkerLabeled

	@brief Groups corresponding isotope-labeled features in a feature map.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureLinkerLabeled \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureFinderCentroided @n (or another feature detection algorithm) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
		</tr>
	</table>
</CENTER>


	This tool provides an algorithm for grouping corresponding features in isotope-labeled experiments. For more details and algorithm-specific parameters (set in the ini file) see "Detailed Description" in the @ref OpenMS::FeatureGroupingAlgorithmLabeled "algorithm documentation".

	FeatureLinkerLabeled takes one feature map (featureXML file) and stores the corresponding features in a consensus map (consensusXML file). Feature maps can be created from MS experiments (peak data) using one of the FeatureFinder TOPP tools.

	@see @ref TOPP_FeatureLinkerUnlabeled @ref TOPP_FeatureLinkerUnlabeledQT

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FeatureLinkerLabeled.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_FeatureLinkerLabeled.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerLabeled
  : public TOPPFeatureLinkerBase
{

public:
	TOPPFeatureLinkerLabeled()
		: TOPPFeatureLinkerBase("FeatureLinkerLabeled", "Groups corresponding isotope-labeled features in a feature map.")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "Input file", true);
		setValidFormats_("in", StringList::create("featureXML"));
		registerOutputFile_("out", "<file>", "", "Output file", true);
		setValidFormats_("out", StringList::create("consensusXML"));
		registerSubsection_("algorithm", "Algorithm parameters section");
	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		FeatureGroupingAlgorithmLabeled algo;
		Param p = algo.getParameters();
		return p;
	}

	ExitCodes main_(int , const char**)
	{
		FeatureGroupingAlgorithmLabeled algo;
		return TOPPFeatureLinkerBase::common_main_(&algo, true);
	}
};


int main(int argc, const char** argv)
{
  TOPPFeatureLinkerLabeled tool;
  return tool.main(argc, argv);
}

/// @endcond
