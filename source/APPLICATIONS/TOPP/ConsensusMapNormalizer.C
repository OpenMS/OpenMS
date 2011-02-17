// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithm.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ConsensusMapNormalizer ConsensusMapNormalizer

	@brief Normalizes maps of one consensusXML file.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ ConsensusMapNormalizer \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_FeatureLinker </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_SeedListGenerator </td>
		</tr>
	</table>
</CENTER>

This tool normalizes the maps of one consensusXML file by comparing every single map to the map with the highest number of features. Therefore it calculates the pairwise intensity ratios of all features of the two maps. Then it calculates the mean of all ratios within the given threshold, so that outliers will not be considered for the following normalization. Finally it adjusts the intensities of all maps by the corresponding mean of the ratios.

Consensus maps can be created from feature maps (featureXML files) using the @ref TOPP_FeatureLinker.

<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ConsensusMapNormalizer.cli

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPConsensusMapNormalizer
  : public TOPPBase
{

public:
	TOPPConsensusMapNormalizer()
		: TOPPBase("ConsensusMapNormalizer","Normalizes maps of one consensusXML file")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "input file");
		setValidFormats_("in", StringList::create("consensusXML"));
		registerOutputFile_("out", "<file>", "", "output file");
		setValidFormats_("out", StringList::create("consensusXML"));
		addEmptyLine_();
		registerDoubleOption_("ratio_threshold", "<ratio>", 0.67, "threshold for the ratio", false);
		setMinFloat_("ratio_threshold", 0.001);
		setMaxFloat_("ratio_threshold", 1.0);
	}

	ExitCodes main_(int , const char**)
	{
		String in = getStringOption_("in");
		double ratio_threshold = getDoubleOption_("ratio_threshold");

		ConsensusXMLFile infile;
		infile.setLogType(log_type_);
		ConsensusMap map;
		infile.load(in, map);

		map.sortBySize();

		//map normalization
		vector<double> results = ConsensusMapNormalizerAlgorithm::computeCorrelation(map, ratio_threshold);
		ConsensusMapNormalizerAlgorithm::normalizeMaps(map, results);

		String out = getStringOption_("out");
		infile.store(out, map);

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPConsensusMapNormalizer tool;
  return tool.main(argc,argv);
}

/// @endcond
