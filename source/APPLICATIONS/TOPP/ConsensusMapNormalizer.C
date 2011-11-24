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

	@brief Normalization of intensities in a set of maps using robust regression.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ ConsensusMapNormalizer \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
		</tr>
	</table>
</CENTER>
 
The tool normalizes the intensities of a set of maps (consensusXML file) using robust regression. Maps are normalized pair-wise relative to the map with the most features. Given two maps, peptide featues are classified as non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) or outliers. From the non-outliers an average intensity ratio is calculated and used for normalization.

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
		registerDoubleOption_("ratio_threshold", "<ratio>", 0.67, "The parameter is used to distinguish between non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) and outliers.", false);
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

		//annotate output with data processing info
		addDataProcessing_(map, getProcessingInfo_(DataProcessing::NORMALIZATION));

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
