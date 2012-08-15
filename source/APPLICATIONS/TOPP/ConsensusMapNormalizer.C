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
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher, Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>

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
 
The tool normalizes the intensities of a set of maps (consensusXML file). The following normalization algorithms are available:

- Robust regression: Maps are normalized pair-wise relative to the map with the most features. Given two maps, peptide featues are classified as non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) or outliers. From the non-outliers an average intensity ratio is calculated and used for normalization.

- Median correction: The median of all maps is set to the median of the map with the most features.

- Quantile normalization: Performs an exact quantile normalization if the number of features is equal across all maps. Otherwise, an approximate quantile normalization using resampling is applied.

<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ConsensusMapNormalizer.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_ConsensusMapNormalizer.html

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
    registerStringOption_("algorithm_type", "<type>", "robust_regression", "The normalization algorithm that is applied.", false, false);
    setValidStrings_("algorithm_type", StringList::create("robust_regression,median,quantile"));
    registerDoubleOption_("ratio_threshold", "<ratio>", 0.67, "Only for 'robust_regression': the parameter is used to distinguish between non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) and outliers.", false);
		setMinFloat_("ratio_threshold", 0.001);
		setMaxFloat_("ratio_threshold", 1.0);
	}

	ExitCodes main_(int , const char**)
	{
		String in = getStringOption_("in");
    String out = getStringOption_("out");
    String algo_type = getStringOption_("algorithm_type");
    double ratio_threshold = getDoubleOption_("ratio_threshold");

		ConsensusXMLFile infile;
		infile.setLogType(log_type_);
		ConsensusMap map;
		infile.load(in, map);

    //map normalization
    if (algo_type == "robust_regression")
    {
      map.sortBySize();
      vector<double> results = ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(map, ratio_threshold);
      ConsensusMapNormalizerAlgorithmThreshold::normalizeMaps(map, results);
    }
    else if (algo_type == "median")
    {
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map);
    }
    else if (algo_type == "quantile")
    {
      ConsensusMapNormalizerAlgorithmQuantile::normalizeMaps(map);
    }
    else
    {
      cerr << "Unknown algorithm type  '" << algo_type.c_str() << "'." << endl;
      return ILLEGAL_PARAMETERS;
    }

    //annotate output with data processing info and save output file
    addDataProcessing_(map, getProcessingInfo_(DataProcessing::NORMALIZATION));
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
