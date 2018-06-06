// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

class TOPPConsensusMapNormalizer :
  public TOPPBase
{

public:
  TOPPConsensusMapNormalizer() :
    TOPPBase("ConsensusMapNormalizer", "Normalizes maps of one consensusXML file")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    addEmptyLine_();
    registerStringOption_("algorithm_type", "<type>", "robust_regression", "The normalization algorithm that is applied. 'robust_regression' scales each map by a fator computed from the ratios of non-differential background features (as determined by the ratio_threshold parameter), 'quantile' performs quantile normalization, 'median' scales all maps to the same median intensity, 'median_shift' shifts the median instead of scaling (WARNING: if you have regular, log-normal MS data, 'median_shift' is probably the wrong choice. Use only if you know what you're doing!)", false, false);
    setValidStrings_("algorithm_type", ListUtils::create<String>("robust_regression,median,median_shift,quantile"));
    registerDoubleOption_("ratio_threshold", "<ratio>", 0.67, "Only for 'robust_regression': the parameter is used to distinguish between non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) and outliers.", false);
    setMinFloat_("ratio_threshold", 0.001);
    setMaxFloat_("ratio_threshold", 1.0);
    registerStringOption_("accession_filter", "<regexp>", "", "Use only features with accessions (partially) matching this regular expression for computing the normalization factors. Useful, e.g., if you have known house keeping proteins in your samples. When this parameter is empty or the regular expression matches the empty string, all features are used (even those without an ID). No effect if quantile normalization is used.", false, true);
    registerStringOption_("description_filter", "<regexp>", "", "Use only features with description (partially) matching this regular expression for computing the normalization factors. Useful, e.g., if you have known house keeping proteins in your samples. When this parameter is empty or the regular expression matches the empty string, all features are used (even those without an ID). No effect if quantile normalization is used.", false, true);
  }

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String algo_type = getStringOption_("algorithm_type");
    String acc_filter = getStringOption_("accession_filter");
    String desc_filter = getStringOption_("description_filter");
    double ratio_threshold = getDoubleOption_("ratio_threshold");

    ConsensusXMLFile infile;
    infile.setLogType(log_type_);
    ConsensusMap map;
    infile.load(in, map);

    //map normalization
    if (algo_type == "robust_regression")
    {
      map.sortBySize();
      vector<double> results = ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(map, ratio_threshold, acc_filter, desc_filter);
      ConsensusMapNormalizerAlgorithmThreshold::normalizeMaps(map, results);
    }
    else if (algo_type == "median")
    {
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, acc_filter, desc_filter);
    }
    else if (algo_type == "median_shift")
    {
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SHIFT, acc_filter, desc_filter);
    }
    else if (algo_type == "quantile")
    {
      if (acc_filter != "" || desc_filter != "")
      {
        LOG_WARN << endl << "NOTE: Accession / description filtering is not supported in quantile normalization mode. Ignoring filters." << endl << endl;
      }
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


int main(int argc, const char ** argv)
{
  TOPPConsensusMapNormalizer tool;
  return tool.main(argc, argv);
}

/// @endcond
