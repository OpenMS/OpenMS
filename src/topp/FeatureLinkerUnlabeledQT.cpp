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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Steffen Sass $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>

#include "FeatureLinkerBase.cpp"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_FeatureLinkerUnlabeledQT FeatureLinkerUnlabeledQT

    @brief Groups corresponding features from multiple maps using a QT clustering approach.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ FeatureLinkerUnlabeledQT \f$ \longrightarrow \f$</td>
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

    Reference:\n
		Weisser <em>et al.</em>: <a href="http://dx.doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

    This tool provides an algorithm for grouping corresponding features in
    multiple runs of label-free experiments. For more details and
    algorithm-specific parameters (set in the ini file) see "Detailed
    Description" in the @ref OpenMS::FeatureGroupingAlgorithmQT "algorithm
    documentation".

    FeatureLinkerUnlabeledQT takes several feature maps (featureXML files) and
    stores the corresponding features in a consensus map (consensusXML file).
    Feature maps can be created from MS experiments (peak data) using one of
    the FeatureFinder TOPP tools.

    @see @ref TOPP_FeatureLinkerUnlabeled @ref TOPP_FeatureLinkerLabeled

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_FeatureLinkerUnlabeledQT.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FeatureLinkerUnlabeledQT.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerUnlabeledQT :
  public TOPPFeatureLinkerBase
{

public:
  TOPPFeatureLinkerUnlabeledQT() :
    TOPPFeatureLinkerBase("FeatureLinkerUnlabeledQT", "Groups corresponding features from multiple maps.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    TOPPFeatureLinkerBase::registerOptionsAndFlags_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    FeatureGroupingAlgorithmQT algo;
    Param p = algo.getParameters();
    return p;
  }

  ExitCodes main_(int, const char **) override
  {
    FeatureGroupingAlgorithmQT algo;
    return TOPPFeatureLinkerBase::common_main_(&algo);
  }

};


int main(int argc, const char ** argv)
{
  TOPPFeatureLinkerUnlabeledQT tool;
  return tool.main(argc, argv);
}

/// @endcond
