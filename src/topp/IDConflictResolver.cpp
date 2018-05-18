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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_IDConflictResolver IDConflictResolver

    @brief Resolves ambiguous annotations of features with peptide identifications.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDConflictResolver \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
        </tr>
        <tr>
          <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping algorithm) </td>
          <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
        </tr>
    </table>
    </CENTER>

    The peptide identifications are filtered so that only one identification
    with a single hit (with the best score) is associated to each feature. (If
    two IDs have the same best score, either one of them may be selected.)

    The the filtered identifications are added to the vector of unassigned peptides
    and also reduced to a single best hit.

    This step may be useful before applying @ref TOPP_ProteinQuantifier
    "ProteinQuantifier", because features with ambiguous annotation are not
    considered for the quantification.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_IDConflictResolver.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_IDConflictResolver.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDConflictResolver :
  public TOPPBase
{
public:

  TOPPIDConflictResolver() :
    TOPPBase("IDConflictResolver", "Resolves ambiguous annotations of features with peptide identifications")
  {
  }

private:
  // registerFlag_("f_and_c:make_unique", "Filter for unique sequence/charge state combinations. In a single map, peptides with the same sequence and charge states may appear. If this flag is set only the feature with the highest intensity will pass this filter. (Features without identification will pass automatically, see also: `id:remove_unannotated_features`).", false);
  // TODO: Move this to /home/admin/Code/OpenMS-MFreidank/OpenMS/src/openms/source/ANALYSIS/ID/IDConflictResolverAlgorithm.cpp?
  static FeatureMap makeUnique_(const FeatureMap& feature_map) 
  {
    // Copy (unique) features from `feature_set` over into a fresh FeatureMap.
    FeatureMap unique_features = feature_map;
    unique_features.clear(false);

    // A std::map tracking the set of unique features.
    // Uniqueness criterion/key is a pair <charge, sequence> for each feature.
    typedef std::map<std::pair<Int, AASequence>, const Feature*> FeatureSet;
    FeatureSet feature_set;

    // 1. create a std::map `feature_set` mapping pairs <charge, sequence> to a pointer to 
    // the feature with the highest intensity for this sequence.
    for (const Feature& feature : feature_map) 
    {
      const std::vector<PeptideIdentification>& pep_ids = feature.getPeptideIdentifications();

      if (!pep_ids.empty())
      {
        if (pep_ids.size() != 1) 
        {
          throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Features may contain at most one identification. Run IDConflictResolver first to remove ambiguities!");
        }

        // Assumption: first hit returned by `getHits()` has always highest search engine score.
        // TODO: Is this assumption reasonable or do we need to sort the hits first?
        const std::vector<PeptideHit>& hits = pep_ids.front().getHits();
        if (!hits.empty()) 
        {
          const PeptideHit& highest_score_hit = hits.front();
          // Pair <charge, sequence> of charge of the new feature and the sequence of its highest scoring peptide hit.
          // TODO: const&?
          const std::pair<Int, AASequence> pair = std::make_pair(feature.getCharge(), highest_score_hit.getSequence());

          // If a <charge, sequence> pair is not yet in the FeatureSet or new feature `fm_it` 
          // has higher intensity than its counterpart `feature_set[<charge, sequence>]`
          // store a pointer to `fm_it` in `feature_set`.
          const FeatureSet::iterator feature_in_set = feature_set.find(pair);
          if (feature_in_set != feature_set.end()) 
          {
            if (feature_in_set->second->getIntensity() < feature.getIntensity()) 
            {
              // Replace feature in the set.
              feature_in_set->second = &(feature);
            }
          }
          else 
          {
            // Feature is not yet in our set -- add it.
            feature_set[pair] = &(feature);
          }
        }
      }
      else
      {
        // Maintain features without peptide identifications.
        unique_features.push_back(feature);
      }
    }
    // 2. copy features from `feature_set` into our new FeatureMap.
    for (auto const& element : feature_set) 
    {
      unique_features.push_back(*(element.second);
    }
    return unique_features;
  }



protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (data annotated with identifications)");
    setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "Output file (data with one peptide identification per feature)");
    setValidFormats_("out", ListUtils::create<String>("featureXML,consensusXML"));
  }

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out");
    FileTypes::Type in_type = FileHandler::getType(in);
    if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap features;
      FeatureXMLFile().load(in, features);
      IDConflictResolverAlgorithm::resolve(features);
      addDataProcessing_(features,
                         getProcessingInfo_(DataProcessing::FILTERING));
      FeatureXMLFile().store(out, features);
    }
    else // consensusXML
    {
      ConsensusMap consensus;
      ConsensusXMLFile().load(in, consensus);
      IDConflictResolverAlgorithm::resolve(consensus);
      addDataProcessing_(consensus,
                         getProcessingInfo_(DataProcessing::FILTERING));
      ConsensusXMLFile().store(out, consensus);
    }
    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPIDConflictResolver tool;
  return tool.main(argc, argv);
}

/// @endcond
