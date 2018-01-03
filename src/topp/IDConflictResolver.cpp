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

protected:

  // compare peptide IDs by score of best hit (hits must be sorted first!)
  // (note to self: the "static" is necessary to avoid cryptic "no matching
  // function" errors from gcc when the comparator is used below)
  static bool compareIDsSmallerScores_(const PeptideIdentification & left,
                          const PeptideIdentification & right)
  {
    // if any of them is empty, the other is considered "greater"
    // independent of the score in the first hit
    if (left.getHits().empty()) return true;
    if (right.getHits().empty()) return false;
    if (left.getHits()[0].getScore() < right.getHits()[0].getScore())
    {
      return true;
    }
    return false;
  }

  void resolveConflict_(vector<PeptideIdentification> & peptides)
  {
    if (peptides.empty()) return;

    for (vector<PeptideIdentification>::iterator pep_id = peptides.begin();
         pep_id != peptides.end(); ++pep_id)
    {
      pep_id->sort();
    }
    vector<PeptideIdentification>::iterator pos;
    if (peptides[0].isHigherScoreBetter())     // find highest-scoring ID
    {
      pos = max_element(peptides.begin(), peptides.end(), compareIDsSmallerScores_);
    }
    else  // find lowest-scoring ID
    {
      pos = min_element(peptides.begin(), peptides.end(), compareIDsSmallerScores_);
    }
    peptides[0] = *pos;
    peptides.resize(1);

    if (!peptides[0].getHits().empty())
    {
      // remove all but the best hit:
      vector<PeptideHit> best_hit(1, peptides[0].getHits()[0]);
      peptides[0].setHits(best_hit);
    }
  }

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
      for (FeatureMap::Iterator feat_it = features.begin();
           feat_it != features.end(); ++feat_it)
      {
        resolveConflict_(feat_it->getPeptideIdentifications());
      }
      addDataProcessing_(features,
                         getProcessingInfo_(DataProcessing::FILTERING));
      FeatureXMLFile().store(out, features);
    }
    else     // consensusXML
    {
      ConsensusMap consensus;
      ConsensusXMLFile().load(in, consensus);
      for (ConsensusMap::Iterator cons_it = consensus.begin();
           cons_it != consensus.end(); ++cons_it)
      {
        resolveConflict_(cons_it->getPeptideIdentifications());
      }
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
