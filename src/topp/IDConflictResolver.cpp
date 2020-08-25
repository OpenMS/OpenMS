// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (data annotated with identifications)");
    setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "Output file (data with one peptide identification per feature)");
    setValidFormats_("out", ListUtils::create<String>("featureXML,consensusXML"));
    registerStringOption_("resolve_between_features", "<resolve_between_features>", "off", "A map may contain multiple features with both identical (possibly modified i.e. not stripped) sequence and charge state. The feature with the 'highest intensity' is very likely the most reliable one. When switched on, the filter removes the sequence annotation from the lower intensity features, thereby resolving the multiplicity. Only the most reliable features for each (possibly modified i.e. not stripped) sequence maintain annotated with this peptide sequence.", false);
    setValidStrings_("resolve_between_features", ListUtils::create<String>("off,highest_intensity"));
  }

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out");
    String resolve_between_features = getStringOption_("resolve_between_features");
    
    FileTypes::Type in_type = FileHandler::getType(in);
    
    if (in_type == FileTypes::FEATUREXML) // featureXML
    {
      FeatureMap features;
      FeatureXMLFile().load(in, features);
      
      IDConflictResolverAlgorithm::resolve(features);
      
      if (resolve_between_features=="highest_intensity")
      {
        IDConflictResolverAlgorithm::resolveBetweenFeatures(features);
      }
      
      addDataProcessing_(features, getProcessingInfo_(DataProcessing::FILTERING));
      FeatureXMLFile().store(out, features);
    }
    else // consensusXML
    {
      ConsensusMap consensus;
      ConsensusXMLFile().load(in, consensus);
      
      IDConflictResolverAlgorithm::resolve(consensus);
      
      if (resolve_between_features=="highest_intensity")
      {
        IDConflictResolverAlgorithm::resolveBetweenFeatures(consensus);
      }
      
      addDataProcessing_(consensus, getProcessingInfo_(DataProcessing::FILTERING));
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
