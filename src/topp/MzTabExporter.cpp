// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>

#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_MzTabExporter MzTabExporter

   @brief This application converts several %OpenMS XML formats (featureXML, consensusXML, and idXML) to mzTab.

  <CENTER>
    <table>
     <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
         <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MzTabExporter \f$ \longrightarrow \f$</td>
     <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> Any tool producing one of the input formats </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> External tools (MS Excel, OpenOffice, Notepad)</td>
    </tr>
   </table>
  </CENTER>

  See the mzTab specification for details on the format.

  @experimental This algorithm and underlying format is work in progress and might change.

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MzTabExporter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_MzTabExporter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

namespace OpenMS
{
  class TOPPMzTabExporter :
    public TOPPBase
  {
public:
    TOPPMzTabExporter() :
      TOPPBase("MzTabExporter", "Exports various XML formats to an mzTab file.")
    {
    }

protected:

    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "", "Input files used to generate the mzTab file.", false);
      setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML,idXML,mzid"));
      registerOutputFile_("out", "<file>", "", "Output file (mzTab)", true);
      setValidFormats_("out", ListUtils::create<String>("mzTab"));
      registerFlag_("first_run_inference_only", "(idXML/mzid only): Does the first IdentificationRun in the file "
                                                "only represent inference results? Read peptide information only "
                                                "from second to last runs.");
          }

    ExitCodes main_(int, const char**) override
    {
      // parameter handling
      String in = getStringOption_("in");
      FileTypes::Type in_type = FileHandler().getType(in);

      String out = getStringOption_("out");

      MzTab mztab;

      if (in_type == FileTypes::FEATUREXML)
      {
        // For featureXML we export a "Summary Quantification" file. This means we don't need to report feature quantification values at the assay level
        // but only at the (single) study variable variable level.

        // load featureXML
        FeatureMap feature_map;
        FeatureXMLFile f;
        f.load(in, feature_map);

        // calculate coverage
        vector<PeptideIdentification> pep_ids;
        vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();        

        for (Size i = 0; i < feature_map.size(); ++i) // collect all (assigned and unassigned to a feature) peptide ids
        {
          vector<PeptideIdentification> pep_ids_bf = feature_map[i].getPeptideIdentifications();
          pep_ids.insert(pep_ids.end(), pep_ids_bf.begin(), pep_ids_bf.end());
        }

        pep_ids.insert(pep_ids.end(), feature_map.getUnassignedPeptideIdentifications().begin(), feature_map.getUnassignedPeptideIdentifications().end());

        try // might throw Exception::MissingInformation()
        {
          for (Size i = 0; i < prot_ids.size(); ++i)
          {
            prot_ids[i].computeCoverage(pep_ids);
          }
        }
        catch (Exception::MissingInformation& e)
        {
          LOG_WARN << "Non-critical exception: " << e.what() << "\n";
        }
        feature_map.setProteinIdentifications(prot_ids);

        mztab = MzTab::exportFeatureMapToMzTab(feature_map, in);
      }

      // export identification data from idXML
      if (in_type == FileTypes::IDXML)
      {
        String document_id;
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;
        IdXMLFile().load(in, prot_ids, pep_ids, document_id);
        mztab = MzTab::exportIdentificationsToMzTab(prot_ids, pep_ids, in, getFlag_("first_run_inference_only"));
      }

      // export identification data from mzIdentML
      if (in_type == FileTypes::MZIDENTML)
      {
        String document_id;
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;
        MzIdentMLFile().load(in, prot_ids, pep_ids);
        mztab = MzTab::exportIdentificationsToMzTab(prot_ids, pep_ids, in, getFlag_("first_run_inference_only"));
      }

      // export quantification data
      if (in_type == FileTypes::CONSENSUSXML)
      {
        ConsensusMap consensus_map;
        ConsensusXMLFile c;
        c.load(in, consensus_map);
        mztab = MzTab::exportConsensusMapToMzTab(consensus_map, in, true, true);
      }

      MzTabFile().store(out, mztab);
      return EXECUTION_OK;
    }
  };
} //namespace OpenMS

#pragma clang diagnostic pop

int main(int argc, const char** argv)
{
  TOPPMzTabExporter t;
  return t.main(argc, argv);
}

/// @endcond

