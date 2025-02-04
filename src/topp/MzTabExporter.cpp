// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
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
    <th ALIGN = "center"> pot. predecessor tools </td>
       <td VALIGN="middle" ROWSPAN=2> &rarr; MzTabExporter &rarr;</td>
   <th ALIGN = "center"> potential successor tools </td>
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
      registerFlag_("first_run_inference_only", "Does the first IdentificationRun in the file "
                                                "only represent (protein) inference results? If so, read peptide information only "
                                                "from second to last runs.", true);
      registerFlag_("export_all_psms", "Export all PSMs instead of only the best per spectrum", true);
      registerStringList_("opt_columns", "<mods>", {"subfeatures"}, "Add optional columns which are not part of the mzTab standard.", false);
      setValidStrings_("opt_columns", {"subfeatures"});
    }

    ExitCodes main_(int, const char**) override
    {
      // parameter handling
      String in = getStringOption_("in");
      FileTypes::Type in_type = FileHandler::getType(in);

      String out = getStringOption_("out");

      StringList optional_columns = getStringList_("opt_columns");
      bool export_subfeatures = ListUtils::contains(optional_columns, "subfeatures");

      MzTab mztab;

      if (in_type == FileTypes::FEATUREXML)
      {
        // For featureXML we export a "Summary Quantification" file. This means we don't need to report feature quantification values at the assay level
        // but only at the (single) study variable variable level.

        // load featureXML
        FeatureMap feature_map;
        FileHandler().loadFeatures(in, feature_map, {FileTypes::FEATUREXML});

        // calculate coverage
        vector<PeptideIdentification> pep_ids;
        vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();        

        for (Size i = 0; i < feature_map.size(); ++i) // collect all (assigned and unassigned to a feature) peptide ids
        {
          const vector<PeptideIdentification>& pep_ids_bf = feature_map[i].getPeptideIdentifications();
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
          OPENMS_LOG_WARN << "Non-critical exception: " << e.what() << "\n";
        }
        feature_map.setProteinIdentifications(prot_ids);

        mztab = MzTab::exportFeatureMapToMzTab(feature_map, in);
      }

      // export identification data from idXML
      if (in_type == FileTypes::IDXML)
      {
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;
        FileHandler().loadIdentifications(in, prot_ids, pep_ids, {FileTypes::IDXML});

        MzTabFile().store(out,
          prot_ids,
          pep_ids,
          getFlag_("first_run_inference_only"),
          false,
          getFlag_("export_all_psms"));
        return EXECUTION_OK;
      }

      // export identification data from mzIdentML
      if (in_type == FileTypes::MZIDENTML)
      {
        String document_id;
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;
        FileHandler().loadIdentifications(in, prot_ids, pep_ids, {FileTypes::MZIDENTML});

        MzTabFile().store(out,
	        prot_ids,
          pep_ids,
          getFlag_("first_run_inference_only"),
          false,
          getFlag_("export_all_psms"));
        return EXECUTION_OK;
      }

      // export quantification data
      if (in_type == FileTypes::CONSENSUSXML)
      {
        ConsensusMap consensus_map;
        FileHandler().loadConsensusFeatures(in, consensus_map, {FileTypes::CONSENSUSXML});
        IDFilter::removeEmptyIdentifications(consensus_map); // MzTab stream exporter currently doesn't support IDs with empty hits.
        MzTabFile().store(out,
           consensus_map,
           getFlag_("first_run_inference_only"), 
           true, 
           true, 
           export_subfeatures,
           false,
           getFlag_("export_all_psms")); // direct stream to disc
        return EXECUTION_OK;
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

