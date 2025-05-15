// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Andreas Bertsch, Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <algorithm>

#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/StopWatch.h>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ProteinInference ProteinInference

@brief Computes a protein identification score based on an aggregation of scores of identified peptides.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> &rarr; ProteinInterference &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_CometAdapter (or other ID engines)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_PeptideIndexer </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
        </tr>
    </table>
</CENTER>

This tool counts and aggregates the scores of peptide sequences that match a protein accession. Only the top PSM for a peptide is used.
By default it also annotates the number of peptides used for the calculation (metavalue "nr_found_peptides") and
can be used for further filtering. 0 probability peptides are counted but ignored in aggregation method "multiplication".

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ProteinInference.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ProteinInference.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPProteinInference :
  public TOPPBase
{
public:
  TOPPProteinInference() :
    TOPPBase("ProteinInference", "Protein inference based on an aggregation of the scores of the identified peptides.")
    {}

protected:

  void registerOptionsAndFlags_() override
  {
    //TODO allow consensusXML version
    registerInputFileList_("in", "<file>", StringList(), "input file(s)");
    setValidFormats_("in", ListUtils::create<String>("idXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("idXML,consensusXML"));
    registerStringOption_("out_type", "<file>", "", "output file type", false);
    setValidStrings_("out_type", ListUtils::create<String>("idXML,consensusXML"));

    //TODO add function to merge based on replicates only. Needs additional exp. design file then.
    registerStringOption_("merge_runs", "<choice>", "all",
                          "If your idXML contains multiple runs, merge them beforehand? Otherwise performs inference separately per run.", false);
    setValidStrings_("merge_runs", ListUtils::create<String>("no,all"));

    registerStringOption_("protein_fdr",
                          "<option>",
                          "false",
                          "Additionally calculate the target-decoy FDR on protein-level after inference", false, false);
    setValidStrings_("protein_fdr", {"true","false"});

    registerStringOption_("conservative_fdr",
                          "<option>",
                          "true",
                          "Use (D+1)/(T) instead of (D+1)/(T+D) for reporting protein FDRs.", false, true);
    setValidStrings_("conservative_fdr", {"true","false"});

    registerStringOption_("picked_fdr",
                          "<option>",
                          "true",
                          "Use picked protein FDRs.", false, true);
    setValidStrings_("picked_fdr", {"true","false"});
    registerStringOption_("picked_decoy_string",
                          "<decoy_string>",
                          "",
                          "If using picked protein FDRs, which decoy string was used? Leave blank for auto-detection.", false, true);
    registerStringOption_("picked_decoy_prefix",
                          "<option>",
                          "prefix",
                          "If using picked protein FDRs, was the decoy string a prefix or suffix? Ignored during auto-detection.", false, true);
    setValidStrings_("picked_decoy_prefix", {"prefix","suffix"});

    // If we support more psms per spectrum, it should be done in the Algorithm class first
    /*registerIntOption_("nr_psms_per_spectrum", "<choice>", 1,
                          "The number of top scoring PSMs per spectrum to consider. 0 means all.", false);
    setMinInt_("nr_psms_per_spectrum", 0);*/

    addEmptyLine_();

    Param merger_with_subsection;
    merger_with_subsection.insert("Merging:", IDMergerAlgorithm().getDefaults());
    registerFullParam_(merger_with_subsection);

    Param algo_with_subsection;
    algo_with_subsection.insert("Algorithm:", BasicProteinInferenceAlgorithm().getDefaults());
    registerFullParam_(algo_with_subsection);
  }


  ExitCodes main_(int, const char**) override
  {
    StopWatch sw;
    sw.start();
    StringList in = getStringList_("in");
    // Merging if specifically asked or multiple files given. If you want to not merge
    // and use multiple files, use a loop
    bool merge_runs = getStringOption_("merge_runs") == "all" || in.size() > 1;
    String out = getStringOption_("out");
    String out_type = getStringOption_("out_type");
    // load identifications
    OPENMS_LOG_INFO << "Loading input..." << std::endl;

    FileTypes::Type in_type = FileHandler::getType(in[0]);

    if (!in.empty() && in_type == FileTypes::CONSENSUSXML)
    {
      if (FileHandler::getTypeByFileName(out) != FileTypes::CONSENSUSXML &&
      FileTypes::nameToType(out_type) != FileTypes::CONSENSUSXML)
      {
        OPENMS_LOG_FATAL_ERROR << "Error: Running on consensusXML requires output as consensusXML. Please change the "
                                  "output type.\n";
      }


      if (in.size() > 1)
      {
        OPENMS_LOG_FATAL_ERROR << "Error: Multiple inputs only supported for idXML\n";
      }

      ConsensusMapMergerAlgorithm cmerge;
      ConsensusMap cmap;
      OPENMS_LOG_INFO << "Loading input..." << std::endl;
      FileHandler().loadConsensusFeatures(in[0], cmap, {FileTypes::CONSENSUSXML});
      OPENMS_LOG_INFO << "Loading input took " << sw.toString() << std::endl;
      sw.clear();

      OPENMS_LOG_INFO << "Merging IDs across runs..." << std::endl;
      cmerge.mergeAllIDRuns(cmap);
      OPENMS_LOG_INFO << "Merging IDs across runs took " << sw.toString() << std::endl;
      sw.clear();

      OPENMS_LOG_INFO << "Aggregating protein scores..." << std::endl;
      BasicProteinInferenceAlgorithm pi;
      pi.setParameters(getParam_().copy("Algorithm:", true));
      pi.run(cmap, cmap.getProteinIdentifications()[0], true);
      OPENMS_LOG_INFO << "Aggregating protein scores took " << sw.toString() << std::endl;
      sw.clear();

      bool calc_protFDR = getStringOption_("protein_fdr") == "true";
      if (calc_protFDR)
      {
        OPENMS_LOG_INFO << "Calculating target-decoy q-values..." << std::endl;
        FalseDiscoveryRate fdr;
        Param fdrparam = fdr.getParameters();
        fdrparam.setValue("conservative", getStringOption_("conservative_fdr"));
        fdrparam.setValue("add_decoy_proteins","true");
        fdr.setParameters(fdrparam);
        if (getStringOption_("picked_fdr") == "true")
        {
          fdr.applyPickedProteinFDR(cmap.getProteinIdentifications()[0], getStringOption_("picked_decoy_string"), getStringOption_("picked_decoy_prefix") == "prefix");
        }
        else
        {
          fdr.applyBasic(cmap.getProteinIdentifications()[0], true);
        }
      }

      OPENMS_LOG_INFO << "Storing output..." << std::endl;
      sw.start();
      // write output
      FileHandler().storeConsensusFeatures(out, cmap, {FileTypes::CONSENSUSXML});
      OPENMS_LOG_INFO << "Storing output took " << sw.toString() << std::endl;
      sw.stop();

    }
    else //----------- IdXML --------------------------
    {
      vector<ProteinIdentification> inferred_protein_ids{1};
      vector<PeptideIdentification> inferred_peptide_ids;

      FileHandler f;
      if (merge_runs)
      {
        //TODO allow keep_best_pepmatch_only option during merging (Peptide-level datastructure would help a lot,
        // otherwise you need to build a map of peptides everytime you want to quickly check if the peptide is already
        // present)
        //TODO allow experimental design aware merging
        IDMergerAlgorithm merger{String("all_merged")};
        merger.setParameters(getParam_().copy("Merging:", true));

        for (const auto &idfile : in)
        {
          vector<ProteinIdentification> protein_ids;
          vector<PeptideIdentification> peptide_ids;
          f.loadIdentifications(idfile, protein_ids, peptide_ids, {FileTypes::IDXML});
          merger.insertRuns(std::move(protein_ids), std::move(peptide_ids));
        }
        merger.returnResultsAndClear(inferred_protein_ids[0], inferred_peptide_ids);
      }
      else
      {
        f.loadIdentifications(in[0], inferred_protein_ids, inferred_peptide_ids, {FileTypes::IDXML});
      }
      OPENMS_LOG_INFO << "Loading input took " << sw.toString() << std::endl;
      sw.reset();

      // groups will be reannotated or scores will not make sense anymore -> delete
      inferred_protein_ids[0].getIndistinguishableProteins().clear();

      OPENMS_LOG_INFO << "Aggregating protein scores..." << std::endl;
      BasicProteinInferenceAlgorithm pi;
      pi.setParameters(getParam_().copy("Algorithm:", true));
      pi.run(inferred_peptide_ids, inferred_protein_ids);
      OPENMS_LOG_INFO << "Aggregating protein scores took " << sw.toString() << std::endl;
      sw.clear();

      bool calc_protFDR = getStringOption_("protein_fdr") == "true";
      if (calc_protFDR)
      {
        OPENMS_LOG_INFO << "Calculating target-decoy q-values..." << std::endl;
        FalseDiscoveryRate fdr;
        Param fdrparam = fdr.getParameters();
        fdrparam.setValue("conservative", getStringOption_("conservative_fdr"));
        fdrparam.setValue("add_decoy_proteins","true");
        fdr.setParameters(fdrparam);
        if (getStringOption_("picked_fdr") == "true")
        {
          fdr.applyPickedProteinFDR(inferred_protein_ids[0], getStringOption_("picked_decoy_string"), getStringOption_("picked_decoy_prefix") == "prefix");
        }
        else
        {
          fdr.applyBasic(inferred_protein_ids[0], true);
        }
      }

      OPENMS_LOG_INFO << "Storing output..." << std::endl;
      sw.start();
      // write output
      FileHandler().storeIdentifications(out, inferred_protein_ids, inferred_peptide_ids, {FileTypes::IDXML});
      OPENMS_LOG_INFO << "Storing output took " << sw.toString() << std::endl;
      sw.stop();
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPProteinInference tool;
  return tool.main(argc, argv);
}

/// @endcond
