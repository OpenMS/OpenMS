// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Dorrestein Lab - University of California San Diego - https://dorresteinlab.ucsd.edu/$
// $Authors: Abinesh Sarvepalli and Louis Felix Nothias$
// $Contributors: Fabian Aicheler and Oliver Alka from Oliver Kohlbacher's group at Tubingen University$
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/GNPSMetaValueFile.h>
#include <OpenMS/FORMAT/GNPSMGFFile.h>
#include <OpenMS/FORMAT/GNPSQuantificationFile.h>
#include <OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

//----------------------------------------------------------
// Doxygen docu
//----------------------------------------------------------
/**
@page TOPP_GNPSExport GNPSExport
@brief Export MS/MS data in .MGF format for GNPS (http://gnps.ucsd.edu).

GNPS (Global Natural Products Social Molecular Networking, http://gnps.ucsd.edu) is an open-access knowledge base for community-wide organization and sharing of raw, processed or identified tandem mass (MS/MS) spectrometry data. The GNPS web-platform makes it possible to perform spectral library search against public MS/MS spectral libraries, as well as to perform various data analysis such as MS/MS molecular networking, network annotation propagation, and the Dereplicator-based annotation. The GNPS manuscript is available here: https://www.nature.com/articles/nbt.3597
This tool was developed for the Feature Based Molecular Networking (FBMN) (https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/) and Ion Identity Molecular Networking (IIMN) (https://ccms-ucsd.github.io/GNPSDocumentation/fbmn-iin/) workflows.

Please cite:
Nothias, L.-F., Petras, D., Schmid, R. et al. [Feature-based molecular networking in the GNPS analysis environment](https://www.nature.com/articles/s41592-020-0933-6). Nat. Methods 17, 905–908 (2020).

In brief, after running an OpenMS metabolomics pipeline, the <b>GNPSExport</b> TOPP tool can be used on the consensusXML file and the mzML files to generate the files needed for FBMN and IIMN.
Those files are:
- A <b>MS/MS spectral data file</b> (.MGF format).
- A <b>feature quantification table</b> (.TXT format). (https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/#feature-quantification-table)
- A <b>supplementary pairs table</b> (.CSV format) required for IIMN. (https://ccms-ucsd.github.io/GNPSDocumentation/fbmn-iin/#supplementary-pairs)
- A <b>meta value table</b> (.TSV format). (https://ccms-ucsd.github.io/GNPSDocumentation/metadata/)

A representative OpenMS-GNPS workflow would use the following OpenMS TOPP tools sequentially:
- Input mzML files
- Run the @ref TOPP_FeatureFinderMetabo tool on the mzML files.
- Run MetaboliteAdductDecharger on the featureXML files (optional, for Ion Identity Molecular Networking).
- Run the @ref TOPP_MapAlignerPoseClustering tool on the featureXML files.
@code
	MapAlignerPoseClustering -in FFM_inputFile0.featureXML FFM_inputFile1.featureXML -out MapAlignerPoseClustering_inputFile0.featureXML MapAlignerPoseClustering_inputFile1.featureXML -trafo_out MapAlignerPoseClustering_inputFile0.trafoXML MapAlignerPoseClustering_inputFile1.trafoXML
@endcode
- Run the @ref TOPP_MapRTTransformer tool on the mzML files to transform retention times based on the feature map alignment by @ref TOPP_MapAlignerPoseClustering.
@code
  MapRTTransformer -in inputFile0.mzML -out MapRTTransformer_inputFile0.mzML -trafo_in MapAlignerPoseClustering_inputFile0.trafoXML
  MapRTTransformer -in inputFile1.mzML -out MapRTTransformer_inputFile1.mzML -trafo_in  MapAlignerPoseClustering_inputFile1.trafoXML
@endcode
- Run the @ref TOPP_IDMapper tool on the featureXML and mzML files.
@code
  IDMapper -id emptyfile.idXML -in MapAlignerPoseClustering_inputFile0.featureXML -spectra:in MapRTTransformer_inputFile0.mzML -out IDMapper_inputFile0.featureXML
	IDMapper -id emptyfile.idXML -in MapAlignerPoseClustering_inputFile1.featureXML -spectra:in MapRTTransformer_inputFile1.mzML -out IDMapper_inputFile1.featureXML
@endcode
- Run the @ref TOPP_MetaboliteAdductDecharger tool on the featureXML files.
- Run the @ref TOPP_FeatureLinkerUnlabeledKD tool or FeatureLinkerUnlabeledQT, on the featureXML files and output a consensusXML file.
@code
  	FeatureLinkerUnlabeledKD -in IDMapper_inputFile0.featureXML IDMapper_inputFile1.featureXML -out FeatureLinkerUnlabeledKD.consensusXML
@endcode
- Run the @ref TOPP_FileFilter on the consensusXML file to keep only consensusElements with at least MS/MS scan (peptide identification). 
@code
  	FileFilter -id:remove_unannotated_features -in FeatureLinkerUnlabeledKD.consensusXML -out FileFilter.consensusXML
@endcode
- Run the @ref TOPP_GNPSExport on the "filtered consensusXML file" to export an .MGF file. For each consensusElement in the consensusXML file, the GNPSExport command produces one representative consensus MS/MS spectrum (named peptide annotation in OpenMS jargon) which is appended in the MS/MS spectral file (.MGF file).
(Note that the parameters for the spectral file generation are defined in the GNPSExport INI parameters file, available here: https://ccms-ucsd.github.io/GNPSDocumentation/openms_gnpsexport/GNPSExport.ini
@code 
	GNPSExport -in_cm filtered.consensusXML -in_mzml MapRTTransformer_inputFile0.mzML MapRTTransformer_inputFile1.mzML -out GNPSExport_output.mgf -out_quantification FeatureQuantificationTable.txt -out_pairs SupplementaryPairsTable.csv -out_meta_values MetaValues.tsv
@endcode
- Upload your files to GNPS and run the Feature-Based Molecular Networking workflow. Instructions can be found here: https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/

The GitHub page for the ProteoSAFe workflow and the OpenMS python wrappers is available here: https://github.com/Bioinformatic-squad-DorresteinLab/openms-gnps-workflow
An online version of the OpenMS-GNPS pipeline for FBMN running on CCMS server (http://proteomics.ucsd.edu/) is available here: https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms/
  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_GNPSExport.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_GNPSExport.html
 */

class TOPPGNPSExport : public TOPPBase
{
public:
  TOPPGNPSExport() : TOPPBase(
    "GNPSExport",
    "Export representative consensus MS/MS scan per consensusElement into a .MGF file format.\nSee the documentation on https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms",
    true,
    {
      {
        "Nothias L.F. et al.", // authors
        "Feature-based Molecular Networking in the GNPS Analysis Environment", // title
        "bioRxiv 812404 (2019)", // when_where
        "10.1101/812404" // doi
      }
    }
  ) {}

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_cm", "<file>", "", "Input consensusXML file containing only consensusElements with \"peptide\" annotations.");
    setValidFormats_("in_cm", {"consensusXML"});

    registerInputFileList_("in_mzml", "<files>", ListUtils::create<String>(""), "Original mzml files containing the ms2 spectra (aka peptide annotation). \nMust be in order that the consensusXML file maps the original mzML files.");
    setValidFormats_("in_mzml", {"mzML"});

    registerOutputFile_("out", "<file>", "", "Output MGF file.");
    setValidFormats_("out", {"mgf"});

    registerOutputFile_("out_quantification", "<file>", "", "Output feature quantification table.");
    setValidFormats_("out_quantification", {"txt"});

    registerOutputFile_("out_pairs", "<file>", "", "Output supplementary pairs table for IIMN.", false);
    setValidFormats_("out_pairs", {"csv"});

    registerOutputFile_("out_meta_values", "<file>", "", "Output meta value file.", false);
    setValidFormats_("out_meta_values", {"tsv"});

    addEmptyLine_();

    registerFullParam_(GNPSMGFFile().getDefaults());
  }

  // the main function is called after all parameters are read
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String consensus_file_path(getStringOption_("in_cm"));
    StringList mzml_file_paths = getStringList_("in_mzml");
    String out(getStringOption_("out"));
    String out_quantification(getStringOption_("out_quantification"));
    String out_pairs(getStringOption_("out_pairs"));
    String out_meta(getStringOption_("out_meta_values"));

    // load ConsensusMap from file
    ConsensusMap cm;
    FileHandler().loadConsensusFeatures(consensus_file_path, cm, {FileTypes::CONSENSUSXML});

    // if at least one of the features has an annotation for Constants::UserParam::IIMN_LINKED_GROUPS, annotate ConsensusMap for IIMN
    for (const auto& f: cm)
    {
      if (f.metaValueExists(Constants::UserParam::IIMN_LINKED_GROUPS))
      {
        IonIdentityMolecularNetworking::annotateConsensusMap(cm);
        break;
      }
    }


    GNPSMGFFile gnps;
    gnps.setLogType(log_type_);
    gnps.setParameters(getParam_()); // copy tool parameter to library class/algorithm
    gnps.store(consensus_file_path, mzml_file_paths, out);

    if (!out_pairs.empty()) IonIdentityMolecularNetworking::writeSupplementaryPairTable(cm, out_pairs);
    if (!out_quantification.empty()) GNPSQuantificationFile::store(cm, out_quantification);
    if (!out_meta.empty()) GNPSMetaValueFile::store(cm, out_meta);
    
    return EXECUTION_OK;
  }
};

// the actual main functioned needed to create an executable
int main (int argc, const char** argv)
{
  TOPPGNPSExport tool;
  return tool.main(argc, argv);
}
