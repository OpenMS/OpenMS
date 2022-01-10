// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Dorrestein Lab - University of California San Diego - https://dorresteinlab.ucsd.edu/$
// $Authors: Abinesh Sarvepalli and Louis Felix Nothias$
// $Contributors: Fabian Aicheler and Oliver Alka from Oliver Kohlbacher's group at Tubingen University$
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/GNPSMGFFile.h>

using namespace OpenMS;
using namespace std;

//----------------------------------------------------------
// Doxygen docu
//----------------------------------------------------------
/**
  @page TOPP_GNPSExport GNPSExport
  @brief Export MS/MS data in .MGF format for GNPS (http://gnps.ucsd.edu).
GNPS (Global Natural Products Social Molecular Networking, http://gnps.ucsd.edu) is an open-access knowledge base for community-wide organization and sharing of raw, processed or identified tandem mass (MS/MS) spectrometry data. The GNPS web-platform makes possible to perform spectral library search against public MS/MS spectral libraries, as well as to perform various data analysis such as MS/MS molecular networking, network annotation propagation (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006089), and the Dereplicator-based annotation (https://www.nature.com/articles/nchembio.2219). The GNPS manuscript is available here: https://www.nature.com/articles/nbt.3597
This tool was developed for the Feature Based Molecular Networking (FBMN) workflow on GNPS (https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash2.jsp)
Please cite our preprint:
Nothias, LF., Petras, D., Schmid, R. et al. Feature-based molecular networking in the GNPS analysis environment.
Nat Methods 17, 905–908 (2020). https://doi.org/10.1038/s41592-020-0933-6
See the FBMN workflow documentation here (https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)

In brief, after running an OpenMS "metabolomics" pipeline, the GNPSExport, together with the TextExporter TOPP tool, can be used
on the consensusXML file and the mzML files to generate the files needed for FBMN.
These two files are:
	- The MS/MS spectral data file (.MGF format) which is generated  with the GNPSExport util.
	- The feature quantification table (.TXT format) which is generated with the TextExport util.
  
For each consensusElement in the consensusXML file, the GNPSExport produces one representative consensus
MS/MS spectrum (named peptide annotation in OpenMS jargon) which is appended in the MS/MS spectral file (.MGF file).
An example command is available and described below.
Note that the parameters for the spectral file generation are defined in the GNPSExport INI parameters file, [available with that link](openms_gnpsexport/GNPSExport.ini)).
Representative command:
@code
GNPSExport -ini iniFile-GNPSExport.ini -in_cm filefilter.consensusXML -in_mzml inputFile0.mzML inputFile1.mzML -out GNPSExport_output.mgf
@endcode

Requirements:
	- The IDMapper needs to be run on the featureXML files in order to associate MS2 scan(s) (peptide annotations) with each
	feature for FBMN. An empty idXML or mzid (peptide annotation format) file is needed as an input.
	- The FileFilter has to be run on the consensusXML file, prior to the GNPSExport, in order to remove consensusElements
	without MS2 scans (peptide annotation).
Parameters:
	- Binning (ms2_bin_size): Defines the binning width of fragment ions during the merging of eligible MS/MS spectra.
	- Cosine Score Threshold (merged_spectra:cos_similarity): Defines the necessary pairwise cosine similarity with the highest precursor intensity MS/MS scan.
  - Output Type (output_type):
Options for outputting GNPSExport spectral processing are:
    -# [RECOMMENDED] merged_spectra
      For each consensusElement, the GNPSExport will merge all the eligible MS/MS scans into one representative consensus MS/MS spectrum.
      Eligible MS/MS scans have a pairwise cosine similarity with the MS/MS scan of highest precursor intensity above the Cosine Similarity Threshold.
	    The fragment ions of merged MS/MS scans are binned in m/z (or Da) range defined by the Binning width parameter.
      .
	  -# Most intense: most_intense - For each consensusElement, the GNPSExport will output the most intense MS/MS scan (with the highest precursor ion intensity) as consensus MS/MS spectrum.
      .
Note that mass accuracy and the retention time window for the pairing between MS/MS scans and a LC-MS feature
or consensusElement is defined at the IDMapper tool step.
A representative OpenMS-GNPS workflow would sequentially use these OpenMS TOPP tools:
  1. Input mzML files
  2. Run the @ref TOPP_FeatureFinderMetabo tool on the mzML files.
  3. Run the @ref TOPP_MapAlignerPoseClustering tool on the featureXML files.
  	MapAlignerPoseClustering -in FFM_inputFile0.featureXML FFM_inputFile1.featureXML -out MapAlignerPoseClustering_inputFile0.featureXML MapAlignerPoseClustering_inputFile1.featureXML
  4. Run the @ref TOPP_IDMapper tool on the featureXML and mzML files.
  	IDMapper -id emptyfile.idXML -in MapAlignerPoseClustering_inputFile0.featureXML -spectra:in MapAlignerPoseClustering_inputFile0.mzML -out IDMapper_inputFile0.featureXML
	IDMapper -id emptyfile.idXML -in MapAlignerPoseClustering_inputFile1.featureXML -spectra:in MapAlignerPoseClustering_inputFile1.mzML -out IDMapper_inputFile1.featureXML
  5. Run the @ref TOPP_MetaboliteAdductDecharger on the featureXML files.
  6. Run the @ref TOPP_FeatureLinkerUnlabeledKD tool or FeatureLinkerUnlabeledQT, on the featureXML files and output a consensusXML file.
  	FeatureLinkerUnlabeledKD -in IDMapper_inputFile0.featureXML IDMapper_inputFile1.featureXML -out FeatureLinkerUnlabeledKD.consensusXML
  7. Run the @ref TOPP_FileFilter on the consensusXML file to keep only consensusElements with at least MS/MS scan (peptide identification). 
  	FileFilter -id:remove_unannotated_features -in FeatureLinkerUnlabeledKD.consensusXML -out FileFilter.consensusXML
  8. Run the @ref TOPP_GNPSExport on the "filtered consensusXML file" to export an .MGF file.
  	GNPSExport -ini iniFile-GNPSExport.ini -in_cm filtered.consensusXML -in_mzml inputFile0.mzML inputFile1.mzML -out GNPSExport_output.mgf
  9. Run the @ref TOPP_TextExporter on the "filtered consensusXML file" to export an .TXT file.
  	TextExporter -in FileFilter.consensusXML -out FeatureQuantificationTable.txt
  10. Upload your files to GNPS and run the Feature-Based Molecular Networking workflow. Instructions are here:
https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/
The GitHub for that ProteoSAFe workflow and an OpenMS python wrappers is available here:
https://github.com/Bioinformatic-squad-DorresteinLab/openms-gnps-workflow
An online version of the OpenMS-GNPS pipeline for FBMN running on CCMS server (http://proteomics.ucsd.edu/) is available on GNPS:
https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms/
GNPS (Global Natural Products Social Molecular Networking, https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash2.jsp)
is an open-access knowledge base for community-wide organization and sharing of raw, processed
or identified tandem mass (MS/MS) spectrometry data.
The GNPS web-platform makes possible to perform spectral library search against public MS/MS spectral libraries,
as well as to perform various data analysis such as MS/MS molecular networking, Network Annotation Propagation
Network Annotation Propagation (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006089)
and the DEREPLICATOR (https://www.nature.com/articles/nchembio.2219)
The GNPS paper is available here (https://www.nature.com/articles/nbt.3597)
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
    "Tool to export representative consensus MS/MS scan per consensusElement into a .MGF file format.\nSee the documentation on https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking_with_openms",
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
    setValidFormats_("in_cm", ListUtils::create<String>("consensusXML"));

    registerInputFileList_("in_mzml", "<files>", ListUtils::create<String>(""), "Original mzml files containing the ms2 spectra (aka peptide annotation). \nMust be in order that the consensusXML file maps the original mzML files.");
    setValidFormats_("in_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Output MGF file");
    setValidFormats_("out", ListUtils::create<String>("mgf"));

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

    GNPSMGFFile gnps;
    gnps.setLogType(log_type_);
    gnps.setParameters(getParam_()); // copy tool parameter to library class/algorithm
    gnps.run(consensus_file_path, mzml_file_paths, out);

    return EXECUTION_OK;
  }
};

// the actual main functioned needed to create an executable
int main (int argc, const char** argv)
{
  TOPPGNPSExport tool;
  return tool.main(argc, argv);
}
/// @endcond
