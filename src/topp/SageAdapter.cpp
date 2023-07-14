// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <regex>

#include <QStringList>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_SageAdapter SageAdapter

    @brief Identifies peptides in MS/MS spectra via sage.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; SageAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em Sage must be installed before this wrapper can be used.     

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_SageAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_SageAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


/*
./bin/SageAdapter -in /nfs/wsi/abi/projects/proteomics/yasset_iPRG2015/JD_06232014_sample1_A.mzML -out ~/tmp -sage_executable ~/OMS/sage/sage-v0.11.2-x86_64-unknown-linux-musl/sage -database /nfs/wsi/abi/projects/proteomics/yasset_iPRG2015/iPRG2015_decoy.fasta -config_file ~/OMS/sage/sage-v0.11.2-x86_64-unknown-linux-musl/config.json
*/

class TOPPSageAdapter :
  public SearchEngineBase
{
public: 
  TOPPSageAdapter() :
    SearchEngineBase("SageAdapter", "Annotates MS/MS spectra using Sage.", true,
             {
                 {"Michael Lazear",
                 "https://lazear.github.io/sage/",
                 "",
                 ""}
             })
  {
  }

protected:
  // values here are sage parameter that can be changed via TOPP tool parameter.
  // They will be pasted into the config_template below. E.g. bucket_size at tag ##bucket_size##
  size_t bucket_size = 32768;
  size_t min_len = 5; 
  size_t max_len = 50; 
  size_t missed_cleavages = 2;
  double fragment_min_mz = 200.0;
  double fragment_max_mz = 2000.0;
  double peptide_min_mass = 500.0;
  double peptide_max_mass = 5000.0;
  size_t min_ion_index = 2;
  size_t max_variable_mods = 2;
  std::string precursor_tol_unit = "ppm";
  double precursor_tol_left = -6.0;
  double precursor_tol_right = 6.0;
  std::string fragment_tol_unit = "ppm";
  double fragment_tol_left = -20.0;
  double fragment_tol_right = 20.0;
  IntList isotope_errors = {-1, 3};
  size_t min_matched_peaks = 6;
  size_t report_psms = 1;
  
  /* TODO: enzyme_details have format:
    "cleave_at": "KR",
    "restrict": "P",
    "c_terminal": true
    and will be filled from OpenMS enzyme information
  */

  /* TODO: static_mods have format:
          "^": 304.207,
          "K": 304.207,
          "C": 57.0215
    and will be filled from OpenMS mod information
  */

  /* TODO: variable_mods have format:
          "M": [15.9949],
          "^Q": [-17.026549],
          "^E": [-18.010565],
          "$": [49.2, 22.9],
          "[": 42.0,
          "]": 111.0
    and will be filled from OpenMS mod information
  */

  std::string config_template = R"(
    {
      "database": {
        "bucket_size": ##bucket_size##,
        "enzyme": {
          "missed_cleavages": 2,
          "min_len": ##min_len##,
          "max_len": ##max_len##,
          ##enzyme_details##
        },
        "fragment_min_mz": ##fragment_min_mz##,
        "fragment_max_mz": ##fragment_max_mz##2000.0,
        "peptide_min_mass": ##peptide_min_mass##500.0,
        "peptide_max_mass": ##peptide_max_mass##5000.0,
        "ion_kinds": ["b", "y"],
        "min_ion_index": ##min_ion_index##,
        "static_mods": {
          ##static_mods##
        },
        "variable_mods": {
          ##variable_mods##
        },
        "max_variable_mods": ##max_variable_mods##,
        "generate_decoys": false,
      },
      "precursor_tol": {
        "##precursor_tol_unit##": [
          ##precursor_tol_left##,
          ##precursor_tol_right##
        ]
      },
      "fragment_tol": {
        "##fragment_tol_unit##": [
        ##fragment_tol_left##,
        ##fragment_tol_right##
        ]
      },
      "isotope_errors": [
        ##isotope_errors##
      ],
      "deisotope": false,
      "chimera": false,
      "wide_window": false,
      "predict_rt": false,
      "min_peaks": 15,
      "max_peaks": 150,
      "min_matched_peaks": ##min_matched_peaks##,
      "report_psms": ##report_psms##,
      ]       
    }
  )";

  std::tuple<std::string, std::string, std::string> getVersionNumber_(const std::string& multi_line_input)
  {
      std::regex version_regex("Version ([0-9]+)\\.([0-9]+)\\.([0-9]+)");

      std::sregex_iterator it(multi_line_input.begin(), multi_line_input.end(), version_regex);
      std::smatch match = *it;
      std::cout << "Found Sage version string: " << match.str() << std::endl;      
          
      return make_tuple(it->str(1), it->str(2), it->str(3)); // major, minor, patch
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", { "mzML" } );

    registerOutputFile_("out", "<file>", "", "Single output file containing all search results.", true, false);
    setValidFormats_("out", { "idXML" } );

    registerInputFile_("database", "<file>", "", "FASTA file", true, false, {"skipexists"});
    setValidFormats_("database", { "FASTA" } );
    registerInputFile_("sage_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      "sage", // this is the name on ALL platforms currently...
      "The Sage executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", true, false, {"is_executable"});

    registerInputFile_("config_file", "<file>", "", "Default Sage config file.", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("config_file", ListUtils::create<String>("json"));

    registerStringOption_("decoy_prefix", "<prefix>", "DECOY_", "Prefix on protein accession used to distinguish decoy from target proteins.", false, false);
    registerIntOption_("batch_size", "<int>", 0, "Number of files to load and search in parallel (default = # of CPUs/2)", false, false);
    
    // register peptide indexing parameter (with defaults for this search engine) TODO: check if search engine defaults are needed
    registerPeptideIndexingParameter_(PeptideIndexing().getParameters());     
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    // do this early, to see if Sage is installed
    String sage_executable = getStringOption_("sage_executable");
    String proc_stdout, proc_stderr;
    TOPPBase::ExitCodes exit_code = runExternalProcess_(sage_executable.toQString(), QStringList() << "--help", proc_stdout, proc_stderr, "");
    auto major_minor_patch = getVersionNumber_(proc_stdout);
    String sage_version = "sage (" + std::get<0>(major_minor_patch) + "." + std::get<1>(major_minor_patch) + "." + std::get<2>(major_minor_patch) + ")";
    
    //-------------------------------------------------------------
    // run sage
    //-------------------------------------------------------------
    StringList input_files = getStringList_("in");
    String output_file = getStringOption_("out");
    String output_folder = File::isDirectory(output_file) ? output_file : File::path(output_file);
    String fasta_file = getStringOption_("database");
    String config = getStringOption_("config_file");
    int batch = getIntOption_("batch_size");
    String decoy_prefix = getStringOption_("decoy_prefix");

    QStringList arguments;
    arguments << config.toQString() 
              << "-f" << fasta_file.toQString() 
              << "-o" << output_folder.toQString() 
              << "--write-pin";
    if (batch >= 1) arguments << "--batch-size" << QString(batch);
    for (auto s : input_files) arguments << s.toQString();

    // Sage execution with the executable and the arguments StringList
    exit_code = runExternalProcess_(sage_executable.toQString(), arguments);
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    //-------------------------------------------------------------
    // writing IdXML output
    //-------------------------------------------------------------

    // read the sage output
    OPENMS_LOG_INFO << "Reading sage output..." << std::endl;
    StringList filenames;
    StringList extra_scores = {"ln(delta_next)", "ln(delta_best)", "matched_peaks", 
       "longest_b", "longest_y", "longest_y_pct",
       "ln(matched_intensity_pct)", "scored_candidates", "ln(-poisson)"};
    vector<PeptideIdentification> peptide_identifications = PercolatorInfile::load(
      output_folder + "/results.sage.pin",
      true,
      "ln(hyperscore)",
      extra_scores,
      filenames,
      decoy_prefix);

    // TODO: split / merge results and create idXMLs
    vector<ProteinIdentification> protein_identifications(1, ProteinIdentification());

    writeDebug_("write idXMLFile", 1);    
    
    protein_identifications[0].setPrimaryMSRunPath(filenames);    
    protein_identifications[0].setSearchEngineVersion(sage_version);

    // protein_identifications[0].getSearchParameters().enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(num_enzyme_termini[getStringOption_("num_enzyme_termini")]);
    protein_identifications[0].getSearchParameters().db = getStringOption_("database");
    
    // add extra scores for percolator rescoring
    vector<String> percolator_features = { "score" };
    for (auto s : extra_scores) percolator_features.push_back(s);
    protein_identifications[0].getSearchParameters().setMetaValue("extra_features",  ListUtils::concatenate(percolator_features, ","));

    // write all (!) parameters as metavalues to the search parameters
    if (!protein_identifications.empty())
    {
      DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_identifications[0].getSearchParameters(), this->getToolPrefix());
    }

    // if "reindex" parameter is set to true will perform reindexing
    if (auto ret = reindex_(protein_identifications, peptide_identifications); ret != EXECUTION_OK) return ret;

    IdXMLFile().store(output_file, protein_identifications, peptide_identifications);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPSageAdapter tool;
  return tool.main(argc, argv);
}

/// @endcond
