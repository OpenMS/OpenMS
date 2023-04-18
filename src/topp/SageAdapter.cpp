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
// $Maintainer: Chris Bielow $
// $Authors: Leon Bichmann, Timo Sachsenberg $
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

std::tuple<std::string, std::string, std::string> getVersionNumber_(const std::string& multi_line_input)
{
    std::regex version_regex("Version ([0-9]+)\\.([0-9]+)\\.([0-9]+)");

    std::sregex_iterator it(multi_line_input.begin(), multi_line_input.end(), version_regex);
    std::smatch match = *it;
    std::cout << "Found version " << match.str() << std::endl;      
        
    return make_tuple(it->str(1), it->str(2), it->str(3)); // major, minor, patch
}

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", { "mzML" } );

    registerOutputFile_("out", "<file>", "", "Path to the output file.", true, false);
    setValidFormats_("out", { "idXML" } );

    registerInputFile_("database", "<file>", "", "FASTA file", true, false, {"skipexists"});
    setValidFormats_("database", { "FASTA" } );
    registerInputFile_("sage_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      "sage", // this is the name on ALL platforms currently...
      "The Sage executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", true, false, {"is_executable"});

    registerInputFile_("config_file", "<file>", "", "Default Sage config file.", false, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("config_file", ListUtils::create<String>("json"));

    registerStringOption_("decoy_prefix", "<prefix>", "DECOY_", "Prefix on protein accession used to distinguish decoy from target proteins.");
    registerIntOption_("batch_size", "<int>", 0, "Number of files to load and search in parallel (default = # of CPUs/2)", false, false);
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
    String output_folder = File::path(output_file);
    String fasta_file = getStringOption_("database");
    String config = getStringOption_("config_file");
    int batch = getIntOption_("batch_size");

    QStringList arguments;
    arguments << config.toQString() << "-f" << fasta_file.toQString() << "-o" << output_folder.toQString() << "--write-pin";
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

    // read the sage pin and pout
    StringList filenames;
    vector<PeptideIdentification> peptide_identifications = PercolatorInfile::load("results.sage.pin", true, "ln(hyperscore)", filenames, decoy_prefix);

    // TODO: split / merge results and create idXMLs    
    vector<ProteinIdentification> protein_identifications;

    writeDebug_("write idXMLFile", 1);    
    
    protein_identifications[0].setPrimaryMSRunPath(filenames);    
    protein_identifications[0].setSearchEngineVersion(sage_version);

    // protein_identifications[0].getSearchParameters().enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(num_enzyme_termini[getStringOption_("num_enzyme_termini")]);
    protein_identifications[0].getSearchParameters().db = getStringOption_("database");

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
