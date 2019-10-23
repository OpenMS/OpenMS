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

#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#ifdef _OPENMP
  #include <omp.h>
#endif


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_SimpleSearchEngine SimpleSearchEngine

    @brief Identifies peptides in MS/MS spectra.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ SimpleSearchEngine \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em This search engine is mainly for educational/benchmarking/prototyping use cases.
    It lacks behind in speed and/or quality of results when compared to state-of-the-art search engines.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_SimpleSearchEngine.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_SimpleSearchEngine.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class SimpleSearchEngine :
    public TOPPBase
{
  public:
    SimpleSearchEngine() :
      TOPPBase("SimpleSearchEngine",
        "Annotates MS/MS spectra using SimpleSearchEngine.",
        false)
    {
    }

  protected:
    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "", "input file ");
      setValidFormats_("in", ListUtils::create<String>("mzML"));

      registerInputFile_("database", "<file>", "", "input file ");
      setValidFormats_("database", ListUtils::create<String>("fasta"));

      registerOutputFile_("out", "<file>", "", "output file ");
      setValidFormats_("out", ListUtils::create<String>("idXML"));

      // put search algorithm parameters at Search: subtree of parameters
      Param search_algo_params_with_subsection;
      search_algo_params_with_subsection.insert("Search:", SimpleSearchEngineAlgorithm().getDefaults());
      registerFullParam_(search_algo_params_with_subsection);
    }

    ExitCodes main_(int, const char**) override
    {
      String in = getStringOption_("in");
      String database = getStringOption_("database");
      String out = getStringOption_("out");

      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);

      vector<ProteinIdentification> protein_ids;
      vector<PeptideIdentification> peptide_ids;

      SimpleSearchEngineAlgorithm sse;
      sse.setParameters(getParam_().copy("Search:", true));
      //TODO ??? Why not use the TOPPBase ExitCodes?
      // same for OpenPepXL etc. Otherwise please write a proper mapping.
      SimpleSearchEngineAlgorithm::ExitCodes e = sse.search(in, database, protein_ids, peptide_ids);
      if (e != SimpleSearchEngineAlgorithm::ExitCodes::EXECUTION_OK)
      {
        return TOPPBase::ExitCodes::INTERNAL_ERROR;
      }

      // MS path already set in algorithm. Overwrite here so we get something testable
      if (getFlag_("test"))
      {
        // if test mode set, add file without path so we can compare it
        protein_ids[0].setPrimaryMSRunPath({"file://" + File::basename(in)});
      }

      IdXMLFile().store(out, protein_ids, peptide_ids);

      return EXECUTION_OK;
    }

};

int main(int argc, const char** argv)
{
  SimpleSearchEngine tool;
  return tool.main(argc, argv);
}
