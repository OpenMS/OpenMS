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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Andreas Bertsch, Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <set>
#include <unordered_set>
#include <algorithm>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>

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
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ ProteinInterference \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
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

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    @todo Integrate Top x support, integrate parsimony approach from @ref OpenMS::PSProteinInference class
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
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    //TODO add function to merge based on replicates only. Needs additional exp. design file then.
    registerStringOption_("merge_runs", "<choice>", "no",
                          "If your idXML contains multiple runs, merge them beforehand?", false);
    setValidStrings_("merge_runs", ListUtils::create<String>("no,all"));

    registerStringOption_("annotate_indist_groups", "<choice>", "true",
        "If you want to annotate indistinguishable protein groups,"
        " either for reporting or for group based quant. later. Only works with a single ID run in the file.", false);
    setValidStrings_("annotate_indist_groups", ListUtils::create<String>("true,false"));

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
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    // load identifications
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    std::cout << "Loading input..." << std::endl;
    IdXMLFile().load(in, prot_ids, pep_ids);
    std::cout << "Loading input took " << sw.toString() << std::endl;
    sw.reset();

    bool merge_runs = getStringOption_("merge_runs") == "all";
    if (merge_runs)
    {
      std::cout << "Merging runs..." << std::endl;
      IDMergerAlgorithm merger;
      merger.insertRun(prot_ids, pep_ids);
      prot_ids.resize(1);
      merger.returnResultsAndClear(prot_ids[0], pep_ids);
      std::cout << "Merging runs took " << sw.toString() << std::endl;
      sw.reset();
    }

    std::cout << "Aggregating protein scores..." << std::endl;
    BasicProteinInferenceAlgorithm pi;
    pi.setParameters(getParam_().copy("Algorithm:", true));
    pi.run(pep_ids, prot_ids);
    std::cout << "Aggregating protein scores took " << sw.toString() << std::endl;
    sw.clear();

    bool annotate_indist_groups = getStringOption_("annotate_indist_groups") == "true";
    if (annotate_indist_groups)
    {
      if (prot_ids.size() > 1)
      {
        throw OpenMS::Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, prot_ids.size());
      }
      //TODO you could actually also do the aggregation/inference as well as the resolution on the Graph structure
      // but it is quite fast right now.
      IDBoostGraph ibg{prot_ids[0], pep_ids};
      ibg.buildGraph(0);
      sw.start();
      //TODO allow computation without splitting into components. Might be worthwhile in some cases
      std::cout << "Splitting into connected components..." << std::endl;
      ibg.computeConnectedComponents();
      std::cout << "Splitting into connected components took " << sw.toString() << std::endl;
      sw.clear();
      ibg.annotateIndistProteins(true);
    }

    std::cout << "Storing output..." << std::endl;
    sw.start();
    // write output
    IdXMLFile().store(out, prot_ids, pep_ids);
    std::cout << "Storing output took " << sw.toString() << std::endl;
    sw.stop();

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPProteinInference tool;
  return tool.main(argc, argv);
}

/// @endcond
