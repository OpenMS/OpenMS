// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <vector>

using namespace OpenMS;
using namespace std;

class TOPPBayesianProteinInference :
public TOPPBase
{
public:
  TOPPBayesianProteinInference() :
  TOPPBase("BayesianProteinInference", "Runs a Bayesian protein inference.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file>", StringList(), "Input: identification results");
    setValidFormats_("in", ListUtils::create<String>("idXML,consensusXML,fasta"));
    registerInputFile_("exp_design", "<file>", "", "Input: experimental design");
    setValidFormats_("exp_design", ListUtils::create<String>("tsv"));
    //TODO make required of course
    registerOutputFile_("out", "<file>", "", "Output: identification results with scored/grouped proteins");
    setValidFormats_("out", ListUtils::create<String>("idXML,consensusXML"));
    registerFlag_("separate_runs", "Process multiple protein identification runs in the input separately,"
                                   " don't merge them. Merging results in loss of descriptive information"
                                   " of the single protein identification runs.", false);
  }

  ExitCodes main_(int, const char**) override
  {
    StringList files = getStringList_("in");
    IdXMLFile idXMLf;
    IDMergerAlgorithm merger{};
    for (String& file : files)
    {
      vector<ProteinIdentification> prots;
      vector<PeptideIdentification> peps;
      idXMLf.load(file,prots,peps);
      merger.insertRun(prots,peps);
    }
    vector<ProteinIdentification> mergedprots{1};
    vector<PeptideIdentification> mergedpeps;
    merger.returnResultsAndClear(mergedprots[0], mergedpeps);
    idXMLf.store(getStringOption_("out"),mergedprots,mergedpeps);
    BayesianProteinInferenceAlgorithm bpi1;
    bpi1.inferPosteriorProbabilities(mergedprots, mergedpeps);
    return ExitCodes::EXECUTION_OK;

    /* That was test code for the old merger
    ProteinIdentification pi1{};
    ProteinIdentification pi2{};
    ProteinIdentification pi3{};
    ProteinIdentification pi4{};
    pi1.setPrimaryMSRunPath({"a.mzml", "b.mzml", "c.mzml"});
    pi2.setPrimaryMSRunPath({"d.mzml", "c.mzml", "e.mzml"});
    pi3.setPrimaryMSRunPath({"g.mzml", "f.mzml", "e.mzml"});
    pi4.setPrimaryMSRunPath({"h.mzml", "e.mzml"});
    vector<ProteinIdentification> pis{pi1,pi2,pi3,pi4};
    vector<PeptideIdentification> peps1{};
    IDMergerAlgorithm merger1;
    vector<vector<PeptideIdentification>> newpeps1{};
    merger1.mergeIDRunsAndSplitPeptides(pis, peps1, {{0,0},{1,0},{2,1},{3,1}}, newpeps1);

    return ExitCodes::EXECUTION_OK;


    //TODO remove until ExitCode and allow a (merged) ConsensusXML as input
    vector<PeptideIdentification> peps;
    vector<ProteinIdentification> prots;
    ConsensusMap cmap;
    ConsensusXMLFile cXML;
    cXML.load(getStringOption_("in"), cmap);
    IDMergerAlgorithm merger;
    ExperimentalDesignFile expFile;
    ExperimentalDesign expDesign = expFile.load(getStringOption_("exp_design"), false);
    merger.mergeProteinsAcrossFractionsAndReplicates(cmap, expDesign);
    cXML.store(getStringOption_("out"), cmap);

    return ExitCodes::EXECUTION_OK;
     */

    // Some thoughts about how to leverage info from different runs.
    //Fractions: Always merge (not much to leverage, maybe agreement at borders)
    // - Think about only allowing one/the best PSM per peptidoform across fractions
    //Replicates: Use matching ID and quant, also a always merge
    //Samples: In theory they could yield different proteins/pep-protein-associations
    // 3 options:
    // - don't merge: -> don't leverage peptide quant profiles (or use them repeatedly -> same as second opt.?)
    // - merge and assume same proteins: -> You can use the current graph and weigh the associations
    //   based on deviation in profiles
    // - merge and don't assume same proteins: -> We need an extended graph, that has multiple versions
    //   of the proteins for every sample

    vector<PeptideIdentification> peps;
    vector<ProteinIdentification> prots;
    IdXMLFile idXML;
    idXML.load(getStringOption_("in"), prots, peps);
    //TODO filter unmatched proteins and peptides before!
    //TODO check t+d annotations

    BayesianProteinInferenceAlgorithm bpi;
    bpi.inferPosteriorProbabilities(prots, peps);
    idXML.store(getStringOption_("out"),prots, peps);

    return ExitCodes::EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPBayesianProteinInference tool;

  return tool.main(argc, argv);
}
