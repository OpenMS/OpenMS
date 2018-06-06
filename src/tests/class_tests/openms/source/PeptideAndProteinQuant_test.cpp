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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>

using namespace OpenMS;
using namespace std;


START_TEST(PeptideAndProteinQuant, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideAndProteinQuant* ptr = nullptr;
PeptideAndProteinQuant* nullPointer = nullptr;
START_SECTION((PeptideAndProteinQuant()))
  ptr = new PeptideAndProteinQuant();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeptideAndProteinQuant()))
  delete ptr;
END_SECTION

PeptideAndProteinQuant quantifier_features;
PeptideAndProteinQuant quantifier_consensus;
PeptideAndProteinQuant quantifier_identifications;
Param params;
params.setValue("include_all", "true");
quantifier_features.setParameters(params);
quantifier_consensus.setParameters(params);
quantifier_identifications.setParameters(params);

START_SECTION((void readQuantData(FeatureMap& features)))
{
  FeatureMap features;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.featureXML"), features);
  TEST_EQUAL(quantifier_features.getPeptideResults().empty(), true);
  quantifier_features.readQuantData(features);
  quantifier_features.quantifyPeptides();
  TEST_EQUAL(quantifier_features.getPeptideResults().empty(), false);
}
END_SECTION

START_SECTION((void readQuantData(ConsensusMap& consensus)))
{
  ConsensusMap consensus;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.consensusXML"), consensus);
  TEST_EQUAL(quantifier_consensus.getPeptideResults().empty(), true);
  quantifier_consensus.readQuantData(consensus);
  quantifier_consensus.quantifyPeptides();
  TEST_EQUAL(quantifier_consensus.getPeptideResults().empty(), false);
}
END_SECTION

START_SECTION((void readQuantData(vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides)))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.idXML"), proteins, peptides);
  TEST_EQUAL(quantifier_identifications.getPeptideResults().empty(), true);
  quantifier_identifications.readQuantData(proteins, peptides);
  quantifier_identifications.quantifyPeptides();
  TEST_EQUAL(quantifier_identifications.getPeptideResults().empty(), false);
}
END_SECTION

START_SECTION((void quantifyPeptides(const std::vector<PeptideIdentification>& peptides = std::vector<PeptideIdentification>())))
{
  NOT_TESTABLE // tested together with the "readQuantData" methods
}
END_SECTION

START_SECTION((void quantifyProteins(const ProteinIdentification& proteins = ProteinIdentification())))
{
  TEST_EQUAL(quantifier_features.getProteinResults().empty(), true);
  quantifier_features.quantifyProteins();
  TEST_EQUAL(quantifier_features.getProteinResults().empty(), false);

  TEST_EQUAL(quantifier_consensus.getProteinResults().empty(), true);
  quantifier_consensus.quantifyProteins();
  TEST_EQUAL(quantifier_consensus.getProteinResults().empty(), false);

  TEST_EQUAL(quantifier_identifications.getProteinResults().empty(), true);
  quantifier_identifications.quantifyProteins();
  TEST_EQUAL(quantifier_identifications.getProteinResults().empty(), false);
}
END_SECTION

START_SECTION((const Statistics& getStatistics()))
{
  PeptideAndProteinQuant::Statistics stats;

  stats = quantifier_features.getStatistics();
  TEST_EQUAL(stats.n_samples, 1);
  TEST_EQUAL(stats.quant_proteins, 2);
  TEST_EQUAL(stats.too_few_peptides, 1);
  TEST_EQUAL(stats.quant_peptides, 5);
  TEST_EQUAL(stats.total_peptides, 7);
  TEST_EQUAL(stats.quant_features, 7);
  TEST_EQUAL(stats.total_features, 8);
  TEST_EQUAL(stats.blank_features, 0);
  TEST_EQUAL(stats.ambig_features, 1);

  stats = quantifier_consensus.getStatistics();
  TEST_EQUAL(stats.n_samples, 3);
  TEST_EQUAL(stats.quant_proteins, 1);
  TEST_EQUAL(stats.too_few_peptides, 0);
  TEST_EQUAL(stats.quant_peptides, 4);
  TEST_EQUAL(stats.total_peptides, 4);
  TEST_EQUAL(stats.quant_features, 9);
  TEST_EQUAL(stats.total_features, 9);
  TEST_EQUAL(stats.blank_features, 0);
  TEST_EQUAL(stats.ambig_features, 0);

  stats = quantifier_identifications.getStatistics();
  TEST_EQUAL(stats.n_samples, 2);
  TEST_EQUAL(stats.quant_proteins, 10);
  TEST_EQUAL(stats.too_few_peptides, 10);
  TEST_EQUAL(stats.quant_peptides, 13);
  TEST_EQUAL(stats.total_peptides, 13);
  TEST_EQUAL(stats.quant_features, 18);
  TEST_EQUAL(stats.total_features, 18);
  TEST_EQUAL(stats.blank_features, 0);
  TEST_EQUAL(stats.ambig_features, 0);
}
END_SECTION

START_SECTION((const PeptideQuant& getPeptideResults()))
{
  PeptideAndProteinQuant::PeptideQuant pep_quant;
  PeptideAndProteinQuant::PeptideData pep_data;

  pep_quant = quantifier_features.getPeptideResults();
  TEST_EQUAL(pep_quant.size(), 7);
  pep_data = pep_quant[AASequence::fromString("AAAAA")];
  TEST_EQUAL(pep_data.abundances.size(), 1);
  TEST_EQUAL(pep_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 3333);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.id_count, 2);
  pep_data = pep_quant[AASequence::fromString("CCCCC")];
  TEST_EQUAL(pep_data.abundances.size(), 2);
  TEST_EQUAL(pep_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 7777);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.id_count, 2);
  pep_data = pep_quant[AASequence::fromString("EEEEE")];
  TEST_EQUAL(pep_data.abundances.size(), 1);
  TEST_EQUAL(pep_data.total_abundances.size(), 0);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.id_count, 1);
  pep_data = pep_quant[AASequence::fromString("GGGGG")];
  TEST_EQUAL(pep_data.abundances.size(), 1);
  TEST_EQUAL(pep_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 7777);
  TEST_EQUAL(pep_data.accessions.size(), 2);
  TEST_EQUAL(pep_data.id_count, 1);

  pep_quant = quantifier_consensus.getPeptideResults();
  TEST_EQUAL(pep_quant.size(), 4);
  pep_data = pep_quant[AASequence::fromString("AAA")];
  TEST_EQUAL(pep_data.abundances.size(), 1);
  TEST_EQUAL(pep_data.total_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 1000);
  TEST_REAL_SIMILAR(pep_data.total_abundances[2], 1000);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.id_count, 1);
  pep_data = pep_quant[AASequence::fromString("CCC")];
  TEST_EQUAL(pep_data.abundances.size(), 1);
  TEST_EQUAL(pep_data.total_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 200);
  TEST_REAL_SIMILAR(pep_data.total_abundances[1], 200);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.id_count, 1);
  pep_data = pep_quant[AASequence::fromString("EEE")];
  TEST_EQUAL(pep_data.abundances.size(), 1);
  TEST_EQUAL(pep_data.total_abundances.size(), 3);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 30);
  TEST_REAL_SIMILAR(pep_data.total_abundances[1], 30);
  TEST_REAL_SIMILAR(pep_data.total_abundances[2], 30);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.id_count, 1);
  pep_data = pep_quant[AASequence::fromString("GGG")];
  TEST_EQUAL(pep_data.abundances.size(), 1);
  TEST_EQUAL(pep_data.total_abundances.size(), 2);
  TEST_REAL_SIMILAR(pep_data.total_abundances[0], 4);
  TEST_REAL_SIMILAR(pep_data.total_abundances[1], 4);
  TEST_EQUAL(pep_data.accessions.size(), 1);
  TEST_EQUAL(pep_data.id_count, 1);
}
END_SECTION

START_SECTION((const ProteinQuant& getProteinResults()))
{
  PeptideAndProteinQuant::ProteinQuant prot_quant;
  PeptideAndProteinQuant::ProteinData prot_data;

  prot_quant = quantifier_features.getProteinResults();
  TEST_EQUAL(prot_quant.size(), 2);
  prot_data = prot_quant["Protein0"];
  TEST_EQUAL(prot_data.abundances.size(), 3);
  TEST_EQUAL(prot_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(prot_data.total_abundances[0], 4711);
  TEST_EQUAL(prot_data.id_count, 6);
  prot_data = prot_quant["Protein1"];
  TEST_EQUAL(prot_data.abundances.size(), 1);
  TEST_EQUAL(prot_data.total_abundances.size(), 1);
  TEST_REAL_SIMILAR(prot_data.total_abundances[0], 8888);
  TEST_EQUAL(prot_data.id_count, 2);

  prot_quant = quantifier_consensus.getProteinResults();
  TEST_EQUAL(prot_quant.size(), 1);
  prot_data = prot_quant["Protein"];
  TEST_EQUAL(prot_data.abundances.size(), 4);
  TEST_EQUAL(prot_data.total_abundances.size(), 3);
  TEST_REAL_SIMILAR(prot_data.total_abundances[0], 200);
  TEST_REAL_SIMILAR(prot_data.total_abundances[1], 30);
  TEST_REAL_SIMILAR(prot_data.total_abundances[2], 515);
  TEST_EQUAL(prot_data.id_count, 4);
}
END_SECTION

START_SECTION(([PeptideAndProteinQuant::PeptideData] PeptideData()))
{
  PeptideAndProteinQuant::PeptideData data;
  TEST_EQUAL(data.abundances.empty(), true);
  TEST_EQUAL(data.total_abundances.empty(), true);
  TEST_EQUAL(data.accessions.empty(), true);
  TEST_EQUAL(data.id_count, 0);
}
END_SECTION

START_SECTION(([PeptideAndProteinQuant::ProteinData] ProteinData()))
{
  PeptideAndProteinQuant::ProteinData data;
  TEST_EQUAL(data.abundances.empty(), true);
  TEST_EQUAL(data.total_abundances.empty(), true);
  TEST_EQUAL(data.id_count, 0);
}
END_SECTION

START_SECTION(([PeptideAndProteinQuant::Statistics] Statistics()))
{
  PeptideAndProteinQuant::Statistics stats;
  TEST_EQUAL(stats.n_samples, 0);
  TEST_EQUAL(stats.quant_proteins, 0);
  TEST_EQUAL(stats.too_few_peptides, 0);
  TEST_EQUAL(stats.quant_peptides, 0);
  TEST_EQUAL(stats.total_peptides, 0);
  TEST_EQUAL(stats.quant_features, 0);
  TEST_EQUAL(stats.total_features, 0);
  TEST_EQUAL(stats.blank_features, 0);
  TEST_EQUAL(stats.ambig_features, 0);
}
END_SECTION

// testing various averaging strategies
START_SECTION((const ProteinQuant& getProteinResults()))
{
  FeatureMap f;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ProteinQuantifier_input.featureXML"), f);
  
  PeptideAndProteinQuant quantifier;
  PeptideAndProteinQuant::ProteinQuant quant;
  PeptideAndProteinQuant::ProteinData protein;
  Param parameters;
  parameters.setValue("top", 0);
  
  parameters.setValue("average", "median");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 4711);

  parameters.setValue("average", "mean");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 5273.666666);

  parameters.setValue("average", "weighted_mean");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 5927.82624360028);

  parameters.setValue("average", "sum");
  quantifier.setParameters(parameters);
  quantifier.readQuantData(f);
  quantifier.quantifyPeptides();
  quantifier.quantifyProteins();
  quant = quantifier.getProteinResults();
  protein = quant["Protein0"];
  TEST_REAL_SIMILAR(protein.total_abundances[0], 15821);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



