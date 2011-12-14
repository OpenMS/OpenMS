// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>

using namespace OpenMS;
using namespace std;


START_TEST(PeptideAndProteinQuant, "$Id PeptideAndProteinQuant_test.C 7941 2011-02-11 22:31:39Z hendrikweisser $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideAndProteinQuant* ptr = 0;
PeptideAndProteinQuant* nullPointer = 0;
START_SECTION((PeptideAndProteinQuant()))
	ptr = new PeptideAndProteinQuant();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeptideAndProteinQuant()))
	delete ptr;
END_SECTION

PeptideAndProteinQuant quantifier_features;
PeptideAndProteinQuant quantifier_consensus;
Param params;
params.setValue("include_all", "true");
quantifier_features.setParameters(params);
quantifier_consensus.setParameters(params);

START_SECTION((void quantifyPeptides(FeatureMap<>& features)))
{
	FeatureMap<> features;
	FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("../TOPP/ProteinQuantifier_input.featureXML"), features);
	TEST_EQUAL(quantifier_features.getPeptideResults().empty(), true);
	quantifier_features.quantifyPeptides(features);
	TEST_EQUAL(quantifier_features.getPeptideResults().empty(), false);
}
END_SECTION

START_SECTION((void quantifyPeptides(ConsensusMap& consensus)))
{
	ConsensusMap consensus;
	ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("../TOPP/ProteinQuantifier_input.consensusXML"), consensus);
	TEST_EQUAL(quantifier_consensus.getPeptideResults().empty(), true);
	quantifier_consensus.quantifyPeptides(consensus);
	TEST_EQUAL(quantifier_consensus.getPeptideResults().empty(), false);
}
END_SECTION

START_SECTION((void quantifyProteins(const ProteinIdentification& proteins=ProteinIdentification())))
{
	TEST_EQUAL(quantifier_features.getProteinResults().empty(), true);
	quantifier_features.quantifyProteins();
	TEST_EQUAL(quantifier_features.getProteinResults().empty(), false);

	TEST_EQUAL(quantifier_consensus.getProteinResults().empty(), true);
	quantifier_consensus.quantifyProteins();
	TEST_EQUAL(quantifier_consensus.getProteinResults().empty(), false);
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
}
END_SECTION

START_SECTION((const PeptideQuant& getPeptideResults()))
{
	PeptideAndProteinQuant::PeptideQuant pep_quant;
	PeptideAndProteinQuant::PeptideData pep_data;

	pep_quant = quantifier_features.getPeptideResults();
	TEST_EQUAL(pep_quant.size(), 7);
	pep_data = pep_quant[AASequence("AAAAA")];
	TEST_EQUAL(pep_data.abundances.size(), 1);
	TEST_EQUAL(pep_data.total_abundances.size(), 1);
	TEST_REAL_SIMILAR(pep_data.total_abundances[0], 3333);
	TEST_EQUAL(pep_data.accessions.size(), 1);
	TEST_EQUAL(pep_data.id_count, 2);
	pep_data = pep_quant[AASequence("CCCCC")];
	TEST_EQUAL(pep_data.abundances.size(), 2);
	TEST_EQUAL(pep_data.total_abundances.size(), 1);
	TEST_REAL_SIMILAR(pep_data.total_abundances[0], 7777);
	TEST_EQUAL(pep_data.accessions.size(), 1);
	TEST_EQUAL(pep_data.id_count, 2);
	pep_data = pep_quant[AASequence("EEEEE")];
	TEST_EQUAL(pep_data.abundances.size(), 1);
	TEST_EQUAL(pep_data.total_abundances.size(), 0);
	TEST_EQUAL(pep_data.accessions.size(), 1);
	TEST_EQUAL(pep_data.id_count, 1);
	pep_data = pep_quant[AASequence("GGGGG")];
	TEST_EQUAL(pep_data.abundances.size(), 1);
	TEST_EQUAL(pep_data.total_abundances.size(), 1);
	TEST_REAL_SIMILAR(pep_data.total_abundances[0], 7777);
	TEST_EQUAL(pep_data.accessions.size(), 2);
	TEST_EQUAL(pep_data.id_count, 1);

	pep_quant = quantifier_consensus.getPeptideResults();
	TEST_EQUAL(pep_quant.size(), 4);
	pep_data = pep_quant[AASequence("AAA")];
	TEST_EQUAL(pep_data.abundances.size(), 1);
	TEST_EQUAL(pep_data.total_abundances.size(), 2);
	TEST_REAL_SIMILAR(pep_data.total_abundances[0], 1000);
	TEST_REAL_SIMILAR(pep_data.total_abundances[2], 1000);
	TEST_EQUAL(pep_data.accessions.size(), 1);
	TEST_EQUAL(pep_data.id_count, 1);
	pep_data = pep_quant[AASequence("CCC")];
	TEST_EQUAL(pep_data.abundances.size(), 1);
	TEST_EQUAL(pep_data.total_abundances.size(), 2);
	TEST_REAL_SIMILAR(pep_data.total_abundances[0], 200);
	TEST_REAL_SIMILAR(pep_data.total_abundances[1], 200);
	TEST_EQUAL(pep_data.accessions.size(), 1);
	TEST_EQUAL(pep_data.id_count, 1);
	pep_data = pep_quant[AASequence("EEE")];
	TEST_EQUAL(pep_data.abundances.size(), 1);
	TEST_EQUAL(pep_data.total_abundances.size(), 3);
	TEST_REAL_SIMILAR(pep_data.total_abundances[0], 30);
	TEST_REAL_SIMILAR(pep_data.total_abundances[1], 30);
	TEST_REAL_SIMILAR(pep_data.total_abundances[2], 30);
	TEST_EQUAL(pep_data.accessions.size(), 1);
	TEST_EQUAL(pep_data.id_count, 1);
	pep_data = pep_quant[AASequence("GGG")];
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
	TEST_REAL_SIMILAR(prot_data.total_abundances[0], 5555);
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



