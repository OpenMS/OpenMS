// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Petra Gutenbrunner $
// $Authors: Petra Gutenbrunner $
// --------------------------------------------------------------------------

//! [doxygen_snippet_Identification]

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // Create new protein identification object corresponding to a single search

  // Each ProteinIdentification object stores a vector of protein hits
  vector<ProteinHit> protein_hits;
  ProteinHit protein_hit;
  protein_hit.setAccession("MyAccession");
  protein_hit.setSequence("PEPTIDEPEPTIDEPEPTIDEPEPTIDER");
  protein_hit.setScore(1.0);
  protein_hits.push_back(protein_hit);

  ProteinIdentification protein_id;
  protein_id.setIdentifier("Identifier");
  protein_id.setHits(protein_hits);

  DateTime now = DateTime::now();
  String date_string = now.getDate();
  protein_id.setDateTime(now);

  // Example of possible search parameters
  ProteinIdentification::SearchParameters search_parameters;
  search_parameters.db = "database";
  search_parameters.charges = "+2";
  protein_id.setSearchParameters(search_parameters);

  // Some search engine meta data
  protein_id.setSearchEngineVersion("v1.0.0");
  protein_id.setSearchEngine("SearchEngine");
  protein_id.setScoreType("HyperScore");

  vector<ProteinIdentification> protein_ids;
  protein_ids.push_back(protein_id);

  // Iterate over protein identifications and protein hits
  for (const auto& prot : protein_ids)
  {
    for (const auto& hit : prot.getHits())
    {
      cout << "Protein hit accession: " << hit.getAccession() << '\n';
      cout << "Protein hit sequence: " << hit.getSequence() << '\n';
      cout << "Protein hit score: " << hit.getScore() << '\n';
    }
  }

  // Create new peptide identifications
  vector<PeptideIdentification> peptide_ids;
  PeptideIdentification peptide_id;

  peptide_id.setRT(1243.56);
  peptide_id.setMZ(440.0);
  peptide_id.setScoreType("ScoreType");
  peptide_id.setHigherScoreBetter(false);
  peptide_id.setIdentifier("Identifier");

  // define additional meta value for the peptide identification
  peptide_id.setMetaValue("AdditionalMetaValue", "Value");

  // add PeptideHit to a PeptideIdentification
  vector<PeptideHit> peptide_hits;
  PeptideHit peptide_hit;
  peptide_hit.setScore(1.0);
  peptide_hit.setRank(1);
  peptide_hit.setCharge(2);
  peptide_hit.setSequence(AASequence::fromString("DLQM(Oxidation)TQSPSSLSVSVGDR"));
  peptide_hits.push_back(peptide_hit);

  // add second best PeptideHit to the PeptideIdentification
  peptide_hit.setScore(1.5);
  peptide_hit.setRank(2);
  peptide_hit.setCharge(2);
  peptide_hit.setSequence(AASequence::fromString("QLDM(Oxidation)TQSPSSLSVSVGDR"));
  peptide_hits.push_back(peptide_hit);

  // add PeptideHit to PeptideIdentification
  peptide_id.setHits(peptide_hits);

  // add PeptideIdentification
  peptide_ids.push_back(peptide_id);

  // We could now store the identification data in an idXML file
  // FileHandler().storeIdentifications(outfile, protein_ids, peptide_ids);
  // And load it back with
  // FileHandler().loadIdentifications(outfile, protein_ids, peptide_ids);

  // Iterate over PeptideIdentification
  for (const auto& peptide_id : peptide_ids)
  {
    // Peptide identification values
    cout << "Peptide ID m/z: " << peptide_id.getMZ() << '\n';
    cout << "Peptide ID rt: " << peptide_id.getRT() << '\n';
    cout << "Peptide ID score type: " << peptide_id.getScoreType() << '\n';

    // PeptideHits
    for (const auto& scored_hit : peptide_id.getHits())
    {
      cout << " - Peptide hit rank: " << scored_hit.getRank() << '\n';
      cout << " - Peptide hit sequence: " << scored_hit.getSequence().toString() << '\n';
      cout << " - Peptide hit score: " << scored_hit.getScore() << '\n';
    }
  }
}


//! [doxygen_snippet_Identification]
