// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Petra Gutenbrunner $
// $Authors: Petra Gutenbrunner $
// --------------------------------------------------------------------------

//! [Identification]

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // Create new protein identification object corresponding to a single search
  vector<ProteinIdentification> protein_ids;
  ProteinIdentification protein_id;
  
  protein_id.setIdentifier("Identifier");

  // Each ProteinIdentification object stores a vector of protein hits
  vector<ProteinHit> protein_hits;
  ProteinHit protein_hit = ProteinHit();
  protein_hit.setAccession("MyAccession");
  protein_hit.setSequence("PEPTIDEPEPTIDEPEPTIDEPEPTIDER");
  protein_hit.setScore(1.0);
  protein_hits.push_back(protein_hit);
  
  protein_id.setHits(protein_hits);    

  DateTime now = DateTime::now();
  String date_string = now.getDate();
  protein_id.setDateTime(now);
    
  // Example of possible search parameters
  ProteinIdentification::SearchParameters search_parameters;
  search_parameters.db = "database";
  search_parameters.charges = "+2";  
  protein_id.setSearchParameters(search_parameters);
 
  // Some seach engine meta data 
  protein_id.setSearchEngineVersion("v1.0.0");
  protein_id.setSearchEngine("SearchEngine");
  protein_id.setScoreType("HyperScore");
  
  protein_ids.push_back(protein_id);

  // Iterate over protein identifications and protein hits
  for (auto it = protein_ids.begin(); it != protein_ids.end(); ++it)
  {
    for (auto hit = it->getHits().begin(); hit < it->getHits().end(); ++hit)
    {
      cout << "Protein hit accession: " << hit->getAccession() << endl;
      cout << "Protein hit sequence: " << hit->getSequence() << endl;
      cout << "Protein hit score: " << hit->getScore() << endl;
    }
  }

  // Create new peptide identifications
  vector<PeptideIdentification> peptide_ids;
  PeptideIdentification peptide_id;

  peptide_id.setRT(1243.56);
  peptide_id.setMZ(440.0);
  peptide_id.setScoreType("ScoreType");
  peptide_id.setHigherScoreBetter(false);
  peptide_id.setIdentifier("Method");
  
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
  
  // add PeptideHit to PeptideIdentification
  peptide_id.setHits(peptide_hits);    
  
  // add PeptideIdentification
  peptide_ids.push_back(peptide_id);
  
  // We could now store the identification data in an idXML file  
  // IdXMLFile().store(out, protein_ids, peptide_ids);
  // And load it back with
  // IdXMLFile().load(id, protein_ids, peptide_ids);
   
  // Iterate over PeptideIdentification
  for (auto peptide_id = peptide_ids.begin(); peptide_id != peptide_ids.end(); 
    ++peptide_id)
  {
    // Peptide identification values
    cout << "Peptide ID m/z: " << peptide_id->getMZ() << endl;
    cout << "Peptide ID rt: " << peptide_id->getRT() << endl;
    cout << "Peptide ID score type: " << peptide_id->getScoreType() << endl;
    
    // PeptideHits
    for (auto hit = peptide_id->getHits().begin(); 
      hit < peptide_id->getHits().end(); 
      ++hit)
    {
      PeptideHit scored_hit = *hit;
      cout << "Peptide hit sequence: " << scored_hit.getSequence().toString() << endl;
      cout << "Peptide hit score: " << scored_hit.getScore() << endl;
    }
  }
  // ...
  return 0;
}


//! [Identification]
