// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // Example 1: Writing and reading compressed IdXML files
  {
    IdXMLFile id_file;
    vector<ProteinIdentification> protein_ids;
    PeptideIdentificationList peptide_ids;
    String document_id = "example_id";
    
    // ... populate protein_ids and peptide_ids with data ...
    
    // Write uncompressed file
    id_file.store("output.idXML", protein_ids, peptide_ids, document_id);
    
    // Write gzip-compressed file (automatically detected by .gz extension)
    id_file.store("output.idXML.gz", protein_ids, peptide_ids, document_id);
    
    // Write bzip2-compressed file (automatically detected by .bz2 extension)
    id_file.store("output.idXML.bz2", protein_ids, peptide_ids, document_id);
    
    // Reading works the same way - compression is automatically detected
    vector<ProteinIdentification> loaded_protein_ids;
    PeptideIdentificationList loaded_peptide_ids;
    String loaded_document_id;
    
    // Load from gzip-compressed file
    id_file.load("output.idXML.gz", loaded_protein_ids, loaded_peptide_ids, loaded_document_id);
    
    cout << "Successfully wrote and read compressed IdXML files!" << endl;
  }
  
  // Example 2: Writing and reading compressed ConsensusXML files
  {
    ConsensusXMLFile consensus_file;
    ConsensusMap consensus_map;
    
    // ... populate consensus_map with data ...
    
    // Write gzip-compressed consensusXML
    consensus_file.store("consensus.consensusXML.gz", consensus_map);
    
    // Load from compressed file
    ConsensusMap loaded_map;
    consensus_file.load("consensus.consensusXML.gz", loaded_map);
    
    cout << "Successfully wrote and read compressed ConsensusXML files!" << endl;
  }
  
  return 0;
}
