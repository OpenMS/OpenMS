//! [InputRead_1]

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

//! [InputRead_1]
//! [InputRead_2]

// read the protein database
vector<FASTAFile::FASTAEntry> db;
FASTAFile().load(in, db);

// read the identification file (contains both protein as well as peptide identifications) 
vector<ProteinIdentification> protein_identifications;
vector<PeptideIdentification> peptide_identifications;
IdXMLFile().load(ids, protein_identifications, peptide_identifications);

//! [InputRead_2]
