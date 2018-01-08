//! [functionality_1]

// store protein accessions 
set<String> id_accessions;
// loop over all identified spectra
for (Size i = 0; i != peptide_identifications.size(); ++i)
{
  const PeptideIdentification& id = peptide_identifications[i];
  const vector<PeptideHit>& hits = id.getHits();
  // loop over every PSM of the current spectrum
  for (Size k = 0; k != hits.size(); ++k) // for every PSM
  {
    const vector<PeptideEvidence>& evidences = hits[k].getPeptideEvidences();
    // loop over every sequence to protein mapping
    for (Size m = 0; m != evidences.size(); ++m)
    {
      const String& id_accession = evidences[m].getProteinAccession();
      id_accessions.insert(id_accession); // add accession to set
    }
  }
}

//! [functionality_1]
//! [functionality_2]

// add method functionality 
vector<FASTAFile::FASTAEntry> db_new;
 
for (Size i = 0; i != db.size() ; ++i)
{
   const String& fasta_accession = db[i].identifier;
   const bool found = id_accessions.find(fasta_accession) != id_accessions.end();
   if ( (found && whitelist) || (!found && !whitelist) )
   {
      db_new.push_back(db[i]);
   }
}

//! [functionality_2]

//! [output]
FASTAFile().store(out, db_new);
//! [output]

