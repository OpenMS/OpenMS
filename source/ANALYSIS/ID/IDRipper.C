// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Immanuel Luhn $
// $Authors: Immanuel Luhn $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDRipper.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <QDir>

using std::vector;
using std::map;
using std::pair;

//using namespace std;

namespace OpenMS
{

  IDRipper::IDRipper()
    : DefaultParamHandler("IDRipper"){}


  IDRipper::IDRipper(const IDRipper& cp)
    : DefaultParamHandler(cp){}

  IDRipper::~IDRipper(){}

  IDRipper& IDRipper::operator= (const IDRipper& rhs)
  {
    if (this == &rhs) return *this;

    DefaultParamHandler::operator= (rhs);
    updateMembers_();

    return *this;
  }

  void IDRipper::rip(map<String, pair< vector<ProteinIdentification>, vector<PeptideIdentification> > >& ripped, vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides)
  {
    // Collect all protein hits
    vector<ProteinHit> all_protein_hits;
    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin(); prot_it != proteins.end(); ++prot_it)
    {
      // remove file origin
      prot_it->removeMetaValue("file_origin");
      vector<ProteinHit>& protein_hits  = prot_it->getHits();
      all_protein_hits.insert(all_protein_hits.end(), protein_hits.begin(), protein_hits.end());
    }

    //store protein and peptides identifications for each file origin

    for ( vector<PeptideIdentification>::iterator pep_it = peptides.begin(); pep_it != peptides.end(); ++pep_it )
    {
      // try to get file_origin, if not present ignore peptide identification
      const String& file_origin = pep_it->getMetaValue("file_origin").toString();
      // QFileInfo fi("/tmp/archive.tar.gz");
      // QString name = fi.fileName(); --> name = "archive.tar.gz"
      const String file_ = QFileInfo(file_origin.toQString()).fileName().toStdString();

      //remove file origin
      pep_it->removeMetaValue("file_origin");

      //TODO LOG that file_origin was not as expected
      if ( file_.empty() ) continue;

      // try to get peptide hits for peptide identification
      const vector<PeptideHit>& peptide_hits = pep_it->getHits();
      if ( peptide_hits.empty() ) continue;

      // collect all protein accesions that are stored in the peptide hits
      vector<String> protein_accessions;
      getProteinAccessions_(protein_accessions, peptide_hits);

      // returns all protein hits that are associated with the given peptide hits
      vector<ProteinHit> protein2accessions;
      getProteinHits_(protein2accessions, all_protein_hits, protein_accessions);

      // search for the protein identification of the peptide identification
      ProteinIdentification prot_ident;
      getProteinIdentification_(prot_ident, *pep_it, proteins);
      // TODO catch case that ProteinIdentification prot_ident is not found in the for-loop


      map<String, pair< vector<ProteinIdentification>, vector<PeptideIdentification> > >::iterator it = ripped.find(file_);
      // If file_origin already exists
      if ( it != ripped.end() )
      {
        vector<ProteinIdentification>& prot_tmp = it->second.first;
        bool flag = true;
        //what to do if there is one then more protein identification, can this occur at all?
        for ( vector<ProteinIdentification>::iterator it2 = prot_tmp.begin(); it2 != prot_tmp.end(); ++it2)
        {
          // ProteinIdentification is already there, just add protein hits
          if ( prot_ident.getIdentifier().compare(it2->getIdentifier()) == 0 )
          {
            for ( vector<ProteinHit>::const_iterator prot_it = protein2accessions.begin(); prot_it != protein2accessions.end(); ++prot_it)
            {
              it2->insertHit(*prot_it);
            }
            flag = false;
            break;
          }
        }
        // if it was not found
        if ( flag )
        {
          prot_ident.setHits(protein2accessions);
          prot_tmp.push_back(prot_ident);
        }
        vector<PeptideIdentification>& pep_tmp = it->second.second;
        pep_tmp.push_back(*pep_it);
      }
      else // otherwise create new entry for file_origin
      {
        // create protein identification, TODO parameters
        vector<ProteinIdentification> protein_idents;
        // only use the protein hits that are needed for the peptide identification
        prot_ident.setHits(protein2accessions);
        protein_idents.push_back(prot_ident);

        //create new peptide identification
        vector<PeptideIdentification> peptide_idents;
        peptide_idents.push_back(*pep_it);

        //create and insert new map entry
        ripped.insert(make_pair(file_, make_pair(protein_idents, peptide_idents)));
      }
    }
  }

  void IDRipper::getProteinHits_( vector<ProteinHit>& result, const vector<ProteinHit>& protein_hits, const vector<String>& protein_accessions)
  {
    for ( vector<String>::const_iterator it = protein_accessions.begin(); it < protein_accessions.end(); ++it)
    {
      for ( vector<ProteinHit>::const_iterator prot_it = protein_hits.begin(); prot_it != protein_hits.end(); ++prot_it)
      {
        if ( prot_it->getAccession().compare(*it) == 0 )
        {
          result.push_back(*prot_it);
        }
      }
    }
  }

  void IDRipper::getProteinAccessions_(vector<String>& result, const vector<PeptideHit>& peptide_hits)
  {
    for ( vector<PeptideHit>::const_iterator it = peptide_hits.begin(); it != peptide_hits.end(); ++it)
    {
      const vector<String>& protein_accessions_tmp = it->getProteinAccessions();
      result.insert(result.end(), protein_accessions_tmp.begin(), protein_accessions_tmp.end());
    }
  }

  void IDRipper::getProteinIdentification_(ProteinIdentification& result, PeptideIdentification pep_ident, std::vector<ProteinIdentification>& prot_idents)
  {
    const String& identifier = pep_ident.getIdentifier();

    for ( vector<ProteinIdentification>::iterator prot_it = prot_idents.begin(); prot_it != prot_idents.end(); ++prot_it )
    {
      if ( identifier.compare(prot_it->getIdentifier() ) == 0)
      {
        result = *prot_it;
        break;
      }
    }
  }

} // namespace OpenMS

