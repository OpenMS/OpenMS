// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Immanuel Luhn $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDRipper.h>

#include <QDir>

using std::vector;
using std::map;
using std::pair;

//using namespace std;

namespace OpenMS
{

  IDRipper::IDRipper() :
    DefaultParamHandler("IDRipper")
  {
  }

  IDRipper::IDRipper(const IDRipper& cp) :
    DefaultParamHandler(cp)
  {
  }

  IDRipper::~IDRipper()
  {
  }

  IDRipper& IDRipper::operator=(const IDRipper& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    updateMembers_();

    return *this;
  }

  void IDRipper::rip(map<String, pair<vector<ProteinIdentification>, vector<PeptideIdentification> > >& ripped, vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides)
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

    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      // try to get file_origin, if not present ignore peptide identification
      const String& file_origin = pep_it->getMetaValue("file_origin").toString();
      // QFileInfo fi("/tmp/archive.tar.gz");
      // QString name = fi.fileName(); --> name = "archive.tar.gz"
      const String file_ = QFileInfo(file_origin.toQString()).fileName().toStdString();

      //remove file origin
      pep_it->removeMetaValue("file_origin");

      //TODO LOG that file_origin was not as expected
      if (file_.empty())
        continue;

      // try to get peptide hits for peptide identification
      const vector<PeptideHit>& peptide_hits = pep_it->getHits();
      if (peptide_hits.empty())
        continue;

      // collect all protein accessions that are stored in the peptide hits
      vector<String> protein_accessions;
      getProteinAccessions_(protein_accessions, peptide_hits);

      // returns all protein hits that are associated with the given peptide hits
      vector<ProteinHit> protein2accessions;
      getProteinHits_(protein2accessions, all_protein_hits, protein_accessions);

      // search for the protein identification of the peptide identification
      ProteinIdentification prot_ident;
      getProteinIdentification_(prot_ident, *pep_it, proteins);
      // TODO catch case that ProteinIdentification prot_ident is not found in the for-loop


      map<String, pair<vector<ProteinIdentification>, vector<PeptideIdentification> > >::iterator it = ripped.find(file_);
      // If file_origin already exists
      if (it != ripped.end())
      {
        vector<ProteinIdentification>& prot_tmp = it->second.first;
        bool flag = true;
        //what to do if there is one then more protein identification, can this occur at all?
        for (vector<ProteinIdentification>::iterator it2 = prot_tmp.begin(); it2 != prot_tmp.end(); ++it2)
        {
          // ProteinIdentification is already there, just add protein hits
          if (prot_ident.getIdentifier().compare(it2->getIdentifier()) == 0)
          {
            for (vector<ProteinHit>::const_iterator prot_it = protein2accessions.begin(); prot_it != protein2accessions.end(); ++prot_it)
            {
              it2->insertHit(*prot_it);
            }
            flag = false;
            break;
          }
        }
        // if it was not found
        if (flag)
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

  void IDRipper::getProteinHits_(vector<ProteinHit>& result, const vector<ProteinHit>& protein_hits, const vector<String>& protein_accessions)
  {
    for (vector<String>::const_iterator it = protein_accessions.begin(); it < protein_accessions.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator prot_it = protein_hits.begin(); prot_it != protein_hits.end(); ++prot_it)
      {
        if (prot_it->getAccession().compare(*it) == 0)
        {
          result.push_back(*prot_it);
        }
      }
    }
  }

  void IDRipper::getProteinAccessions_(vector<String>& result, const vector<PeptideHit>& peptide_hits)
  {
    for (vector<PeptideHit>::const_iterator it = peptide_hits.begin(); it != peptide_hits.end(); ++it)
    {
      std::set<String> protein_accessions = it->extractProteinAccessionsSet();
      result.insert(result.end(), protein_accessions.begin(), protein_accessions.end());
    }
  }

  void IDRipper::getProteinIdentification_(ProteinIdentification& result, PeptideIdentification pep_ident, std::vector<ProteinIdentification>& prot_idents)
  {
    const String& identifier = pep_ident.getIdentifier();

    for (vector<ProteinIdentification>::iterator prot_it = prot_idents.begin(); prot_it != prot_idents.end(); ++prot_it)
    {
      if (identifier.compare(prot_it->getIdentifier()) == 0)
      {
        result = *prot_it;
        break;
      }
    }
  }

} // namespace OpenMS
