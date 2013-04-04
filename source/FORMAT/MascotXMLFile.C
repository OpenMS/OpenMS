// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/SYSTEM/File.h>


using namespace xercesc;
using namespace std;

namespace OpenMS
{

  MascotXMLFile::MascotXMLFile() :
    Internal::XMLFile()
  {
  }

  void MascotXMLFile::load(const String & filename,
                           ProteinIdentification & protein_identification,
                           vector<PeptideIdentification> & id_data,
                           const RTMapping & rt_mapping)
  {
    map<String, vector<AASequence> > peptides;

    load(filename, protein_identification, id_data, peptides, rt_mapping);
  }

  void MascotXMLFile::load(const String & filename,
                           ProteinIdentification & protein_identification,
                           vector<PeptideIdentification> & id_data,
                           map<String, vector<AASequence> > & peptides,
                           const RTMapping & rt_mapping)
  {
    //clear
    protein_identification = ProteinIdentification();
    id_data.clear();

    Internal::MascotXMLHandler handler(protein_identification, id_data, filename, peptides, rt_mapping);
    parse_(filename, &handler);

    // Since the mascot xml can contain "peptides" without sequences the identifications
    // without any real peptide hit are removed
    {
      vector<PeptideIdentification> filtered_hits;
      filtered_hits.reserve(id_data.size());
      vector<PeptideIdentification>::iterator id_it = id_data.begin();

      while (id_it != id_data.end())
      {
        const vector<PeptideHit> & peptide_hits = id_it->getHits();
        if (peptide_hits.empty() || (peptide_hits.size() == 1 && peptide_hits[0].getSequence() == ""))
        {
          //std::cerr << "removing ID: " << std::distance(id_data.begin(), id_it) << "\n";
        }
        else
        {
          filtered_hits.push_back(*id_it);
        }
        ++id_it;
      }
      id_data.swap(filtered_hits);
    }

    // argh!
    // since Mascot xml 2.2 tends to repeat the first hit (yes it appears twice, we delete one of them)
    for (vector<PeptideIdentification>::iterator it = id_data.begin(); it != id_data.end(); ++it)
    {
      vector<PeptideHit> peptide_hits = it->getHits();
      // check if equal, except for rank
      if (peptide_hits.size() > 1 &&
          peptide_hits[0].getScore() == peptide_hits[1].getScore() &&
          peptide_hits[0].getSequence() == peptide_hits[1].getSequence() &&
          peptide_hits[0].getCharge() == peptide_hits[1].getCharge() /* &&
                  peptide_hits[0].getProteinAccessions() == peptide_hits[1].getProteinAccessions() &&
                    peptide_hits[0].getAABefore() == peptide_hits[1].getAABefore() &&
                    peptide_hits[0].getAAAfter() == peptide_hits[1].getAAAfter()*/)
      {
        // erase first hit
        peptide_hits.erase(peptide_hits.begin() + 1);
        it->setHits(peptide_hits);
      }
    }
  }

} // namespace OpenMS
