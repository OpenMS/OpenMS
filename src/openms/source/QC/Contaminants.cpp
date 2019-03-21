// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow $
// $Authors: Dominik Schmitz, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/QC/Contaminants.h>
#include <include/OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <include/OpenMS/METADATA/ProteinIdentification.h>
//#include <include/OpenMS/CHEMISTRY/DigestionEnzymeDB.h>
#include <include/OpenMS/CHEMISTRY/ProteaseDB.h>

using namespace OpenMS;
using namespace std;



void Contaminants::compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants)
{
  if (digested_db_.empty())
  {
    vector<ProteinIdentification> protein_identifications;

    vector<PeptideIdentification> identifications;
    PeptideIdentification peptide_identification;
    DateTime date_time = DateTime::now();
    String date_time_string = date_time.get();
    peptide_identification.setIdentifier("In-silico_digestion" + date_time_string);

    ProteinIdentification protein_identification;

    protein_identifications.push_back(ProteinIdentification());

    ProteinIdentification::SearchParameters search_parameters;
    String enzyme = features.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme.getName();
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    UInt missed_cleavages (features.getProteinIdentifications()[0].getSearchParameters().missed_cleavages);
    digestor.setMissedCleavages(missed_cleavages);
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme);

    PeptideHit temp_peptide_hit;
    PeptideEvidence temp_pe;

    protein_identifications[0].setSearchParameters(search_parameters);
    protein_identifications[0].setDateTime(date_time);
    protein_identifications[0].setSearchEngine("In-silico digestion");
    protein_identifications[0].setIdentifier("In-silico_digestion" + date_time_string);

    Size dropped_by_length(0); // stats for removing candidates
    Size fasta_out_count(0);

    for (const FASTAFile::FASTAEntry& fe : contaminants)
    {
      vector<AASequence> current_digest;
      if (enzyme == "none")
      {
        current_digest.push_back(AASequence::fromString(fe.sequence));
      }
      else
      {
        dropped_by_length += digestor.digest(AASequence::fromString(fe.sequence), current_digest);
      }

      String id = fe.identifier;
      for (auto const& s : current_digest)
      {
          ++fasta_out_count;
          id = fe.identifier + "_" + String(fasta_out_count); // used switch case "BOTH" 
          digested_db_.push_back(FASTAFile::FASTAEntry(id, fe.description , s.toString())); // FASTA_desc set to true because description is always wanted
          //digested_sequences_.push_back(fe.getSequence())
      }
      //digested_sequences_.push_back(digested_db_.getSequence());
    }
  }
    Int64 total;
    Int64 cont;
    Int64 sum_total;
    Int64 sum_cont;

    if (contaminants.size() == 0)
    {
      std::cerr << "FASTAFile is empty";
    }
    else
    {
      for (auto f : features)
      {
        for (auto pep_id : f.getPeptideIdentifications())
        {
          for (auto pep_hit : pep_id.getHits())
          {
            if (std::find(digested_sequences_.begin(), digested_sequences_.end(), (pep_hit.getSequence())) != digested_sequences_.end())
            {
              int x = 1;
              ++total;
              ++cont;
              sum_total += f.getIntensity();
              sum_cont += f.getIntensity();
              pep_hit.setMetaValue("is_contaminant", x);

            }
            else
            {
              int x = 0;
              ++total;
              sum_total += f.getIntensity();
              pep_hit.setMetaValue("is_contaminant", x);
            }
          }
        }
      }
      results_.push_back(std::make_tuple(cont / total, sum_cont / sum_total));
    }
}

const std::vector<std::tuple<double, double>>& Contaminants::getResults()
{
  return results_;
}
