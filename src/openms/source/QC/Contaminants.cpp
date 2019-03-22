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
//#include <include/OpenMS/CHEMISTRY/ProteaseDB.h>

using namespace OpenMS;
using namespace std;



void Contaminants::compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants)
{
  if (digested_db_.empty())
  {
    String enzyme = features.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme.getName();
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    UInt missed_cleavages (features.getProteinIdentifications()[0].getSearchParameters().missed_cleavages);
    digestor.setMissedCleavages(missed_cleavages);

    Size dropped_by_length(0); // stats for removing candidates


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
//      DigestedProtein temp;
//      temp.id = fe.identifier;
      for (auto const& s : current_digest)
      {
         digested_db_.insert(s.toUnmodifiedString());
          //temp.peptides.push_back(s.toString()); // FASTA_desc set to true because description is always wanted
      }
      //digested_db_.push_back(temp);
    }
  }
    Int64 total;
    Int64 cont;
    double sum_total;
    double sum_cont;

    if (contaminants.empty())
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
            String key = (pep_hit.getSequence().toUnmodifiedString());
            if (digested_db_.find(key) == digested_db_.end())                    //Peptide is not in Contaminant database
            {
              int x = 0;
              ++total;
              sum_total += f.getIntensity();
              pep_hit.setMetaValue("is_contaminant", x);
            }
            else                                                                 //Peptide is Contaminant
            {
              int x = 1;
              ++total;
              ++cont;
              sum_total += f.getIntensity();
              sum_cont += f.getIntensity();
              pep_hit.setMetaValue("is_contaminant", x);
            }
          }
        }
      }
      results_.push_back(std::make_tuple((cont / double(total)), (sum_cont / sum_total)));
    }
}

const std::vector<std::tuple<double, double>>& Contaminants::getResults()
{
  return results_;
}

int main()
{

}
//if (std::find(protein.peptides.begin(), protein.peptides.end(), (pep_hit.getSequence())) != protein.peptides.end())
//{
//int x = 1;
//++total;
//++cont;
//sum_total += f.getIntensity();
//sum_cont += f.getIntensity();
//pep_hit.setMetaValue("is_contaminant", x);

//}
//else
//{
//int x = 0;
//++total;
//sum_total += f.getIntensity();
//pep_hit.setMetaValue("is_contaminant", x);
//}

//           for (auto protein : digested_db_)
//          {
//            for (auto seq : protein.peptides)
//           {
//if (pep_hit.getSequence() == seq)
//{
//int x = 1;
//++total;
//++cont;
//sum_total += f.getIntensity();
//sum_cont += f.getIntensity();
//pep_hit.setMetaValue("is_contaminant", x);
//}
//else
//{
//}
//}
//}

