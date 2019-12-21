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
#include <algorithm>


using namespace std;


namespace OpenMS
{

  void Contaminants::compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants)
  {
    // empty FeatureMap
    if (features.empty())
    {
      OPENMS_LOG_WARN << "FeatureMap is empty" << "\n";
    }
    // empty contaminants database
    if (contaminants.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No contaminants provided.");
    }
    // fill the unordered set once with the digested contaminants database
    if (digested_db_.empty())
    {
      if (features.getProteinIdentifications().empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No proteinidentifications in FeatureMap.");
      }
      ProteaseDigestion digestor;
      String enzyme = features.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme.getName();

      // no enzyme is given
      if (enzyme == "unknown_enzyme")
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No digestion enzyme in FeatureMap detected. No computation possible.");
      }

      digestor.setEnzyme(enzyme);

      // get the missed cleavages for the digestor. If none are given, its default is 0.
      UInt missed_cleavages(features.getProteinIdentifications()[0].getSearchParameters().missed_cleavages);
      digestor.setMissedCleavages(missed_cleavages);

      // digest the contaminants database and add the peptides into the unordered set
      for (const FASTAFile::FASTAEntry& fe : contaminants)
      {
        vector<AASequence> current_digest;
        digestor.digest(AASequence::fromString(fe.sequence), current_digest);

        // fill unordered set digested_db_ with digested sequences
        for (auto const& s : current_digest)
        {
          digested_db_.insert(s.toUnmodifiedString());
        }
      }
    }
    Int64 total = 0;
    Int64 cont = 0;
    double sum_total = 0.0;
    double sum_cont = 0.0;
    Int64 feature_has_no_sequence = 0;

    // Check if peptides of featureMap are contaminants or not and add is_contaminant = 0/1 to the first hit of the peptideidentification.
    // If so, raise contaminants ratio.
    for (auto& f : features)
    {
      if (f.getPeptideIdentifications().empty())
      {
        ++feature_has_no_sequence;
        continue;
      }
      for (auto& id : f.getPeptideIdentifications())
      {
        // peptideidentifications in feature f is not empty
        if (id.getHits().empty())
        {
          ++feature_has_no_sequence;
          continue;
        }

        // the one existing peptideidentification has at least one getHits entry
        PeptideHit& pep_hit = id.getHits()[0];
        String key = (pep_hit.getSequence().toUnmodifiedString());
        this->compare_(key, pep_hit, total, cont, sum_total, sum_cont, f.getIntensity());
      }
    }
    // save the contaminants ratio in object before searching through the unassigned peptideidentifications
    ContaminantsSummary final;
    final.assigned_contaminants_ratio = (cont / double(total));

    final.empty_features.first = feature_has_no_sequence;
    final.empty_features.second = features.size();

    UInt64 utotal = 0;
    UInt64 ucont = 0;

    // Change the assigned contaminants ratio to total contaminants ratio by adding the unassigned.
    // Additionally save the unassigned contaminants ratio and add the is_contaminant = 0/1 to the first hit of the unassigned peptideidentifications.
    for (auto& fu : features.getUnassignedPeptideIdentifications())
    {
      if ( fu.getHits().empty())
      {
        continue;
      }
      auto& fu_hit = fu.getHits()[0];
      String key = (fu_hit.getSequence().toUnmodifiedString());
      ++utotal;

      // peptide is not in contaminant database
      if (!digested_db_.count(key))
      {
        fu_hit.setMetaValue("is_contaminant", 0);
        continue;
      }

      // peptide is contaminant
      ++ucont;
      fu_hit.setMetaValue("is_contaminant", 1);
    }
    total += utotal;
    cont += ucont;

    // save all ratios and the intensity to the object
    final.all_contaminants_ratio = (cont / double(total));
    final.unassigned_contaminants_ratio = (ucont / double(utotal));
    final.assigned_contaminants_intensity_ratio = (sum_cont / sum_total);


    // add the object to the results vector
    results_.push_back(final);
  }
  
  const String& Contaminants::getName() const
  {
    return name_;
  }
  
  const std::vector<Contaminants::ContaminantsSummary>& Contaminants::getResults()
  {
    return results_;
  }
  

  // Check if peptide is in contaminants database or not and add the is_contaminant = 0/1.
  // If so, raise the contaminant ratio.
  void Contaminants::compare_(const String& key,
                              PeptideHit& pep_hit,
                              Int64& total,
                              Int64& cont,
                              double& sum_total,
                              double& sum_cont,
                              double intensity)
  {
    ++total;
    sum_total += intensity;
    // peptide is not in contaminant database
    if (!digested_db_.count(key))
    {
      pep_hit.setMetaValue("is_contaminant", 0);
      return;
    }
    // peptide is contaminant
    ++cont;
    sum_cont += intensity;
    pep_hit.setMetaValue("is_contaminant", 1);
  }

  QCBase::Status Contaminants::requires() const
  {
    return (QCBase::Status(QCBase::Requires::POSTFDRFEAT) | QCBase::Requires::CONTAMINANTS);
  }


}
