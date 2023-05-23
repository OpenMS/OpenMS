// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/MissedCleavages.h>
#include <iostream>

namespace OpenMS
{
  typedef std::map<UInt32, UInt32> MapU32;
  // digests the Sequence in PeptideHit and counts the number of missed cleavages
  void MissedCleavages::get_missed_cleavages_from_peptide_identification_(const ProteaseDigestion& digestor, MapU32& result, const UInt32& max_mc, PeptideIdentification& pep_id)
  {
    if (pep_id.getHits().empty())
    {
      OPENMS_LOG_WARN << "There is a Peptideidentification(RT: " << pep_id.getRT() << ", MZ: " << pep_id.getMZ() << ") without PeptideHits.\n";
      return;
    }
    std::vector<AASequence> digest_output;
    digestor.digest(pep_id.getHits()[0].getSequence(), digest_output);
    UInt32 num_mc = UInt32(digest_output.size() - 1);

    // warn if number of missed cleavages is greater than allowed maximum number of missed cleavages
    if (num_mc > max_mc)
    {
      OPENMS_LOG_WARN << "Observed number of missed cleavages: " << num_mc << " is greater than: " << max_mc
                      << " the allowed maximum number of missed cleavages during MS2-Search in: " << pep_id.getHits()[0].getSequence() << "\n";
    }

    ++result[num_mc];

    pep_id.getHits()[0].setMetaValue("missed_cleavages", num_mc);
  };

  void MissedCleavages::compute(std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)
  {
    MapU32 result {};

    // Exception if ProteinIdentification is empty
    if (prot_ids.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing information in ProteinIdentifications.");
    }

    String enzyme = prot_ids[0].getSearchParameters().digestion_enzyme.getName();
    auto max_mc = prot_ids[0].getSearchParameters().missed_cleavages;

    // Exception if digestion enzyme is not given
    if (enzyme == "unknown_enzyme")
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No digestion enzyme in ID data detected. No computation possible.");
    }

    // create a digestor, which doesn't allow any missed cleavages
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(0);

    for (PeptideIdentification& pep_id : pep_ids)
    {
      get_missed_cleavages_from_peptide_identification_(digestor, result, max_mc, pep_id);
    }

    mc_result_.push_back(result);
  }


  void MissedCleavages::compute(FeatureMap& fmap)
  {
    MapU32 result {};

    bool has_pepIDs = QCBase::hasPepID(fmap);
    if (!has_pepIDs)
    {
      mc_result_.push_back(result);
      return;
    }

    // if the FeatureMap is empty, result is 0
    if (fmap.empty())
    {
      OPENMS_LOG_WARN << "FeatureXML is empty.\n";
      mc_result_.push_back(result);
      return;
    }

    // Exception if ProteinIdentification is empty
    if (fmap.getProteinIdentifications().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing information in ProteinIdentifications.");
    }

    String enzyme = fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme.getName();
    auto max_mc = fmap.getProteinIdentifications()[0].getSearchParameters().missed_cleavages;

    // Exception if digestion enzyme is not given
    if (enzyme == "unknown_enzyme")
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No digestion enzyme in FeatureMap detected. No computation possible.");
    }

    // create a digestor, which doesn't allow any missed cleavages
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(0);

    // small lambda function to apply get_missed_cleavages_from_peptide_identification on pep_ids
    auto l = [&](PeptideIdentification& pep_id) {
      get_missed_cleavages_from_peptide_identification_(digestor, result, max_mc, pep_id);
      return;
    };
    // iterate through all PeptideIdentifications of a given FeatureMap and applies the given lambda function
    fmap.applyFunctionOnPeptideIDs(l);

    mc_result_.push_back(result);
  }


  const String& MissedCleavages::getName() const
  {
    static const String& name = "MissedCleavages";
    return name;
  }


  const std::vector<MapU32>& MissedCleavages::getResults() const
  {
    return mc_result_;
  }


  QCBase::Status MissedCleavages::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS
