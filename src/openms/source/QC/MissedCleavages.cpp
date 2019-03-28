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
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/QC/MissedCleavages.h>
#include <iostream>

namespace OpenMS
{

  void MissedCleavages::compute(FeatureMap& fmap)
  {
    std::map<UInt64, UInt64> result{};

    if(fmap.empty())
    {
      LOG_WARN << "FeatureXML is empty.";
      mc_result_.push_back(result);
      return;
    }
    //ProteinIdentification kann leer sein

    if(fmap.getProteinIdentifications().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing information in ProteinIdentifications.");
    }

    String enzyme = fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme.getName();
    UInt64 max_mc = fmap.getProteinIdentifications()[0].getSearchParameters().missed_cleavages;

    if(enzyme == "unknown_enzyme")
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No digestion enzyme in FeatureMap detected. No computation possible.");
    }



    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(0);


    auto l = [&digestor, &result, &max_mc](PeptideIdentification& pep_id)
    {
      if(pep_id.getHits().empty())
      {
        return;
      }
      std::vector<AASequence> digest_output;
      digestor.digest(pep_id.getHits()[0].getSequence(), digest_output);
      UInt64 num_mc = digest_output.size() - 1;

      if (num_mc > max_mc)
      {
        LOG_WARN << "Number of missed cleavages is greater than the allowed maximum number of missed cleavages.";
      }

      if (result.count(num_mc) != 0)
      {
        ++result[num_mc];
      }
      else
      {
        result.insert(std::pair<UInt64, UInt64>(num_mc, 1));
      }

      pep_id.getHits()[0].setMetaValue("missed_cleavages", num_mc);
    };

    QCBase::iterateFeatureMap(fmap, l);

    mc_result_.push_back(result);
  }


  const std::vector<std::map<UInt64, UInt64>>& MissedCleavages::getResults() const
  {
    return mc_result_;
  }


  QCBase::Status MissedCleavages::requires() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML;
  }
}// namespace OpenMS