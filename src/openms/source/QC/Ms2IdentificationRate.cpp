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
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#include <OpenMS/QC/Ms2IdentificationRate.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>

namespace OpenMS
{
  //computes number of peptide identifications, number of ms2 spectra and ratio
  //data is stored in vector of structs
  void Ms2IdentificationRate::compute(const FeatureMap& feature_map,const MSExperiment& exp, bool force_fdr)
  {
    if (exp.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSExperiment is empty");
    }

    // counts ms2 spectra
    UInt64 ms2_level_counter{};
    for (auto const& spec : exp.getSpectra())
    {
      if (spec.getMSLevel() == 2)
      {
        ++ ms2_level_counter;
      }
    }

    if (ms2_level_counter == 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No MS2 spectra found");
    }

    //counts peptideIdentifications
    UInt64 peptide_identification_counter{};

    auto lam = [force_fdr, &peptide_identification_counter](const PeptideIdentification& pep_id)
    {
      if (pep_id.getHits().empty())
      {
        return;
      }
      if (force_fdr)
      {
        ++peptide_identification_counter;
        return;
      }
      if (!(pep_id.getHits()[0].metaValueExists("target_decoy")))
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FDR was not made. If you want to continue without FDR use -MS2_id_rate:force_no_fdr");
      }
      if (pep_id.getHits()[0].getMetaValue("target_decoy") == "target")
      {
        ++peptide_identification_counter;
        return;
      }
    };

    //iterates through all PeptideIdentifications in FeatureMap, applies lambda function lam to all of them
    QCBase::iterateFeatureMap(feature_map, lam);

    if (ms2_level_counter < peptide_identification_counter)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There are more Identifications than MS2 spectra. Please check your data.");
    }

    //computes ratio
    double ratio = (double) peptide_identification_counter / ms2_level_counter;

    // struct that is made to store results
    IdentificationRateData id_rate_data;

    //store results
    id_rate_data.num_peptide_identification = peptide_identification_counter;
    id_rate_data.num_ms2_spectra = ms2_level_counter;
    id_rate_data.identification_rate = ratio;

    rate_result_.push_back(id_rate_data);
  }

  
  const String& Ms2IdentificationRate::getName() const
  {
    return name_;
  }
  

  const std::vector<OpenMS::Ms2IdentificationRate::IdentificationRateData>& Ms2IdentificationRate::getResults() const
  {
    return rate_result_;
  }


  QCBase::Status Ms2IdentificationRate::requires() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS
