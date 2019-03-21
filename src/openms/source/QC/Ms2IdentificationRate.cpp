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

//#include <OpenMS/QC/Ms2IdentificationRate.h>
#include <include/OpenMS/QC/Ms2IdentificationRate.h>
#include <algorithm>
#include <iostream>



namespace OpenMS
{

  Ms2IdentificationRate::Ms2IdentificationRate() = default;
  Ms2IdentificationRate::~Ms2IdentificationRate() = default;

  Int64 Ms2IdentificationRate::countPeptideId_(std::vector<PeptideIdentification> peptide_id, bool force_fdr)
  {
    Int64 counter{};
    counter = count_if(peptide_id.begin(), peptide_id.end(), [force_fdr] (PeptideIdentification const & x)
    {
      if ( !x.getHits().empty() )
      {
        if (x.getHits()[0].metaValueExists("target_decoy") && !force_fdr)
        {
          if (x.getHits()[0].getMetaValue("target_decoy") == "target")
          {
            return true;
          }
          else
          {
            return false;
          }
        }
        else if (force_fdr)
        {
          return true;
        }
        else
        {
          LOG_ERROR << "Ms2IdentificationRate: no FDR was made. If you want to continue without FDR use -force";
          return false;
        }
      }
      else
      {
        LOG_WARN << "Ms2IdentificationRate: empty hits";
        return false;
      }
    });
    return counter;
  }

  //computes number of peptide identifications, number of ms2 spectra and ratio
  //data is stored in vector of structs
//  void Ms2IdentificationRate::compute(FeatureMap const & feature_map, MSExperiment const & exp, std::string file, bool force_fdr)
  void Ms2IdentificationRate::compute(FeatureMap const & feature_map, MSExperiment const & exp, bool force_fdr)
  {
    //checks if data exists
      if (feature_map.empty())
      {
        LOG_WARN << "Ms2IdentificationRate: FeatureXML is corrupted or empty";
      }
      if (exp.empty())
      {
        LOG_WARN << "Ms2IdentificationRate: mzML is corrupted or empty";
      }

//    if (file == "default")
//    {
//      std::cout << "Ms2IdentificationRate: There is no filename, you can enter it now: " << std::endl;
//      std::cin >> file;
//    }

    //count ms2 spectra
    Int64 ms2_level_counter{};

    for (auto const &spec : exp.getSpectra())
    {
      if (spec.getMSLevel() == 2)
      {
        ++ ms2_level_counter;
      }
    }

    //count peptideIdentifications
    Int64 peptide_identification_counter{};
    peptide_identification_counter += countPeptideId_(feature_map.getUnassignedPeptideIdentifications(), force_fdr);
    for (auto const &f : feature_map)
    {
      peptide_identification_counter += countPeptideId_(f.getPeptideIdentifications(), force_fdr);
    }

    //compute ratio
    double ratio{};
    ratio = (double) peptide_identification_counter / ms2_level_counter;

    //store results
//    id_rate_data_.filename = file;
    id_rate_data_.num_peptide_identification = peptide_identification_counter;
    id_rate_data_.num_ms2_spectra = ms2_level_counter;
    id_rate_data_.identification_rate = ratio;

    rate_result_.push_back(id_rate_data_);
  }


  std::vector<OpenMS::Ms2IdentificationRate::IdentificationRateData> Ms2IdentificationRate::getResults()
  {
    return rate_result_;
  }


  void Ms2IdentificationRate::clear()
  {
    rate_result_.clear();
  }


  QCBase::Status Ms2IdentificationRate::requires() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS