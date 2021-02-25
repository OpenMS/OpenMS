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
// $Maintainer: Tom Waschischeck$
// $Authors: Tom Waschischeck$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <vector>

namespace OpenMS
{
  class FeatureMap;
  class MSExperiment;
  class PeptideIdentification;
  
  class OPENMS_DLLAPI PSMCorrectness : public QCBase
  {
  public:
    enum class ToleranceUnit
    {
      AUTO,
      PPM,
      DA,
      SIZE_OF_TOLERANCEUNIT
    };
    /// strings corresponding to enum ToleranceUnit
    static const std::string names_of_toleranceUnit[];

    /// Default constructor
    PSMCorrectness() = default;

    /// Destructor
    virtual ~PSMCorrectness() = default;

    /**
     * @brief Structure for storing results: average and variance over all PSMs
     */
    struct Statistics
    {
      double average_correctness = 0;
      double variance_correctness = 0;
    };

    void compute(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20);

    void compute(std::vector<PeptideIdentification> & pep_ids, const ProteinIdentification::SearchParameters& search_params, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20);

    /// returns the name of the metric
    const String& getName() const override;
    
    /// returns results
    const std::vector<Statistics>& getResults() const;


    /**
    * @brief Returns the input data requirements of the compute(...) function
    * @return Status for RAWMZML and POSTFDRFEAT
    */
    QCBase::Status requires() const override;

  private:
    /// container that stores results
    std::vector<Statistics> results_{};
  };

} //namespace OpenMS