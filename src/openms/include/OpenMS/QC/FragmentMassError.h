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

#pragma once

#include "OpenMS/QC/QCBase.h"
#include <vector>

namespace OpenMS
{
  class FeatureMap;
  class MSExperiment;
  
  class OPENMS_DLLAPI FragmentMassError : public QCBase
  {
  public:
    enum class ToleranceUnit
    {
      PPM,
      DA,
      AUTO,
      SIZE_OF_TOLERANCEUNIT
    };
    /// strings corresponding to enum ToleranceUnit
    static const std::string names_of_toleranceUnit[];

    /// Default constructor
    FragmentMassError() = default;

    /// Destructor
    virtual ~FragmentMassError() = default;

    /**
     * @brief Structure for storing results: average and variance of all FragmentMassErrors in ppm
     */
    struct FMEStatistics
    {
      double average_ppm = 0;
      double variance_ppm = 0;
    };

    /**
     * @brief computes FragmentMassError in ppm and Dalton (only of the first PeptideHit of each PepID)
     *
     * Stores average and variance of FragmentMassErrors in ppm as a struct in a vector.
     * Each FragmentMassError (in ppm) is stored in the first PeptideHit of the corresponding PeptideIdentification as metavalue "fragment_mass_error_ppm".
     * Each FragmentMassError (in Da) is stored in the first PeptideHit of the corresponding PeptideIdentification as metavalue "fragment_mass_error_da".
     *
     * @param fmap Input FeatureMap for annotation and data for theoretical spectra
     * @param exp Input MSExperiment for MS2 spectra; spectra should be sorted (ascending RT)
     * @param map_to_spectrum Map to find index of spectrum given by meta value at PepID
     * @param tolerance Search window for matching peaks; distance has to be lower than tolerance value
     * @param tolerance_unit Tolerance in ppm or Dalton (if auto was chosen, the unit and value will taken from FeatureMap metadata)
     * @throws Exceptions::MissingInformation If fragment mass tolerance is missing in metadata of FeatureMap
     * @throws Exception::InvalidParameter PeptideID is missing meta value 'spectrum_reference'
     * @throws Exception::IllegalArgument Spectrum for a PepID has ms-level of 1
     * @throws Exception::MissingInformation If no fragmentation method given in a MS2 precursor
     * @throws Exception::InvalidParameter If the fragmentation method is not ECD, ETD, CID or HCD
     */
    void compute(FeatureMap& fmap, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum, ToleranceUnit tolerance_unit = ToleranceUnit::AUTO, double tolerance = 20);

    /// returns the name of the metric
    const String& getName() const override;
    
    /// returns results
    const std::vector<FMEStatistics>& getResults() const;


    /**
    * @brief Returns the input data requirements of the compute(...) function
    * @return Status for RAWMZML and POSTFDRFEAT
    */
    QCBase::Status requires() const override;


  private:
    /// name of the metric
    const String name_ = "FragmentMassError";
    
    /// container that stores results
    std::vector<FMEStatistics> results_{};
  };

} //namespace OpenMS
