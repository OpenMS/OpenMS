//--------------------------------------------------------------------------
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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
  /**
  @brief
  @ingroup Topdown
  */

  class OPENMS_DLLAPI TopDownIsobaricQuantifier : public DefaultParamHandler
  {
  public:

    /// constructor
    TopDownIsobaricQuantifier();

    /// destructor
    ~TopDownIsobaricQuantifier() override = default;

    /// copy constructor
    TopDownIsobaricQuantifier(const TopDownIsobaricQuantifier&);

    /// move constructor
    TopDownIsobaricQuantifier(TopDownIsobaricQuantifier&& other) = default;

    /// assignment operator
    TopDownIsobaricQuantifier& operator=(const TopDownIsobaricQuantifier& other);

    /**
       @brief
       @param exp
       */
    void quantify(const MSExperiment& exp, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features);
    const FLASHDeconvHelperStructs::IsobaricQuantities getQuantities(int scan) const;

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

  private:
    /// The quantification method used for the dataset to be analyzed.
    std::map<String, std::unique_ptr<IsobaricQuantitationMethod>> quant_methods_;

    /// peak group information is stored in here for tracing
    std::map<int, FLASHDeconvHelperStructs::IsobaricQuantities> quantities_;

    void addMethod_(std::unique_ptr<IsobaricQuantitationMethod> ptr)
    {
      std::string internal_name = ptr->getMethodName();
      quant_methods_[internal_name] = std::move(ptr);
    }

  };
} // namespace OpenMS