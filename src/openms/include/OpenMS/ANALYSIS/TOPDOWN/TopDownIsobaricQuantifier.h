// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeon $
// $Authors: Kyowon Jeong $
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

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

  private:
    /// The quantification method used for the dataset to be analyzed.
    std::map<String, std::unique_ptr<IsobaricQuantitationMethod>> quant_methods_;

    void addMethod_(std::unique_ptr<IsobaricQuantitationMethod> ptr)
    {
      std::string internal_name = ptr->getMethodName();
      quant_methods_[internal_name] = std::move(ptr);
    }
  };
} // namespace OpenMS