// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeon $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
  /**
  @brief Isobaric quantification for Top down proteomics. The report ion ratios from multiple MS2 are
   merged when their precursor masses belong to the same deconvolved features (MassFeature instances).
  @ingroup Topdown
  */

  class OPENMS_DLLAPI TopDownIsobaricQuantification : public DefaultParamHandler
  {
  public:
    /// constructor
    TopDownIsobaricQuantification();

    /// destructor
    ~TopDownIsobaricQuantification() override = default;

    /// copy constructor
    TopDownIsobaricQuantification(const TopDownIsobaricQuantification&);

    /// move constructor
    TopDownIsobaricQuantification(TopDownIsobaricQuantification&& other) = default;

    /// assignment operator
    TopDownIsobaricQuantification& operator=(const TopDownIsobaricQuantification& other);

    /**
     * Run quantification
     * @param exp the MS experiment
     * @param deconvolved_spectra deconvolved spectra for which the quantification will be carried out
     * @param mass_features mass features that are used to merge quantification results for the MS2 spectra from the same precursor mass
     */
    void quantify(const MSExperiment& exp, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHHelperClasses::MassFeature>& mass_features);

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

  private:
    /// The quantification method used for the dataset to be analyzed.
    std::map<String, std::unique_ptr<IsobaricQuantitationMethod>> quant_methods_;
    /// retain only fully quantified ratios
    bool only_fully_quantified_ = false;
    void addMethod_(std::unique_ptr<IsobaricQuantitationMethod> ptr)
    {
      std::string internal_name = ptr->getMethodName();
      quant_methods_[internal_name] = std::move(ptr);
    }
  };
} // namespace OpenMS