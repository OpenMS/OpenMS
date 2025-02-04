// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

namespace OpenMS
{
  /**
    @brief iTRAQ 4 plex quantitation to be used with the IsobaricQuantitation.
   
    @htmlinclude OpenMS_ItraqFourPlexQuantitationMethod.parameters
  */
  class OPENMS_DLLAPI ItraqFourPlexQuantitationMethod :
    public IsobaricQuantitationMethod
  {
public:
    /// Default c'tor
    ItraqFourPlexQuantitationMethod();

    /// d'tor
    ~ItraqFourPlexQuantitationMethod() override;

    /// Copy c'tor
    ItraqFourPlexQuantitationMethod(const ItraqFourPlexQuantitationMethod& other);

    /// Assignment operator
    ItraqFourPlexQuantitationMethod& operator=(const ItraqFourPlexQuantitationMethod& rhs);

    /// @brief Methods to implement from IsobaricQuantitationMethod
    /// @{

    const String& getMethodName() const override;

    const IsobaricChannelList& getChannelInformation() const override;

    Size getNumberOfChannels() const override;

    Matrix<double> getIsotopeCorrectionMatrix() const override;

    Size getReferenceChannel() const override;
    
    /// @}

private:
    /// the actual information on the different itraq4plex channels.
    IsobaricChannelList channels_;

    /// The name of the quantitation method.
    static const String name_;
    
    /// The reference channel for this experiment.
    Size reference_channel_;

protected:
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

    /// implemented for DefaultParamHandler
    void updateMembers_() override;
  };
} // namespace

