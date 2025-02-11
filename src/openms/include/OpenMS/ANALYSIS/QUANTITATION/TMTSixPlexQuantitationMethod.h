// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

namespace OpenMS
{
  /**
    @brief TMT 6plex quantitation to be used with the IsobaricQuantitation.

    @htmlinclude OpenMS_TMTSixPlexQuantitationMethod.parameters
  */
  class OPENMS_DLLAPI TMTSixPlexQuantitationMethod :
    public IsobaricQuantitationMethod
  {
public:
    /// Default c'tor
    TMTSixPlexQuantitationMethod();

    /// d'tor
    ~TMTSixPlexQuantitationMethod() override = default;

    /// Copy c'tor
    TMTSixPlexQuantitationMethod(const TMTSixPlexQuantitationMethod& other);

    /// Assignment operator
    TMTSixPlexQuantitationMethod & operator=(const TMTSixPlexQuantitationMethod& rhs);

    /// @brief Methods to implement from IsobaricQuantitationMethod
    /// @{

    const String& getMethodName() const override;

    const IsobaricChannelList& getChannelInformation() const override;

    Size getNumberOfChannels() const override;

    Matrix<double> getIsotopeCorrectionMatrix() const override;

    Size getReferenceChannel() const override;

    /// @}

  private:
    /// the actual information on the different tmt 6plex channels.
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

