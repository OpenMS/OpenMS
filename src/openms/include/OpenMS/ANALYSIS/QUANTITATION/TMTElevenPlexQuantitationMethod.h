// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
    @brief TMT 11plex quantitation to be used with the IsobaricQuantitation.

    @htmlinclude OpenMS_TMTSixPlexQuantitationMethod.parameters
  */
  class OPENMS_DLLAPI TMTElevenPlexQuantitationMethod :
    public IsobaricQuantitationMethod
  {
public:
    /// Default c'tor
    TMTElevenPlexQuantitationMethod();

    /// d'tor
    ~TMTElevenPlexQuantitationMethod() override = default;

    /// Copy c'tor
    TMTElevenPlexQuantitationMethod(const TMTElevenPlexQuantitationMethod& other);

    /// Assignment operator
    TMTElevenPlexQuantitationMethod & operator=(const TMTElevenPlexQuantitationMethod& rhs);

    /// @brief Methods to implement from IsobaricQuantitationMethod
    /// @{

    const String& getMethodName() const override;

    const IsobaricChannelList& getChannelInformation() const override;

    Size getNumberOfChannels() const override;

    Matrix<double> getIsotopeCorrectionMatrix() const override;

    Size getReferenceChannel() const override;

    /// @}

  private:
    /// the actual information on the different tmt11plex channels.
    IsobaricChannelList channels_;

    /// The name of the quantitation method.
    static const String name_;

    /// The reference channel for this experiment.
    Size reference_channel_;

    /// List of available channel names as they are presented to the user
    static const std::vector<std::string> channel_names_;

  protected:
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

    /// implemented for DefaultParamHandler
    void updateMembers_() override;

  };
} // namespace

