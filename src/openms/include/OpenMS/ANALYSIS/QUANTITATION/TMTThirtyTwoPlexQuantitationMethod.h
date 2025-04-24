// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// // --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $, Julianus Pfeuffer $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

namespace OpenMS{
    /**
        @brief TMT 32plex quantitation to be used for the IsobaricQuantitation.
        @htmlinclude OpenMS_TMTThirtyTwoPlexQuantitationMethod.parameters
        */
    class OPENMS_DLLAPI TMTThirtyTwoPlexQuantitationMethod : public IsobaricQuantitationMethod
    {
        public:
            /// Default Constructor
            TMTThirtyTwoPlexQuantitationMethod();
            /// Destructor
            ~TMTThirtyTwoPlexQuantitationMethod() override = default;
            /// Copy Constructor
            TMTThirtyTwoPlexQuantitationMethod(const TMTThirtyTwoPlexQuantitationMethod& other);

            /// Assignment operator
            TMTThirtyTwoPlexQuantitationMethod& operator=(const TMTThirtyTwoPlexQuantitationMethod& rhs);

            /// @brief Methods to implement from IsobaricQuantitationMethod
            /// @{
            const String& getMethodName() const override;
            const IsobaricChannelList& getChannelInformation() const override;
            Size getNumberOfChannels() const override;
            Matrix<double> getIsotopeCorrectionMatrix() const override;
            Size getReferenceChannel() const override;
            /// @}
        private:
            /// The actual information on the different TMT 32plex channels.
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
} // namespace OpenMS