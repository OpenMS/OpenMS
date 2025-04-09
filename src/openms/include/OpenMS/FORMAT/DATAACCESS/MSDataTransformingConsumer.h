// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

namespace OpenMS
{

    /**
      @brief Transforming consumer of MS data

      Is able to work on spectrum/chromatogram on the fly while it is read using a
      lambda function provided by the user.
      The lambda can transform the spectrum/chromatogram in-place or just extract some
      information and store it in locally captured variables. 
      Just make sure the captured variables are still in scope when this consumer is used.

    */
    class OPENMS_DLLAPI MSDataTransformingConsumer :
      public Interfaces::IMSDataConsumer
    {

    public:

      /**
        @brief Constructor
      */
      MSDataTransformingConsumer();

      /// Default destructor
      ~MSDataTransformingConsumer() override;

      /// ignored
      void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) override;

      void consumeSpectrum(SpectrumType& s) override;

      /**
        @brief Sets the lambda function to be called for every spectrum which is passed to this interface

        The lambda can contain locally captured variables, which allow to change some state. Make sure
        that the captured variables are still in scope (i.e. not destroyed at the point of calling the lambda).
        Pass a nullptr if the spectrum should be left unchanged.
      */
      virtual void setSpectraProcessingFunc( std::function<void (SpectrumType&)> f_spec );

      void consumeChromatogram(ChromatogramType& c) override;

      /**
        @brief Sets the lambda function to be called for every chromatogram which is passed to this interface

        The lambda can contain locally captured variables, which allow to change some state. Make sure
        that the captured variables are still in scope (i.e. not destroyed at the point of calling the lambda).
        Pass a nullptr if the chromatogram should be left unchanged.
      */
      virtual void setChromatogramProcessingFunc( std::function<void (ChromatogramType&)> f_chrom );

      /**
        @brief Sets the lambda function to be called when setExperimentalSettings is called via this interface

        The lambda can contain locally captured variables, which allow to change some state. Make sure
        that the captured variables are still in scope (i.e. not destroyed at the point of calling the lambda).
        Pass a nullptr if nothing should happen (default).
      */
      virtual void setExperimentalSettingsFunc( std::function<void (const OpenMS::ExperimentalSettings&)> f_exp_settings );

      /// ignored
      void setExperimentalSettings(const OpenMS::ExperimentalSettings&) override;

    protected:
      std::function<void (SpectrumType&)> lambda_spec_;
      std::function<void (ChromatogramType&)> lambda_chrom_;
      std::function<void (const OpenMS::ExperimentalSettings&)> lambda_exp_settings_;
    };

} //end namespace OpenMS

