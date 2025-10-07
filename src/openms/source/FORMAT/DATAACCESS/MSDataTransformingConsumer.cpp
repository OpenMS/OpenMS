// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>

#include <utility>

namespace OpenMS
{
  /**
        @brief Constructor
      */
      MSDataTransformingConsumer::MSDataTransformingConsumer()
       : lambda_spec_(nullptr),
         lambda_chrom_(nullptr),
         lambda_exp_settings_(nullptr)
      {
      }

      /// Default destructor
      MSDataTransformingConsumer::~MSDataTransformingConsumer()
      = default;

      void MSDataTransformingConsumer::setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */)
      {
      }

      void MSDataTransformingConsumer::consumeSpectrum(SpectrumType& s)
      {
        // apply the given function to it (unless nullptr)
        if (lambda_spec_) lambda_spec_(s);
      }

      void MSDataTransformingConsumer::setSpectraProcessingFunc( std::function<void (SpectrumType&)> f_spec )
      {
        lambda_spec_ = std::move(f_spec);
      }

      void MSDataTransformingConsumer::consumeChromatogram(ChromatogramType & c)
      {
        // apply the given function to it (unless nullptr)
        if (lambda_chrom_) lambda_chrom_(c);
      }

      void MSDataTransformingConsumer::setChromatogramProcessingFunc( std::function<void (ChromatogramType&)> f_chrom )
      {
        lambda_chrom_ = std::move(f_chrom);
      }
      
      void MSDataTransformingConsumer::setExperimentalSettingsFunc( std::function<void (const OpenMS::ExperimentalSettings&)> f_exp_settings )
      {
        lambda_exp_settings_ = std::move(f_exp_settings);
      }

      void MSDataTransformingConsumer::setExperimentalSettings(const OpenMS::ExperimentalSettings& es)
      { 
        // apply the given function to it (unless nullptr)
        if (lambda_exp_settings_)
        {
          lambda_exp_settings_(es);
        }        
      }
} // namespace OpenMS
