// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

namespace OpenMS
{

    /**
      @brief Aggregates spectra by retention time

      This consumer will merge spectra passed to it that have the same
      retention time and will then pass them to the next consumer (see
      Constructor). Spectra are aggregated using
      SpectrumAddition::addUpSpectra() which merges the spectra.

    */
    class OPENMS_DLLAPI MSDataAggregatingConsumer :
      public Interfaces::IMSDataConsumer
    {

      Interfaces::IMSDataConsumer* next_consumer_;
      double previous_rt_;
      bool rt_initialized_;
      SpectrumType s_tmp;
      std::vector<SpectrumType> s_list;

    public:

      /**
        @brief Constructor

        @note This does not transfer ownership of the consumer
      */
      MSDataAggregatingConsumer(Interfaces::IMSDataConsumer* next_consumer) :
        next_consumer_(next_consumer),
        previous_rt_(0.0),
        rt_initialized_(false)
      {}

      /**
        @brief Constructor

        Flushes data to next consumer

        @note It is essential to not delete the underlying next_consumer before
        deleting this object, otherwise we risk a memory error
      */
      ~MSDataAggregatingConsumer() override;

      void setExpectedSize(Size, Size) override {}

      void consumeSpectrum(SpectrumType& s) override;

      void consumeChromatogram(ChromatogramType& c) override;

      void setExperimentalSettings(const OpenMS::ExperimentalSettings&) override {}

    };

} //end namespace OpenMS


