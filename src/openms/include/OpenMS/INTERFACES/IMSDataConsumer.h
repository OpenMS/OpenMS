// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  class MSSpectrum;
  class MSChromatogram;
  class ExperimentalSettings;

namespace Interfaces
{

    /**
      @brief The interface of a consumer of spectra and chromatograms

      The data consumer is able to consume data of type MSSpectrum and
      MSChromatogram and process them (it may modify the spectra). The consumer
      interface may be used when data is generated sequentially (e.g. by
      reading from disc) and needs to be processed as fast as possible without
      ever holding the full set of data in memory.

      The consumer expects to be informed about the number of spectra and
      chromatograms to consume and potentially about the ExperimentalSettings
      @a before_consuming any spectra. This can be critical for consumers who
      write data to disk. Depending on the implementation, an exception may
      occur if the ExperimentalSettings and the size of the experiment are not
      set before consuming any spectra.

      Implementations in OpenMS can be found in OpenMS/FORMAT/DATAACCESS

      @note The member functions setExpectedSize and setExperimentalSettings
      are expected to be called before consuming starts.

    */
    class OPENMS_DLLAPI IMSDataConsumer
    {
    public:
      typedef MSSpectrum SpectrumType;
      typedef MSChromatogram ChromatogramType;

      virtual ~IMSDataConsumer() {}

      /**
        @brief Consume a spectrum

        The spectrum will be consumed by the implementation and possibly modified.

        @note The implementation might not allow to consume spectra and chromatograms in any order

        @param s The spectrum to be consumed
      */
      virtual void consumeSpectrum(SpectrumType& s) = 0;

      /**
        @brief Consume a chromatogram

        The chromatogram will be consumed by the implementation and possibly modified.

        @note The implementation might not allow to consume spectra and chromatograms in any order

        @param c The chromatogram to be consumed
      */
      virtual void consumeChromatogram(ChromatogramType& c) = 0;

      /**
        @brief Set expected size of spectra and chromatograms to be consumed.

        Some implementations might care about the number of spectra and
        chromatograms to be consumed and need to be informed about this
        (usually before consuming starts).

        @note Calling this method is optional but good practice.

        @param expectedSpectra Number of spectra expected
        @param expectedChromatograms Number of chromatograms expected
      */
      virtual void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms) = 0;

      /**
        @brief Set experimental settings (meta-data) of the data to be consumed

        Some implementations might need to know about the meta-data (or the
        context) of the spectra and chromatograms to be consumed. This method
        allows them learn this.

        @note Calling this method is optional but good practice.

        @param exp Experimental settings meta data for the data to be consumed
      */
      virtual void setExperimentalSettings(const ExperimentalSettings& exp) = 0;
    };

    typedef IMSDataConsumer IMSDataConsumer;

} //end namespace Interfaces
} //end namespace OpenMS

