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

#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>


namespace OpenMS
{

    /**
      @brief Transforming and cached writing consumer of MS data

      Is able to transform a spectrum on the fly while it is read using a
      function pointer that can be set on the object. The spectra is then
      cached to disk using the functions provided in CachedMzMLHandler.
    */
    class OPENMS_DLLAPI MSDataCachedConsumer :
      public Internal::CachedMzMLHandler,
      public Interfaces::IMSDataConsumer
    {
      typedef MSSpectrum SpectrumType;
      typedef MSChromatogram ChromatogramType;

    public:

      /**
        @brief Constructor

        Opens the output file and writes the header.

        @param filename The output file name to which data is written
        @param clearData Whether to clear the spectral and chromatogram data
        after writing (only keep meta-data)

        @note Clearing data from spectra and chromatograms also clears float
        and integer data arrays associated with the structure as these are
        written to disk as well.

      */
      MSDataCachedConsumer(const String& filename, bool clearData=true);

      /**
        @brief Destructor

        Closes the output file and writes the footer.
      */
      ~MSDataCachedConsumer() override;

      /**
        @brief Write a spectrum to the output file

        @note May delete data from spectrum (if clearData is set)
      */
      void consumeSpectrum(SpectrumType & s) override;

      /**
        @brief Write a chromatogram to the output file

        @note May delete data from chromatogram (if clearData is set)
      */
      void consumeChromatogram(ChromatogramType & c) override;

      void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) override {;}

      void setExperimentalSettings(const ExperimentalSettings& /* exp */) override {;}

    protected:
      std::ofstream ofs_;
      bool clearData_;
      Size spectra_written_;
      Size chromatograms_written_;

    };

} //end namespace OpenMS

