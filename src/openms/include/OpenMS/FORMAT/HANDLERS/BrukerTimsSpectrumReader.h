// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <memory>
#include <string>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief The spectra BrukerTimsFile::load() gives for a .d directory, listed without reading their peaks.

      open() resolves the configuration as load() does and lists the spectra of load() in its order: their retention
      time, MS level and native ID, and what each is built from. It reads the frame table and analysis.tdf, but no
      peaks. A BrukerTimsSpectrumReader builds any of the spectra on request.

      Some spectra of the list turn out to have no peaks only when they are built, and load() leaves those out: a DDA
      precursor or a DIA window without peaks (or without centroids), and a DDA MS2 spectrum whose ion mobility lies
      outside Config::im_min/im_max. They stay in the list (the reader gives them without peaks), so load() gives the
      spectra of the list without those.

      Implemented in BrukerTimsFile.cpp, next to the loaders whose code it shares. Used by OnDiscBrukerTimsExperiment.
    */
    class BrukerTimsRun
    {
    public:
      /**
        @brief Lists the spectra that BrukerTimsFile::load(@p path, exp, @p config) gives

        @throws Exception::IllegalArgument for Config::calibrate with the linear m/z model: its OLS fit needs all
                precursors of the run decoded
        @throws the exceptions of BrukerTimsFile::load()
      */
      static std::shared_ptr<const BrukerTimsRun> open(const std::string& path, const BrukerTimsFile::Config& config);

      ~BrukerTimsRun();

      /// The number of spectra
      Size size() const;

      /// Retention time (in seconds) of spectrum @p index
      double getRT(Size index) const;

      /// MS level of spectrum @p index
      UInt getMSLevel(Size index) const;

      /// Native ID of spectrum @p index
      const std::string& getNativeID(Size index) const;

      /// The run-level metadata, as load() sets it
      const ExperimentalSettings& getSettings() const;

      /// The spectrum @p index without peaks, from the list alone (see BrukerTimsSpectrumReader::readMetadata())
      void readMetadata(Size index, MSSpectrum& spectrum) const;

      /// Defined in BrukerTimsFile.cpp
      struct Impl;

    private:
      BrukerTimsRun();
      std::unique_ptr<Impl> impl_;
      friend class BrukerTimsSpectrumReader;
    };

    /**
      @brief Builds the spectra of a BrukerTimsRun, one at a time, in any order.

      Each reader opens its own handle on the data (when it reads the first spectrum), so different readers of one run
      can be used from different threads. A single reader must not.
    */
    class BrukerTimsSpectrumReader
    {
    public:
      explicit BrukerTimsSpectrumReader(std::shared_ptr<const BrukerTimsRun> run);
      ~BrukerTimsSpectrumReader();

      BrukerTimsSpectrumReader(const BrukerTimsSpectrumReader&) = delete;
      BrukerTimsSpectrumReader& operator=(const BrukerTimsSpectrumReader&) = delete;

      /// Spectrum @p index as load() gives it, or without peaks if load() leaves it out (see BrukerTimsRun)
      void read(Size index, MSSpectrum& spectrum);

    private:
      struct Impl;
      std::unique_ptr<Impl> impl_;
    };
  } // namespace Internal
} // namespace OpenMS

#endif // WITH_OPENTIMS
