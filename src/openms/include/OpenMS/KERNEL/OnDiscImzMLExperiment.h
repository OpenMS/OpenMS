// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstddef>
#include <cstdint>
#include <memory>

namespace OpenMS
{

  /**
    @brief Random-access, on-disc reader for imzML mass spectrometry imaging datasets.

    @ingroup Kernel

    Analogous to @p OnDiscMSExperiment for indexed mzML files, but designed
    specifically for imzML's two-file layout (.imzML + .ibd).

    @b Workflow:
    @code
      OnDiscImzMLExperiment exp;
      exp.open("tissue.imzML");           // parse XML index only — no peak I/O

      std::cout << exp.getNrSpectra();    // instant, no IBD peak reads
      const auto& meta = exp.getImzMLMeta();

      MSSpectrum s = exp.getSpectrum(42);              // one fseek + fread per array
      MSSpectrum s = exp.getSpectrumAtCoord(3, 7);    // pixel lookup (1-based coords)
      const ImzMLSpectrumIndex& e = exp.getIndex(42); // IBD offsets, no read
    @endcode

    @b Design.
    The PIMPL idiom hides all implementation details (IBD FILE*, index vector,
    coordinate map) so that this header has no dependency on expat, Xerces, or
    the native @p imzml:: parser library.  The index is built by running
    @p ImzMLFile::loadSpectraIndex() during @p open() — no peak arrays
    are decoded until @p getSpectrum() is called.

    @note This class is not thread-safe.  If concurrent reads are needed,
    construct one instance per thread.

    @see ImzMLFile, ImzMLMeta, ImzMLSpectrumIndex, OnDiscMSExperiment
  */
  class OPENMS_DLLAPI OnDiscImzMLExperiment
  {
  public:

    OnDiscImzMLExperiment();
    ~OnDiscImzMLExperiment();           // needs Impl complete — defined in .cpp

    OnDiscImzMLExperiment(const OnDiscImzMLExperiment&)            = delete;
    OnDiscImzMLExperiment& operator=(const OnDiscImzMLExperiment&) = delete;

    OnDiscImzMLExperiment(OnDiscImzMLExperiment&&)            noexcept;
    OnDiscImzMLExperiment& operator=(OnDiscImzMLExperiment&&) noexcept;

    // ------------------------------------------------------------------
    // Open
    // ------------------------------------------------------------------

    /**
      @brief Open an imzML dataset: parse the XML index and open the .ibd.

      The .ibd path is inferred from @p imzml_path unless @p ibd_path is
      provided explicitly.  No peak data is read.

      @param imzml_path  Path to the .imzML file.
      @param ibd_path    Optional override for the .ibd path.
      @throws Exception::FileNotFound if either file cannot be opened.
      @throws Exception::ParseError   on malformed XML.
    */
    void open(const String& imzml_path, const String& ibd_path = "");

    /**
      @brief Close the companion .ibd file and release on-disc resources.

      After @p close(), @p isOpen() returns @c false and spectrum decode calls fail
      until @p open() is invoked again.
    */
    void close() noexcept;

    // ------------------------------------------------------------------
    // State
    // ------------------------------------------------------------------

    bool        isOpen()       const noexcept;
    std::size_t getNrSpectra() const noexcept;
    std::size_t size()         const noexcept { return getNrSpectra(); }

    // ------------------------------------------------------------------
    // Index access — no IBD read, O(1)
    // ------------------------------------------------------------------

    /**
      @brief Return the index entry for spectrum @p i.

      Contains pixel coordinates and byte offsets — no IBD read performed.

      @param i  0-based spectrum index.
      @throws Exception::IndexOverflow if @p i >= getNrSpectra().
    */
    const ImzMLSpectrumIndex& getIndex(std::size_t i) const;

    // ------------------------------------------------------------------
    // On-demand spectrum decode — one fseek + fread per call
    // ------------------------------------------------------------------

    /**
      @brief Decode and return spectrum @p i from the .ibd file.

      The returned MSSpectrum contains peak arrays decoded from the .ibd file
      and pixel coordinates as MetaValues (@p imzml:x, @p imzml:y, @p imzml:z).
      mzML scan metadata (RT, MS level, etc.) is not loaded in on-disc mode.

      @param i  0-based spectrum index.
      @throws Exception::IndexOverflow if @p i >= getNrSpectra().
      @throws Exception::FileNotFound  if the .ibd is not open.
    */
    MSSpectrum getSpectrum(std::size_t i) const;

    /// Sugar for @p getSpectrum(i).
    MSSpectrum operator[](std::size_t i) const { return getSpectrum(i); }

    /**
      @brief Return the spectrum at pixel coordinate (x, y[, z]).

      Builds a coordinate→index map on first call (O(n)), then O(1).

      @param x  Pixel column (1-based).
      @param y  Pixel row    (1-based).
      @param z  Depth slice  (1-based, default 1).
      @throws Exception::ElementNotFound if no spectrum exists at (x, y, z).
    */
    MSSpectrum getSpectrumAtCoord(uint32_t x, uint32_t y, uint32_t z = 1) const;

    // ------------------------------------------------------------------
    // Dataset-level metadata
    // ------------------------------------------------------------------

    /// Returns imaging metadata parsed during open() — no IBD reads.
    const ImzMLMeta& getImzMLMeta() const noexcept;

    /// Shorthand for getImzMLMeta().max_count_x.
    uint32_t gridWidth()  const noexcept;

    /// Shorthand for getImzMLMeta().max_count_y.
    uint32_t gridHeight() const noexcept;

  private:

    struct Impl;
    std::unique_ptr<Impl> pimpl_;
  };

} // namespace OpenMS
