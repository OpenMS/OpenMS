// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ScanWindow.h>
#include <OpenMS/METADATA/IonSource.h>

namespace OpenMS
{
  /**
      @brief Description of the settings a MS Instrument was run with.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI InstrumentSettings :
    public MetaInfoInterface
  {
public:
    /// scan mode
    enum ScanMode
    {
      UNKNOWN,          ///< Unknown scan method
      MASSSPECTRUM,     ///< general spectrum type
      MS1SPECTRUM,      ///< full scan mass spectrum, is a "mass spectrum" @n Synonyms: 'full spectrum', 'Q1 spectrum', 'Q3 spectrum', 'Single-Stage Mass Spectrometry'
      MSNSPECTRUM,      ///< MS2+ mass spectrum, is a "mass spectrum"
      SIM,              ///< Selected ion monitoring scan @n Synonyms: 'Multiple ion monitoring scan', 'SIM scan', 'MIM scan'
      SRM,              ///< Selected reaction monitoring scan @n Synonyms: 'Multiple reaction monitoring scan', 'SRM scan', 'MRM scan'
      CRM,              ///< Consecutive reaction monitoring scan @n Synonyms: 'CRM scan'
      CNG,              ///< Constant neutral gain scan @n Synonyms: 'CNG scan'
      CNL,              ///< Constant neutral loss scan @n Synonyms: 'CNG scan'
      PRECURSOR,        ///< Precursor ion scan
      EMC,              ///< Enhanced multiply charged scan
      TDF,              ///< Time-delayed fragmentation scan
      EMR,              ///< Electromagnetic radiation scan @n Synonyms: 'EMR spectrum'
      EMISSION,         ///< Emission scan
      ABSORPTION,       ///< Absorption scan
      SIZE_OF_SCANMODE
    };

    /// Names of scan modes
    static const std::string NamesOfScanMode[SIZE_OF_SCANMODE];

    /// Constructor
    InstrumentSettings();
    /// Copy constructor
    InstrumentSettings(const InstrumentSettings &) = default;
    /// Move constructor
    InstrumentSettings(InstrumentSettings&&) = default;
    /// Destructor
    ~InstrumentSettings();

    /// Assignment operator
    InstrumentSettings & operator=(const InstrumentSettings &) = default;
    /// Move assignment operator
    InstrumentSettings& operator=(InstrumentSettings&&) & = default;

    /// Equality operator
    bool operator==(const InstrumentSettings & rhs) const;
    /// Equality operator
    bool operator!=(const InstrumentSettings & rhs) const;

    /// returns the scan mode
    ScanMode getScanMode() const;
    /// sets the scan mode
    void setScanMode(ScanMode scan_mode);

    /// return if this scan is a zoom (enhanced resolution) scan
    bool getZoomScan() const;
    /// sets if this scan is a zoom (enhanced resolution) scan
    void setZoomScan(bool zoom_scan);

    /// returns the polarity
    IonSource::Polarity getPolarity() const;
    /// sets the polarity
    void setPolarity(IonSource::Polarity polarity);

    /// returns a const reference to the m/z scan windows
    const std::vector<ScanWindow> & getScanWindows() const;
    /// returns a mutable reference to the m/z scan windows
    std::vector<ScanWindow> & getScanWindows();
    /// sets the m/z scan windows
    void setScanWindows(std::vector<ScanWindow>  scan_windows);

protected:
    ScanMode scan_mode_;
    bool zoom_scan_;
    IonSource::Polarity polarity_;
    std::vector<ScanWindow> scan_windows_;
  };
} // namespace OpenMS

