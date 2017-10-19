// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_INSTRUMENTSETTINGS_H
#define OPENMS_METADATA_INSTRUMENTSETTINGS_H

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
      UNKNOWN,                      ///< Unknown scan method
      MASSSPECTRUM,       ///< general spectrum type
      MS1SPECTRUM,              ///< full scan mass spectrum, is a "mass spectrum" @n Synonyms: 'full spectrum', 'Q1 spectrum', 'Q3 spectrum', 'Single-Stage Mass Spectrometry'
      MSNSPECTRUM,        ///< MS2+ mass spectrum, is a "mass spectrum"
      SIM,                              ///< Selected ion monitoring scan @n Synonyms: 'Multiple ion monitoring scan', 'SIM scan', 'MIM scan'
      SRM,                              ///< Selected reaction monitoring scan @n Synonyms: 'Multiple reaction monitoring scan', 'SRM scan', 'MRM scan'
      CRM,                              ///< Consecutive reaction monitoring scan @n Synonyms: 'CRM scan'
      CNG,                              ///< Constant neutral gain scan @n Synonyms: 'CNG scan'
      CNL,                              ///< Constant neutral loss scan @n Synonyms: 'CNG scan'
      PRECURSOR,                ///< Precursor ion scan
      EMC,                              ///< Enhanced multiply charged scan
      TDF,                              ///< Time-delayed fragmentation scan
      EMR,                              ///< Electromagnetic radiation scan @n Synonyms: 'EMR spectrum'
      EMISSION,                     ///< Emission scan
      ABSORPTION,               ///< Absorption scan
      SIZE_OF_SCANMODE
    };

    /// Names of scan modes
    static const std::string NamesOfScanMode[SIZE_OF_SCANMODE];

    /// Constructor
    InstrumentSettings();
    /// Copy constructor
    InstrumentSettings(const InstrumentSettings & source);
    /// Destructor
    ~InstrumentSettings();

    /// Assignment operator
    InstrumentSettings & operator=(const InstrumentSettings & source);

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

#endif // OPENMS_METADATA_INSTRUMENTSETTINGS_H
