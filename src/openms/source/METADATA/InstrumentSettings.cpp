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

#include <OpenMS/METADATA/InstrumentSettings.h>

using namespace std;

namespace OpenMS
{
  const std::string InstrumentSettings::NamesOfScanMode[] = {"Unknown", "MassSpectrum", "MS1Spectrum", "MSnSpectrum", "SelectedIonMonitoring", "SelectedReactionMonitoring", "ConsecutiveReactionMonitoring", "ConstantNeutralGain", "ConstantNeutralLoss", "Precursor", "EnhancedMultiplyCharged", "TimeDelayedFragmentation", "ElectromagneticRadiation", "Emission", "Absorption"};

  InstrumentSettings::InstrumentSettings() :
    MetaInfoInterface(),
    scan_mode_(UNKNOWN),
    zoom_scan_(false),
    polarity_(IonSource::POLNULL),
    scan_windows_()
  {
  }

  InstrumentSettings::InstrumentSettings(const InstrumentSettings & source) :
    MetaInfoInterface(source),
    scan_mode_(source.scan_mode_),
    zoom_scan_(source.zoom_scan_),
    polarity_(source.polarity_),
    scan_windows_(source.scan_windows_)
  {
  }

  InstrumentSettings::~InstrumentSettings()
  {
  }

  InstrumentSettings & InstrumentSettings::operator=(const InstrumentSettings & source)
  {
    if (&source == this)
      return *this;

    scan_mode_ = source.scan_mode_;
    zoom_scan_  = source.zoom_scan_;
    polarity_ = source.polarity_;
    scan_windows_ = source.scan_windows_;
    MetaInfoInterface::operator=(source);

    return *this;
  }

  bool InstrumentSettings::operator==(const InstrumentSettings & rhs) const
  {
    return scan_mode_ == rhs.scan_mode_ &&
           zoom_scan_  == rhs.zoom_scan_ &&
           polarity_ == rhs.polarity_ &&
           scan_windows_ == rhs.scan_windows_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool InstrumentSettings::operator!=(const InstrumentSettings & rhs) const
  {
    return !(operator==(rhs));
  }

  InstrumentSettings::ScanMode InstrumentSettings::getScanMode() const
  {
    return scan_mode_;
  }

  void InstrumentSettings::setScanMode(InstrumentSettings::ScanMode scan_mode)
  {
    scan_mode_ = scan_mode;
  }

  IonSource::Polarity InstrumentSettings::getPolarity() const
  {
    return polarity_;
  }

  void InstrumentSettings::setPolarity(IonSource::Polarity polarity)
  {
    polarity_ = polarity;
  }

  const std::vector<ScanWindow> & InstrumentSettings::getScanWindows() const
  {
    return scan_windows_;
  }

  std::vector<ScanWindow> & InstrumentSettings::getScanWindows()
  {
    return scan_windows_;
  }

  void InstrumentSettings::setScanWindows(std::vector<ScanWindow>  scan_windows)
  {
    scan_windows_ =  scan_windows;
  }

  bool InstrumentSettings::getZoomScan() const
  {
    return zoom_scan_;
  }

  void InstrumentSettings::setZoomScan(bool zoom_scan)
  {
    zoom_scan_ = zoom_scan;
  }

}
