// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/InstrumentSettings.h>

#include <utility>

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

  InstrumentSettings::~InstrumentSettings() = default;

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
    scan_windows_ =  std::move(scan_windows);
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

