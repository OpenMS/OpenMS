// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/MassAnalyzerVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  MassAnalyzerVisualizer::MassAnalyzerVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<MassAnalyzer>()
  {
    addLabel_("Modify massanalyzer information.");
    addSeparator_();

    addIntLineEdit_(order_, "Order");
    addComboBox_(type_, "Type");
    addComboBox_(res_method_, "Resolution method");
    addComboBox_(res_type_, "Resolution type");
    addComboBox_(scan_dir_, "Scan direction");
    addComboBox_(scan_law_, "Scan law");
    addComboBox_(reflectron_state_, "Reflectron state");

    addDoubleLineEdit_(res_, "Resolution");
    addDoubleLineEdit_(acc_, "Accuracy");
    addDoubleLineEdit_(scan_rate_, "Scan rate (in s)");
    addDoubleLineEdit_(scan_time_, "Scan time (in s)");
    addDoubleLineEdit_(TOF_, "TOF Total path length (in meter)");
    addDoubleLineEdit_(iso_, "Isolation width (in m/z)");
    addDoubleLineEdit_(final_MS_, "Final MS exponent");
    addDoubleLineEdit_(magnetic_fs_, "Magnetic field strength (in T)");

    finishAdding_();
  }

  void MassAnalyzerVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(type_, &temp_.NamesOfAnalyzerType[static_cast<size_t>(temp_.getType())], 1);
      fillComboBox_(res_method_, &temp_.NamesOfResolutionMethod[static_cast<size_t>(temp_.getResolutionMethod())], 1);
      fillComboBox_(res_type_, &temp_.NamesOfResolutionType[static_cast<size_t>(temp_.getResolutionType())], 1);
      fillComboBox_(scan_dir_, &temp_.NamesOfScanDirection[static_cast<size_t>(temp_.getScanDirection())], 1);
      fillComboBox_(scan_law_, &temp_.NamesOfScanLaw[static_cast<size_t>(temp_.getScanLaw())], 1);
      fillComboBox_(reflectron_state_, &temp_.NamesOfReflectronState[static_cast<size_t>(temp_.getReflectronState())], 1);
    }
    else
    {
      fillComboBox_(type_, temp_.NamesOfAnalyzerType, static_cast<int>(MassAnalyzer::AnalyzerType::SIZE_OF_ANALYZERTYPE));
      fillComboBox_(res_method_, temp_.NamesOfResolutionMethod, static_cast<int>(MassAnalyzer::ResolutionMethod::SIZE_OF_RESOLUTIONMETHOD));
      fillComboBox_(res_type_, temp_.NamesOfResolutionType, static_cast<int>(MassAnalyzer::ResolutionType::SIZE_OF_RESOLUTIONTYPE));
      fillComboBox_(scan_dir_, temp_.NamesOfScanDirection, static_cast<int>(MassAnalyzer::ScanDirection::SIZE_OF_SCANDIRECTION));
      fillComboBox_(scan_law_, temp_.NamesOfScanLaw, static_cast<int>(MassAnalyzer::ScanLaw::SIZE_OF_SCANLAW));
      fillComboBox_(reflectron_state_, temp_.NamesOfReflectronState, static_cast<int>(MassAnalyzer::ReflectronState::SIZE_OF_REFLECTRONSTATE));

      type_->setCurrentIndex(static_cast<int>(temp_.getType()));
      res_method_->setCurrentIndex(static_cast<int>(temp_.getResolutionMethod()));
      res_type_->setCurrentIndex(static_cast<int>(temp_.getResolutionType()));
      scan_dir_->setCurrentIndex(static_cast<int>(temp_.getScanDirection()));
      scan_law_->setCurrentIndex(static_cast<int>(temp_.getScanLaw()));
      reflectron_state_->setCurrentIndex(static_cast<int>(temp_.getReflectronState()));
    }

    order_->setText(StringUtils::toStr(temp_.getOrder()).c_str());
    res_->setText(StringUtils::toStr(temp_.getResolution()).c_str());
    acc_->setText(StringUtils::toStr(temp_.getAccuracy()).c_str());
    scan_rate_->setText(StringUtils::toStr(temp_.getScanRate()).c_str());
    scan_time_->setText(StringUtils::toStr(temp_.getScanTime()).c_str());
    TOF_->setText(StringUtils::toStr(temp_.getTOFTotalPathLength()).c_str());
    iso_->setText(StringUtils::toStr(temp_.getIsolationWidth()).c_str());
    final_MS_->setText(StringUtils::toStr(temp_.getFinalMSExponent()).c_str());
    magnetic_fs_->setText(StringUtils::toStr(temp_.getMagneticFieldStrength()).c_str());
  }

  void MassAnalyzerVisualizer::store()
  {
    ptr_->setOrder(order_->text().toInt());
    ptr_->setType((MassAnalyzer::AnalyzerType)type_->currentIndex());
    ptr_->setResolutionMethod((MassAnalyzer::ResolutionMethod)res_method_->currentIndex());
    ptr_->setResolutionType((MassAnalyzer::ResolutionType)res_type_->currentIndex());
    ptr_->setScanDirection((MassAnalyzer::ScanDirection)scan_dir_->currentIndex());
    ptr_->setScanLaw((MassAnalyzer::ScanLaw)scan_law_->currentIndex());
    ptr_->setReflectronState((MassAnalyzer::ReflectronState)reflectron_state_->currentIndex());

    ptr_->setResolution(res_->text().toDouble());
    ptr_->setAccuracy(acc_->text().toDouble());
    ptr_->setScanRate(scan_rate_->text().toDouble());
    ptr_->setScanTime(scan_time_->text().toDouble());
    ptr_->setTOFTotalPathLength(TOF_->text().toDouble());
    ptr_->setIsolationWidth(iso_->text().toDouble());
    ptr_->setFinalMSExponent(final_MS_->text().toInt());
    ptr_->setMagneticFieldStrength(magnetic_fs_->text().toDouble());

    temp_ = (*ptr_);
  }

  void MassAnalyzerVisualizer::undo_()
  {
    update_();
  }

}
