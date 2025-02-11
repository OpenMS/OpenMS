// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/InstrumentSettingsVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  InstrumentSettingsVisualizer::InstrumentSettingsVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<InstrumentSettings>()
  {
    addLabel_("Modify the settings of the instrument.");
    addSeparator_();
    addComboBox_(instrumentsettings_scan_mode_, "Scan mode");
    addBooleanComboBox_(zoom_scan_, "Zoom scan");
    addComboBox_(instrumentsettings_polarity_, "Polarity");

    finishAdding_();
  }

  void InstrumentSettingsVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(instrumentsettings_scan_mode_, &temp_.NamesOfScanMode[temp_.getScanMode()], 1);
      fillComboBox_(instrumentsettings_polarity_, &IonSource::NamesOfPolarity[temp_.getPolarity()], 1);

    }
    else
    {
      fillComboBox_(instrumentsettings_scan_mode_, InstrumentSettings::NamesOfScanMode, InstrumentSettings::SIZE_OF_SCANMODE);
      fillComboBox_(instrumentsettings_polarity_, IonSource::NamesOfPolarity, IonSource::SIZE_OF_POLARITY);


      instrumentsettings_scan_mode_->setCurrentIndex(temp_.getScanMode());
      zoom_scan_->setCurrentIndex(temp_.getZoomScan());
      instrumentsettings_polarity_->setCurrentIndex(temp_.getPolarity());
    }
  }

  void InstrumentSettingsVisualizer::store()
  {
    ptr_->setScanMode((InstrumentSettings::ScanMode)instrumentsettings_scan_mode_->currentIndex());
    ptr_->setZoomScan(zoom_scan_->currentIndex());
    ptr_->setPolarity((IonSource::Polarity)instrumentsettings_polarity_->currentIndex());

    temp_ = (*ptr_);
  }

  void InstrumentSettingsVisualizer::undo_()
  {
    update_();
  }

}
