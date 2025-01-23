// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <ui_TOPPViewOpenDialog.h>

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>


// QT includes
#include <QtWidgets/QButtonGroup>

// STL includes
#include <iostream>

using namespace std;


namespace OpenMS
{

  TOPPViewOpenDialog::TOPPViewOpenDialog(const String & data_name, bool as_window, bool as_2d, bool cutoff, QWidget * parent) :
    QDialog(parent),
    map_as_2d_disabled_(false),
    ui_(new Ui::TOPPViewOpenDialogTemplate)
  {
    ui_->setupUi(this);

    //init map view
    if (!as_2d)
    {
      ui_->d1_->setChecked(true);
      ui_->d1_->setFocus();
    }
    else
    {
      ui_->d2_->setChecked(true);
      ui_->d2_->setFocus();
    }

    //init intensity cutoff
    if (cutoff)
    {
      ui_->intensity_cutoff_->setChecked(true);
    }

    //init open as
    if (!as_window)
    {
      ui_->layer_->setChecked(true);
      ui_->layer_->setFocus();
    }
    else
    {
      ui_->window_->setChecked(true);
      ui_->window_->setFocus();
    }
    connect(ui_->merge_combo_, SIGNAL(activated(int)), ui_->merge_, SLOT(click()));

    //set title
    setWindowTitle((String("Open data options for ") + data_name).toQString());
  }

  TOPPViewOpenDialog::~TOPPViewOpenDialog()
  {
    delete ui_;
  }

  bool TOPPViewOpenDialog::viewMapAs2D() const
  {
    return ui_->d2_->isChecked();
  }

  bool TOPPViewOpenDialog::viewMapAs1D() const
  {
    return ui_->d1_->isChecked();
  }

  bool TOPPViewOpenDialog::isCutoffEnabled() const
  {
    return ui_->intensity_cutoff_->isChecked();
  }

  bool TOPPViewOpenDialog::isDataDIA() const
  {
    return ui_->dia_data_->isChecked();
  }

  bool TOPPViewOpenDialog::openAsNewWindow() const
  {
    return ui_->window_->isChecked();
  }

  void TOPPViewOpenDialog::disableDimension(bool as_2d)
  {
    ui_->d1_->setChecked(!as_2d);
    ui_->d1_->setEnabled(false);
    ui_->d2_->setChecked(as_2d);
    ui_->d2_->setEnabled(false);
    ui_->d3_->setEnabled(false);
    map_as_2d_disabled_ = true;
  }

  void TOPPViewOpenDialog::disableCutoff(bool /* cutoff_on */)
  {
    ui_->intensity_cutoff_->setChecked(false);
  }

  void TOPPViewOpenDialog::disableLocation(bool as_window)
  {
    ui_->window_->setEnabled(false);
    ui_->layer_->setEnabled(false);
    ui_->merge_->setEnabled(false);
    ui_->merge_combo_->setEnabled(false);
    if (as_window)
    {
      ui_->window_->setChecked(true);
    }
    else
    {
      ui_->layer_->setChecked(true);
    }
  }

  void TOPPViewOpenDialog::updateViewMode_(QAbstractButton * button)
  {
    if (button == ui_->layer_ || button == ui_->merge_)
    {
      ui_->d1_->setEnabled(false);
      ui_->d2_->setEnabled(false);
      ui_->d3_->setEnabled(false);
    }
    else if (!map_as_2d_disabled_)
    {
      ui_->d1_->setEnabled(true);
      ui_->d2_->setEnabled(true);
      ui_->d3_->setEnabled(true);
    }
  }

  void TOPPViewOpenDialog::setMergeLayers(const std::map<Size, String> & layers)
  {
    // remove all items
    ui_->merge_combo_->clear();

    if (!layers.empty())
    {
      ui_->merge_->setEnabled(true);
      ui_->merge_combo_->setEnabled(true);
      UInt i = 0;
      for (const auto& it : layers)
      {
        ui_->merge_combo_->insertItem(i++, it.second.toQString(), (int)(it.first));
      }
    }
    else
    {
      ui_->merge_->setEnabled(false);
      ui_->merge_combo_->setEnabled(false);
    }
  }

  Int TOPPViewOpenDialog::getMergeLayer() const
  {
    if (ui_->merge_->isChecked())
    {
      return ui_->merge_combo_->itemData(ui_->merge_combo_->currentIndex()).toInt();
    }

    return -1;
  }

}
