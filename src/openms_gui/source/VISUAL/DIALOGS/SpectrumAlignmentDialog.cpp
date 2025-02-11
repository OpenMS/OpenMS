// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/SpectrumAlignmentDialog.h>
#include <ui_SpectrumAlignmentDialog.h>


#include <OpenMS/VISUAL/Plot1DWidget.h>

// QT includes
#include <QtWidgets/QButtonGroup>

#include <vector>

namespace OpenMS
{
  SpectrumAlignmentDialog::SpectrumAlignmentDialog(Plot1DWidget * parent) :
    layer_indices_1_(),
    layer_indices_2_(),
    ui_(new Ui::SpectrumAlignmentDialogTemplate)
  {
    ui_->setupUi(this);

    QButtonGroup * button_group = new QButtonGroup(this);
    button_group->addButton(ui_->ppm);
    button_group->addButton(ui_->da);
    ui_->da->setChecked(true);

    Plot1DCanvas * cc = parent->canvas();
    for (UInt i = 0; i < cc->getLayerCount(); ++i)
    {
      const auto& layer = cc->getLayer(i);
      if (layer.flipped)
      {
        ui_->layer_list_2->addItem(layer.getName().toQString());
        layer_indices_2_.push_back(i);
      }
      else
      {
        ui_->layer_list_1->addItem(layer.getName().toQString());
        layer_indices_1_.push_back(i);
      }
    }
    // select first item of each list
    if (ui_->layer_list_1->count() > 0)
    {
      ui_->layer_list_1->setCurrentRow(0);
    }
    if (ui_->layer_list_2->count() > 0)
    {
      ui_->layer_list_2->setCurrentRow(0);
    }
  }

  SpectrumAlignmentDialog::~SpectrumAlignmentDialog()
  {
    delete ui_;
  }

  double SpectrumAlignmentDialog::getTolerance() const
  {
    return ui_->tolerance_spinbox->value();
  }

  bool SpectrumAlignmentDialog::isPPM() const
  {
    return ui_->ppm->isChecked();
  }


  Int SpectrumAlignmentDialog::get1stLayerIndex()
  {
    if (ui_->layer_list_1->count() == 0 || ui_->layer_list_1->currentRow() == -1)
    {
      return -1;
    }
    if (layer_indices_1_.size() > (Size)(ui_->layer_list_1->currentRow()))
    {
      return layer_indices_1_[(Size)(ui_->layer_list_1->currentRow())];
    }
    return -1;
  }

  Int SpectrumAlignmentDialog::get2ndLayerIndex()
  {
    if (ui_->layer_list_2->count() == 0 || ui_->layer_list_2->currentRow() == -1)
    {
      return -1;
    }
    if (layer_indices_2_.size() > (Size)(ui_->layer_list_2->currentRow()))
    {
      return layer_indices_2_[(Size)(ui_->layer_list_2->currentRow())];
    }
    return -1;
  }

} // namespace
