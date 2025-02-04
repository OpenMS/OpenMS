// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/FeatureEditDialog.h>
#include <ui_FeatureEditDialog.h>

using namespace std;

namespace OpenMS
{

  FeatureEditDialog::FeatureEditDialog(QWidget * parent) :
    QDialog(parent),
    feature_(),
    ui_(new Ui::FeatureEditDialogTemplate)
  {
    ui_->setupUi(this);
  }

  FeatureEditDialog::~FeatureEditDialog()
  {
    delete ui_;
  }

  void FeatureEditDialog::setFeature(const Feature & feature)
  {
    //copy feature
    feature_ = feature;
    //update widgets according to feature data
    ui_->mz_->setValue(feature_.getMZ());
    ui_->rt_->setValue(feature_.getRT());
    ui_->int_->setValue(feature_.getIntensity());
    ui_->charge_->setValue(feature_.getCharge());
  }

  const Feature & FeatureEditDialog::getFeature() const
  {
    //update feature data according to widget
    feature_.setMZ(ui_->mz_->value());
    feature_.setRT(ui_->rt_->value());
    feature_.setIntensity(ui_->int_->value());
    feature_.setCharge(ui_->charge_->value());

    //return feature
    return feature_;
  }

} // namespace
