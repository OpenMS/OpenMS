// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Plot1DPrefDialog.h>
#include <ui_Plot1DPrefDialog.h>


using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    Plot1DPrefDialog::Plot1DPrefDialog(QWidget * parent) :
      QDialog(parent),
      ui_(new Ui::Plot1DPrefDialogTemplate)
    {
      ui_->setupUi(this);
    }

    Plot1DPrefDialog::~Plot1DPrefDialog()
    {
      delete ui_;
    }

  }   //namespace Internal
} //namspace OpenMS
