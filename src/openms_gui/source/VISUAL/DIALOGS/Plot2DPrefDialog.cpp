// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Plot2DPrefDialog.h>
#include <ui_Plot2DPrefDialog.h>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    Plot2DPrefDialog::Plot2DPrefDialog(QWidget * parent) :
      QDialog(parent),
      ui_(new Ui::Plot2DPrefDialogTemplate)
    {
      ui_->setupUi(this);
    }

    Plot2DPrefDialog::~Plot2DPrefDialog()
    {
      delete ui_;
    }

  }   //namespace Internal
} //namspace OpenMS
