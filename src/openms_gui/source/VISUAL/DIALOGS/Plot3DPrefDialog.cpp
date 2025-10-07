// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Plot3DPrefDialog.h>
#include <ui_Plot3DPrefDialog.h>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    Plot3DPrefDialog::Plot3DPrefDialog(QWidget * parent) :
      QDialog(parent),
      ui_(new Ui::Plot3DPrefDialogTemplate)
    {
      ui_->setupUi(this);
    }

    Plot3DPrefDialog::~Plot3DPrefDialog()
    {
      delete ui_;
    }

  }   //namespace Internal
} //namspace OpenMS
