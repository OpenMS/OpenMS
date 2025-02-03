// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtWidgets/QDialog>

namespace Ui
{
  class Plot1DPrefDialogTemplate;
}

namespace OpenMS
{
  namespace Internal
  {
    ///Preferences dialog for Plot1DWidget
    class OPENMS_GUI_DLLAPI Plot1DPrefDialog :
      public QDialog
    {
      Q_OBJECT

public:
      ///Constructor
      Plot1DPrefDialog(QWidget * parent);
      ~Plot1DPrefDialog() override;
private:
      Ui::Plot1DPrefDialogTemplate* ui_;
    };
  }
}
