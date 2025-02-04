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
  class Plot2DPrefDialogTemplate;
}

namespace OpenMS
{
  namespace Internal
  {
    ///Preferences dialog for Plot2DWidget
    class OPENMS_GUI_DLLAPI Plot2DPrefDialog :
      public QDialog
    {
      Q_OBJECT

public:
      ///Constructor
      Plot2DPrefDialog(QWidget * parent);
      ~Plot2DPrefDialog() override;
private:
      Ui::Plot2DPrefDialogTemplate* ui_;
    };
  }
}
