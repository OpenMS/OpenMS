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
  class Plot3DPrefDialogTemplate;
}

namespace OpenMS
{
  namespace Internal
  {
    ///Preferences dialog for Plot3DWidget
    class OPENMS_GUI_DLLAPI Plot3DPrefDialog :
      public QDialog
    {
      Q_OBJECT

public:
      ///Constructor
      Plot3DPrefDialog(QWidget * parent);
      ~Plot3DPrefDialog() override;
    private:
      Ui::Plot3DPrefDialogTemplate* ui_;
    };
  }
}
