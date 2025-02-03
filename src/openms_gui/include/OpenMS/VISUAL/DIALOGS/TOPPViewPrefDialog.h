// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>


#include <OpenMS/DATASTRUCTURES/Param.h>

#include <QtWidgets/QDialog>

namespace Ui
{
  class TOPPViewPrefDialogTemplate;
}

namespace OpenMS
{
  namespace Internal
  {
    /**
        @brief Preferences dialog for TOPPView

        @ingroup TOPPView_elements
    */
    class OPENMS_GUI_DLLAPI TOPPViewPrefDialog :
      public QDialog
    {
      Q_OBJECT

public:
      TOPPViewPrefDialog(QWidget * parent);
      ~TOPPViewPrefDialog() override;

      /// initialize GUI values with these parameters
      void setParam(const Param& param);

      /// update the parameters given the current GUI state.
      /// Can be used to obtain default parameters and their names.
      Param getParam() const;

protected slots:
      void browseDefaultPath_();
      void browsePluginsPath_();
private:
      Ui::TOPPViewPrefDialogTemplate* ui_;
      mutable Param param_; ///< is updated in getParam()
      Param tsg_param_; ///< params for TheoreticalSpectrumGenerator in the TSG tab
    };
  }
}
