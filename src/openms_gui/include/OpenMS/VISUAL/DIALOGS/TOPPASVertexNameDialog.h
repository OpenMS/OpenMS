// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtWidgets/QDialog>

namespace Ui
{
  class TOPPASVertexNameDialogTemplate;
}

namespace OpenMS
{
  /**
      @brief Dialog which allows to change the name of an input/output vertex

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASVertexNameDialog :
    public QDialog
  {
    Q_OBJECT

public:

    /// Constructor
    TOPPASVertexNameDialog(const QString& name, const QString& input_regex = QString());
    ~TOPPASVertexNameDialog() override;

    /// Returns the name
    QString getName();



  private:
    Ui::TOPPASVertexNameDialogTemplate* ui_;

  };

}
