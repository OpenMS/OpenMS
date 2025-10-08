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

#include <QtCore/QVector>
#include <QtWidgets/QDialog>

namespace Ui
{
  class TOPPASIOMappingDialogTemplate;
}

namespace OpenMS
{
  class TOPPASEdge;

  /**
      @brief Dialog which allows to configure the input/output parameter mapping of an edge.

      This dialog allows to select an output parameter of the source vertex and an input
      parameter of the target vertex. Only valid selections are allowed, i.e. the type
      (file or list of files) and at least one valid file type of either vertex must match.

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASIOMappingDialog :
    public QDialog
  {
    Q_OBJECT

public:

    /// Constructor
    TOPPASIOMappingDialog(TOPPASEdge * parent);
    ~TOPPASIOMappingDialog() override;

public slots:

    /// Called instead of exec() after edge is constructed (in order to avoid showing the dialog if not necessary)
    int firstExec();

protected:

    /// Fills the table
    void fillComboBoxes_();

    /// The edge we are configuring
    TOPPASEdge * edge_;

protected slots:

    /// Called when OK is pressed; checks if the selected parameters are valid
    void checkValidity_();

private:
    Ui::TOPPASIOMappingDialogTemplate* ui_;
  };

}
