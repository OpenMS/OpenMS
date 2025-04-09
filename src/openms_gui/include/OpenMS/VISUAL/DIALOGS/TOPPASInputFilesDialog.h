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
  class TOPPASInputFilesDialogTemplate;
}

namespace OpenMS
{
  namespace Internal
  {
    class InputFileList;
  }

  /**
      @brief Dialog which allows to specify a list of input files

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASInputFilesDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /// Constructor
    TOPPASInputFilesDialog(QWidget* parent)
     : TOPPASInputFilesDialog(QStringList(), "", parent) {}
    TOPPASInputFilesDialog(const QStringList& list, const QString& cwd, QWidget* parent = 0);
    ~TOPPASInputFilesDialog() override;

    void getFilenames(QStringList& files) const;

    const QString& getCWD() const;


private:
    Ui::TOPPASInputFilesDialogTemplate* ui_;
    Internal::InputFileList* ifl_;
  };
  
}

// this is required to allow Ui_SwathTabWidget (auto UIC'd from .ui) to have a TOPPASInputFilesDialog member
using TOPPASInputFilesDialog = OpenMS::TOPPASInputFilesDialog;

