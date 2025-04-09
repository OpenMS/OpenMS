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
  class TOPPASInputFileDialogTemplate;
}

namespace OpenMS
{
  class InputFile;
  /**
      @brief Dialog which allows to specify an input file

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASInputFileDialog :
    public QDialog
  {
    Q_OBJECT
  public:
    /// Constructor
    TOPPASInputFileDialog(const QString& file_name);
    ~TOPPASInputFileDialog();

    /// users can only choose certain filetypes
    void setFileFormatFilter(const QString& fff);

    /// Returns the filename
    QString getFilename() const;

protected slots:

    /// Called when OK is pressed; checks if the selected file is valid
    void checkValidity_();

private:
    Ui::TOPPASInputFileDialogTemplate* ui_;
  };
  
  
} // ns OpenMS

// this is required to allow Ui_TOPPASInputFileDialog (auto UIC'd from .ui) to have a InputFile member
using InputFile = OpenMS::InputFile;

