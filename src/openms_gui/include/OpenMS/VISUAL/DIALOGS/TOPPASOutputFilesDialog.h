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
  class TOPPASOutputFilesDialogTemplate;
}

namespace OpenMS
{
  class OutputDirectory;

  /**
      @brief Dialog which allows to specify the directory for the output files

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASOutputFilesDialog :
    public QDialog
  {
    Q_OBJECT

public:

    /// Constructor
    TOPPASOutputFilesDialog(const QString& dir_name, int num_jobs);
    ~TOPPASOutputFilesDialog() override;

    /// Returns the name of the directory
    QString getDirectory() const;

    /// Returns the maximum number of jobs in the spinbox
    int getNumJobs() const;

public slots:

    /// Lets the user select the directory via a file dialog
    void showFileDialog();

protected slots:

    /// Called when OK is pressed; checks if the selected file is valid
    void checkValidity_();
private:
    Ui::TOPPASOutputFilesDialogTemplate* ui_;
  };

}

using OutputDirectory = OpenMS::OutputDirectory;
