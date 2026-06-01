// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class InputFileTemplate;
}

namespace OpenMS
{
  /**
      @brief A simple widget with a line-edit and a browse button to choose filenames

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI InputFile :
    public QWidget
  {
    Q_OBJECT
  public:
    /// Constructor
    InputFile(QWidget* parent);
    /// Destructor
    ~InputFile();

    /// support drag'n'drop of files from OS window manager
    void dragEnterEvent(QDragEnterEvent* e) override;
    /// support drag'n'drop of files from OS window manager
    void dropEvent(QDropEvent* e) override;
    void dragMoveEvent(QDragMoveEvent* pEvent) override;

    /// Sets the text in the line-edit
    void setFilename(const QString& filename);

    /// Returns the filename currently set in the line-edit
    QString getFilename() const;

    /// Users can only choose certain filetypes, e.g. "Transition sqLite file (*.pqp)"
    void setFileFormatFilter(const QString& fff);

    /// get the CWD (according to most recently added file)
    const QString& getCWD() const;
    /// set the current working directory (for opening files). If the input is not empty, the cwd will not be altered, unless @p force is used
    void setCWD(const QString& cwd, bool force = false);
 
  signals:
    /// emitted when a new file is added (by drag'n'drop or 'Browse' button)
    void updatedCWD(QString new_cwd);
    /// emitted when a new file is added (by drag'n'drop or 'Browse' button)
    void updatedFile(QString new_path);

  public slots:

    /// Lets the user select the file via a file dialog
    void showFileDialog();


  protected:
    /// optional filter during file browsing
    QString file_format_filter_;
    /// the current working directory according to the last file added
    QString cwd_;

  private:
    Ui::InputFileTemplate* ui_;
  };

}
