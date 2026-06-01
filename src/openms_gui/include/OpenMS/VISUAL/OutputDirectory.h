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
  class OutputDirectoryTemplate;
}

namespace OpenMS
{
  /**
      @brief A simple widget with a line-edit and a browse button to choose filenames

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI OutputDirectory :
    public QWidget
  {
    Q_OBJECT
  public:
    /// Constructor
    OutputDirectory(QWidget* parent);
    /// Destructor
    ~OutputDirectory();

    /// Sets the text in the line-edit
    void setDirectory(const QString& dir);

    /// return the directory currently set (does not need to be valid)
    QString getDirectory() const;
    
    /// check if the current directory exists and is writeable
    bool dirNameValid() const;

  signals:
    /// emitted whenever the outputdirectory is changed (also when setDirectory() is used)
    void directoryChanged(const QString& dir);

  public slots:

    /// Lets the user select the file via a file dialog
    void showFileDialog();
    
  private slots:
    /// forward internal textEdit::textChanged to directoryChanged signal
    void textEditChanged_(const QString& new_text);

  private:
    Ui::OutputDirectoryTemplate* ui_;
  };

}
