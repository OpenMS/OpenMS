// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtWidgets/QMdiArea>
#include <QtWidgets/QMainWindow>

class QToolBar;
class QAction;
class QString;
class QFileDialog;

namespace OpenMS
{
  /**
      @brief shows the ParamEditor widget in a QMainWindow with a toolbar
  */
  class OPENMS_GUI_DLLAPI INIFileEditorWindow :
    public QMainWindow
  {
    Q_OBJECT

public:
    /// menu is created here
    INIFileEditorWindow(QWidget * parent = nullptr);
    /// when user closes window a message box asks the user if he wants to save
    void closeEvent(QCloseEvent * event) override;

public slots:
    ///loads the xml-file into a Param object and loads Param into ParamEditor
    bool openFile(const String & filename = "");
    /// saves the users changes in a xml-file if the Param object is valid
    bool saveFile();
    /// like saveFile but with a file dialog to choose a filename
    bool saveFileAs();
    /// if the user changes data in ParamEditor the title shows a '*'
    void updateWindowTitle(bool);

private:
    /// ParamEditor object for visualization
    ParamEditor * editor_;
    /// Param object for storing data
    Param param_;
    /// filename of xml-file to store the Param object
    QString filename_;
    /// path used as next default location of the load/store dialogs
    String current_path_;
  };
}

