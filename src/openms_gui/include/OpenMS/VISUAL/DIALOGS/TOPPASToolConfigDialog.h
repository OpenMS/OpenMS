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

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

class QComboBox;
class QPushButton;
class QRadioButton;
class QString;

#include <QtWidgets/QDialog>

namespace OpenMS
{
  class ParamEditor;

  /**
      @brief TOPP tool configuration dialog

      In the dialog, the user can set the parameters for the tool

      This information can then be used to execute the tool.

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASToolConfigDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /**
        @brief Constructor

        @param parent Qt parent widget
        @param param The param we are editing
        @param default_dir The default directory for loading and storing INI files
        @param tool_name The name of the TOPP tool (used to invoke it on the commandline)
        @param tool_type The type of the tool ('-type' parameter of TOPP tool on the commandline). Leave empty if no type exists.
        @param tool_desc The tool description
        @param hidden_entries List of entries that are used already in edges etc and should not be shown
    */
    TOPPASToolConfigDialog(QWidget * parent, Param & param, const String& default_dir, const String& tool_name, const String& tool_type, const String& tool_desc, const QVector<String>& hidden_entries);
    ///Destructor
    ~TOPPASToolConfigDialog() override;

private:
    /// ParamEditor for reading ini-files
    ParamEditor * editor_;
    /// The param we are editing
    Param * param_;
    /// Param for loading the ini-file
    Param arg_param_;
    /// default-dir of ini-file to open
    String default_dir_;
    /// name of ini-file
    QString filename_;
    /// The name of the tool
    String tool_name_;
    /// The type of the tool
    String tool_type_;
    /// The parameters already explained by in edges
    QVector<String> hidden_entries_;

protected slots:
    /// Slot for OK button
    void ok_();
    /// loads an ini-file into the editor_
    void loadINI_();
    /// stores an ini-file from the editor_
    void storeINI_();
  };

}
