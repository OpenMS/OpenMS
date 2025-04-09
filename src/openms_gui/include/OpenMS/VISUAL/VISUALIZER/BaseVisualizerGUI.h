// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Types.h>

#include <QtWidgets>

class QPushButton;
class QGridLayout;
class QLineEdit;
class QTextEdit;
class QComboBox;
class QListWidget;
class QLabel;

namespace OpenMS
{
  /**
      @brief A base class for all visualizer classes

      This class provides the GUI for part for all visualizers.

      Several additional members are provided by the BaseVisualizer class.
      The two classes cannot be merged, as templates and the Qt meta object compiler cannot be combined.

      Visualizers are mainly used by the MetaDataBrowser.
  */
  class OPENMS_GUI_DLLAPI BaseVisualizerGUI :
    public QWidget
  {
    Q_OBJECT

public:

    ///Constructor
    BaseVisualizerGUI(bool editable = false, QWidget * parent = nullptr);

    /// Returns if the values are editable
    bool isEditable() const;

signals:

    /// Sends a status message, if date is not in proper format.
    void sendStatus(std::string status);

public slots:

    ///Saves the changes made in the GUI to the object.
    virtual void store() = 0;

protected:

    /// Adds a label to the grid layout.
    void addLabel_(const QString& label);
    /// Adds a label to a certain row
    void addLabel_(const QString& label, UInt row);
    /// Adds a line edit field with label to the grid layout.
    void addLineEdit_(QLineEdit * & ptr, QString label);
    /// Adds a line edit field to the grid layout including a int validator
    void addIntLineEdit_(QLineEdit * & ptr, QString label);
    /// Adds a line edit field to the grid layout including a double validator
    void addDoubleLineEdit_(QLineEdit * & ptr, QString label);
    /// Adds a line edit field with label and button to the next free position in the grid.
    void addLineEditButton_(const QString& label, QLineEdit * & ptr1, QPushButton * & ptr2, const QString& buttonlabel);
    /// Adds a list edit field to the grid layout.
    void addListView_(QListWidget * & ptr, QString label);
    /// Adds a text edit field to the grid layout.
    void addTextEdit_(QTextEdit * & ptr, QString label);
    /// Adds a drop-down field to the grid layout.
    void addComboBox_(QComboBox * & ptr, QString label);
    /// Adds a boolean drop-down field to the grid layout ( 'true'=1,  'false'=0 ).
    void addBooleanComboBox_(QComboBox * & ptr, QString label);
    /// Fills a combo box with string @p items (the number of strings is determined by @p item_count).
    void fillComboBox_(QComboBox * & ptr, const std::string * items, int item_count);
    /// Adds vertical spacer.
    void addVSpacer_();
    /// Adds a button to the next free position in the grid.
    void addButton_(QPushButton * & ptr, const QString& label);
    /// Adds two buttons in a row.
    void add2Buttons_(QPushButton * & ptr1, const QString& label1, QPushButton * & ptr2, const QString& label2);
    /// Adds a horizontal line as a separator.
    void addSeparator_();
    /// Adds buttons common to all visualizers
    void finishAdding_();

    ///Undo button
    QPushButton * undo_button_;
    /// The main layout.
    QGridLayout * mainlayout_;
    /// Counter for the current grid row.
    UInt row_;
    /// Edit flag
    bool editable_;
  };


}
