// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Types.h>

#include <QtWidgets/QLineEdit>
#include <QtWidgets/QItemDelegate>
#include <QtWidgets/QTreeWidget>

class QModelIndex;
class QStyleOptionViewItem;
class QAbstractItemModel;
#include <QtCore/qcontainerfwd.h> // for QStringList
class QString;

namespace Ui
{
  class ParamEditorTemplate;
}

namespace OpenMS
{
  class String;
  class Param;
  class ParamEditor;
  /**
      @brief Namespace used to hide implementation details from users.

  */
  namespace Internal
  {

    /**
      @brief Custom QLineEdit which emits a signal when losing focus (such that we can commit its data)

     */
    class OPENMS_GUI_DLLAPI OpenMSLineEdit
      : public QLineEdit
    {
      Q_OBJECT
  public:
      OpenMSLineEdit(QWidget * w)
        :QLineEdit(w)
      {}

signals:
      /// emitted on focusOutEvent
      void lostFocus();


    protected:
        void   focusOutEvent ( QFocusEvent * e ) override;
        void   focusInEvent ( QFocusEvent * e ) override;
    };
    /**
        @brief Internal delegate class for QTreeWidget

        This handles editing of items.
    */
    class OPENMS_GUI_DLLAPI ParamEditorDelegate :
      public QItemDelegate
    {
      Q_OBJECT

public:
      ///Constructor
      ParamEditorDelegate(QObject * parent);
      /// Returns the widget(combobox or QLineEdit) used to edit the item specified by index for editing. Prevents edit operations on nodes' values and types
      QWidget * createEditor(QWidget * parent, const QStyleOptionViewItem & option, const QModelIndex & index) const override;
      /// Sets the data to be displayed and edited by the editor for the item specified by index.
      void setEditorData(QWidget * editor, const QModelIndex & index) const override;
      /// Sets the data for the specified model and item index from that supplied by the editor. If data changed in a cell, that is if it is different from an initial value, then set its background color to yellow and emit the modified signal otherwise make it white
      void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const override;
      /// Updates the editor for the item specified by index according to the style option given.
      void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const override;

      /// true if the underlying tree has an open QLineEdit which has uncommitted data
      bool hasUncommittedData() const;
signals:
      /// signal for showing ParamEditor if the Model data changed
      void modified(bool) const;

protected:
      /// a shortcut to calling commit(), which calls setModelData(); useful for embedded editors, but not for QDialogs etc
      bool eventFilter(QObject* editor, QEvent* event) override;

private slots:
      ///For closing any editor and updating ParamEditor
      void commitAndCloseEditor_();
      ///if cancel in any editor is clicked, the Dialog is closed and changes are rejected
      void closeEditor_();
      /// ... a bit special, because reset uncommited data
      void commitAndCloseLineEdit_();

private:
      /// Not implemented
      ParamEditorDelegate();
      /// used to modify value of output and input files( not for output and input lists)
      mutable QString fileName_;
      /// holds a directory name (for output directories)
      mutable QString dirName_;
      /// true if a QLineEdit is still open and has not committed its data yet (so storing the current param is a bad idea)
      mutable bool has_uncommited_data_;
    };

    /// QTreeWidget that emits a signal whenever a new row is selected
    class OPENMS_GUI_DLLAPI ParamTree :
      public QTreeWidget
    {
      Q_OBJECT

public:
      /// Constructor
      ParamTree(QWidget * parent);
      /// Overloaded edit method to activate F2 use
      bool edit(const QModelIndex & index, EditTrigger trigger, QEvent * event) override;

signals:
      /// Signal that is emitted when a new item is selected
      void selected(const QModelIndex & index);

protected slots:
      /// Reimplemented virtual slot
      void selectionChanged(const QItemSelection & selected, const QItemSelection &) override;
    };

  }

  
  /**
    @brief A GUI for editing or viewing a Param object

    It supports two display modes:
    - normal mode: only the main parameters are displayed, advanced parameters are hidden.
    - advanced mode: all parameters are displayed (only available when advanced parameters are provided)

    @image html ParamEditor.png

    @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI ParamEditor :
    public QWidget
  {
    Q_OBJECT

public:
    /// Role of the entry
    enum
    {
      NODE,                     ///< Section
      NORMAL_ITEM,              ///< Item that is always shown
      ADVANCED_ITEM             ///< Item that is shown only in advanced mode
    };

    /// constructor
    ParamEditor(QWidget* parent = nullptr);
    /// destructor
    ~ParamEditor() override;
    
    /// load method for Param object
    void load(Param& param);
    /// store edited data in Param object
    void store();
    /// Indicates if the data changed since last save
    bool isModified() const;
    /// Clears all parameters
    void clear();

public slots:
    /// Notifies the widget that the content was changed.
    /// Emits the modified(bool) signal if the state changed.
    void setModified(bool is_modified);

signals:
    /// item was edited
    void modified(bool);

protected slots:
    /// Switches between normal and advanced mode
    void toggleAdvancedMode(bool advanced);
    /// Shows the documentation of an item in doc_
    void showDocumentation(const QModelIndex & index);

protected:
    /// recursive helper method for method storeRecursive()
    void storeRecursive_(QTreeWidgetItem * child, String path, std::map<String, String> & section_descriptions);

    /// Pointer to the tree widget
    Internal::ParamTree* tree_;
    /// The data to edit
    Param* param_;
    /// Indicates that the data was modified since last store/load operation
    bool modified_;
    /// Indicates if normal mode or advanced mode is activated
    bool advanced_mode_;

private:
    Ui::ParamEditorTemplate* ui_;
  };


} // namespace OpenMS

