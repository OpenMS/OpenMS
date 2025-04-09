// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#ifndef Q_MOC_RUN
  #include <OpenMS/DATASTRUCTURES/ListUtils.h>
#endif

#include <QDialog>
#include <QListWidget>
#include <QItemDelegate>

class QPushButton;

namespace OpenMS
{
  namespace Internal
  {
    class ListTable;
    class ListEditorDelegate;
  }

  /**
      @brief Editor for editing int, double and string lists (including output and input file lists)
  */
  class OPENMS_GUI_DLLAPI ListEditor :
    public QDialog
  {
    Q_OBJECT

public:
    //types of lists
    enum Type
    {
      INT,
      FLOAT,
      STRING,
      OUTPUT_FILE,
      INPUT_FILE
    };

    ///Constructor
    ListEditor(QWidget * parent = nullptr, const QString& title = "");
    ///returns modified list
    StringList getList() const;
    ///sets list (and its type)that will be modified by user
    void setList(const StringList & list, ListEditor::Type type);
    ///set restrictions for list elements
    void setListRestrictions(const String & restrictions);
    ///set name of type
    void setTypeName(QString name);

private:
    ///List type
    ListEditor::Type type_;
    ///displays the list
    Internal::ListTable * listTable_;
    ///Delegate between view and model
    Internal::ListEditorDelegate * listDelegate_;
    /// button for new Row
    QPushButton * newRowButton_;
    ///button for removing row
    QPushButton * removeRowButton_;
    ///button clicked if modifications are accepted
    QPushButton * OkButton_;
    ///button clicked if modifications are rejected
    QPushButton * CancelButton_;
  };

  /**
      @brief Namespace used to hide implementation details from users.

  */
  namespace Internal
  {
    class OPENMS_GUI_DLLAPI ListTable :
      public QListWidget
    {
      Q_OBJECT

public:

      //Default Constructor
      ListTable(QWidget * parent = nullptr);

      //returns a list_
      StringList getList();

      //sets new list
      void setList(const StringList & list, ListEditor::Type type);

public slots:
      void createNewRow();
      void removeCurrentRow();

private:
      ///List type
      ListEditor::Type type_;
      //everything is internally stored as stringlist
      StringList list_;
    };

    /**
        @brief Internal delegate class

        This handles editing of items.
    */
    class OPENMS_GUI_DLLAPI ListEditorDelegate :
      public QItemDelegate
    {
      Q_OBJECT

public:
      ///Constructor
      ListEditorDelegate(QObject * parent);
      /// not reimplemented
      QWidget * createEditor(QWidget * parent, const QStyleOptionViewItem & option, const QModelIndex & index) const override;
      /// Sets the data to be displayed and edited by the editor for the item specified by index.
      void setEditorData(QWidget * editor, const QModelIndex & index) const override;
      /// Sets the data for the specified model and item index from that supplied by the editor. If data changed in a cell, that is if it is different from an initial value, then set its background color to yellow and emit the modified signal otherwise make it white
      void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const override;
      /// Updates the editor for the item specified by index according to the style option given.
      void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const override;

      //sets Type of List
      void setType(const ListEditor::Type type);
      //sets restrictions for listelements
      void setRestrictions(const String & restrictions);
      ///set name of type
      void setTypeName(QString name);
      ///sets the fileName
      void setFileName(QString name);

private:
      /// Not implemented => private
      ListEditorDelegate();
      ///List type
      ListEditor::Type type_;
      ///restrictions for list elements
      String restrictions_;
      ///type name. used to distinguish output/input from string lists
      QString typeName_;
      ///used to set input and output values in setModelData
      mutable QString file_name_;

    };
  }

} // namespace OpenMS
