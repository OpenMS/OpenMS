// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_PARAMEDITOR_H
#define OPENMS_VISUAL_PARAMEDITOR_H

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/VISUAL/UIC/ui_ParamEditor.h>


#include <QtGui/QItemDelegate>
#include <QtGui/QTreeWidget>

class QModelIndex;
class QStyleOptionViewItem;
class QAbstractItemModel;
class QStringList;
class QString;

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
      QWidget * createEditor(QWidget * parent, const QStyleOptionViewItem & option, const QModelIndex & index) const;
      /// Sets the data to be displayed and edited by the editor for the item specified by index.
      void setEditorData(QWidget * editor, const QModelIndex & index) const;
      /// Sets the data for the specified model and item index from that supplied by the editor. If data changed in a cell, that is if it is different from an initial value, then set its background color to yellow and emit the modified signal otherwise make it white
      void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const;
      /// Updates the editor for the item specified by index according to the style option given.
      void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const;

signals:
      /// signal for showing ParamEditor if the Model data changed
      void modified(bool) const;

protected:
      /// Checks if a @p name is valid for the entry corresponding to @p index (checks if it would be duplicate)
      bool exists_(QString name, QModelIndex index) const;

private slots:
      ///For closing ListEditor and updating ParamEditor
      void commitAndCloseListEditor_();
      ///For closing QcomboBox and updating ParamEditor
      void commitAndCloseComboBox_();
      ///if cancel in ListEditor is clicked Dialog is closed and changes are rejected
      void closeListEditor_();
private:
      /// Not implemented
      ParamEditorDelegate();
      ///used to modify value of output and input files( not for output and input lists)
      mutable QString fileName_;
    };

    /// QTreeWidget that emits a signal whenever a new row is selected
    class OPENMS_GUI_DLLAPI ParamTree :
      public QTreeWidget
    {
      Q_OBJECT

public:
      ///Constructor
      ParamTree(QWidget * parent);
      /// Overloaded edit method to activate F2 use
      bool edit(const QModelIndex & index, EditTrigger trigger, QEvent * event);

signals:
      ///Signal that is emitted when a new item is selected
      void selected(const QModelIndex & index);

protected slots:
      /// Reimplemented virtual slot
      void selectionChanged(const QItemSelection & selected, const QItemSelection &);
    };

  }

  /**
      @brief A GUI for editing or viewing a Param object

      It supports two display modes:
      - normal mode: only the main parameters are displayed, advanced parameters are hidden.
      - advanced mode: all parameters are displayed.

      @image html ParamEditor.png

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI ParamEditor :
    public QWidget,
    public Ui::ParamEditorTemplate
  {
    Q_OBJECT

public:
    /// Role of the entry
    enum
    {
      NODE,                           ///< Section
      NORMAL_ITEM,              ///< Item that is always shown
      ADVANCED_ITEM             ///< Item that is shown only in advanced mode
    };

    /// constructor
    ParamEditor(QWidget * parent = 0);
    /// load method for Param object
    void load(Param & param);
    /// store edited data in Param object
    void store();
    /// Indicates if the data changed since last save
    bool isModified() const;
    /// Clears all parameters
    void clear();

signals:
    /// item was edited
    void modified(bool);

protected slots:
    /// Notifies the widget that the content was changed.
    /// Emits the modified(bool) signal if the state changed.
    void setModified(bool is_modified);
    /// Switches between normal and advanced mode
    void toggleAdvancedMode(bool advanced);
    /// Shows the documentation of an item in doc_
    void showDocumentation(const QModelIndex & index);

protected:
    /// recursive helper method for method storeRecursive()
    void storeRecursive_(QTreeWidgetItem * child, String path, std::map<String, String> & section_descriptions);

    /// Pointer to the tree widget
    Internal::ParamTree * tree_;
    /// The data to edit
    Param * param_;
    /// Indicates that the data was modified since last store/load operation
    bool modified_;
    /// Indicates if normal mode or advanced mode is activated
    bool advanced_mode_;
  };


} // namespace OpenMS

#endif // OPENMS_VISUAL_PARAMEDITOR_H
