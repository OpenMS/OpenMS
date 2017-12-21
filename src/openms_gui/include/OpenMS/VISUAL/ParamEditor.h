// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_PARAMEDITOR_H
#define OPENMS_VISUAL_PARAMEDITOR_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/VISUAL/UIC/ui_ParamEditor.h>
#include <QtGui/QLineEdit>

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
      /// Checks if a @p name is valid for the entry corresponding to @p index (checks if it would be duplicate)
      bool exists_(QString name, QModelIndex index) const;

private slots:
      ///For closing ListEditor and updating ParamEditor
      void commitAndCloseListEditor_();
      ///For closing QcomboBox and updating ParamEditor
      void commitAndCloseComboBox_();
      ///if cancel in ListEditor is clicked Dialog is closed and changes are rejected
      void closeListEditor_();
      /// ...
      void commitAndCloseLineEdit_();
private:
      /// Not implemented
      ParamEditorDelegate();
      /// used to modify value of output and input files( not for output and input lists)
      mutable QString fileName_;
      /// true if a QLineEdit is still open and has not committed its data yet (so storing the current param is a bad idea)
      mutable bool has_uncommited_data_;
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
      bool edit(const QModelIndex & index, EditTrigger trigger, QEvent * event) override;

signals:
      ///Signal that is emitted when a new item is selected
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
    ParamEditor(QWidget * parent = nullptr);
    /// load method for Param object
    void load(Param & param);
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
