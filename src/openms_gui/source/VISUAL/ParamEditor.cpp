// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/VISUAL/ParamEditor.h>
#include <ui_ParamEditor.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/VISUAL/ListEditor.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtWidgets/QMessageBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QShortcut>
#include <QtWidgets/QMenu>
#include <QItemSelection>
#include <QtCore/QStringList>
#include <QtWidgets/QLabel>
#include <QtWidgets/QFileDialog>

#include <stack>
#include <limits>
#include <sstream>

using namespace std;


/*

Description of the data stored in the items:

            | Column 0  | Column 1    | Column 2 | Column 3         |
---------------------------------------------------------------------
DisplayRole | name      | value       | type     | restr. (display) |
UserRole    | NODE/ITEM | description | restr.   |                  |


*/

namespace OpenMS
{
  namespace Internal
  {
    void OpenMSLineEdit::focusInEvent ( QFocusEvent * /* e */)
    {
      //std::cerr << "got focus";
    }

    void OpenMSLineEdit::focusOutEvent ( QFocusEvent * /* e */ )
    {
      //std::cerr << "lost focus";
      emit lostFocus();
    }

    ParamEditorDelegate::ParamEditorDelegate(QObject* parent) :
      QItemDelegate(parent)
    {
    }

    QWidget * ParamEditorDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem&, const QModelIndex& index) const
    {
      Int type = index.sibling(index.row(), 0).data(Qt::UserRole).toInt();

      // only create editor for first column (value column)
      if (index.column() == 1 && type != ParamEditor::NODE)
      {
        has_uncommited_data_ = false; // by default all data is committed

        QString dtype = index.sibling(index.row(), 2).data(Qt::DisplayRole).toString();
        QString restrictions = index.sibling(index.row(), 2).data(Qt::UserRole).toString();
        QString value = index.sibling(index.row(), 1).data(Qt::DisplayRole).toString();
        if (dtype == "string" && restrictions != "")     //Drop-down list for enums
        {
          QComboBox * editor = new QComboBox(parent);
          QStringList list;
          list.append("");
          list += restrictions.split(",");
          editor->addItems(list);
          connect(editor, SIGNAL(activated(int)), this, SLOT(commitAndCloseComboBox_()));
          return editor;
        }
        else if (dtype == "output file")
        {
          QLineEdit * editor = new QLineEdit(parent);
          QString dir = "";        // = index.sibling(index.row(),0).data(Qt::DisplayRole).toString();
          if (File::isDirectory(value) || File::writable(value))
          {
            dir = File::absolutePath(value).toQString();
          }
          fileName_ = QFileDialog::getSaveFileName(editor, tr("Output File"), dir);
          return editor;
        }
        else if (dtype == "input file")
        {
          QLineEdit * editor = new QLineEdit(parent);
          QString dir = "";        // = index.sibling(index.row(),0).data(Qt::DisplayRole).toString();
          if (File::isDirectory(value) || File::exists(value))
          {
            dir = File::absolutePath(value).toQString();
          }
          fileName_ = QFileDialog::getOpenFileName(editor, tr("Input File"), dir);
          return editor;
        }
        else if (dtype == "string list" || dtype == "int list" || dtype == "double list" || dtype == "input file list" || dtype == "output file list")   // for lists
        {
          QString str = "<" + index.sibling(index.row(), 0).data(Qt::DisplayRole).toString() + "> " + "(<" + dtype + ">)";
          ListEditor * editor = new ListEditor(nullptr, str);
          editor->setTypeName(index.sibling(index.row(), 0).data(Qt::DisplayRole).toString());
          editor->setModal(true);
          connect(editor, SIGNAL(accepted()), this, SLOT(commitAndCloseListEditor_()));
          connect(editor, SIGNAL(rejected()), this, SLOT(closeListEditor_()));
          return editor;
        }
        else         //LineEditor for rest
        {
          OpenMSLineEdit * editor = new OpenMSLineEdit(parent);
          editor->setFocusPolicy(Qt::StrongFocus);
          connect(editor, SIGNAL(lostFocus()), this, SLOT(commitAndCloseLineEdit_()));
          //std::cerr << "created ... \n";
          has_uncommited_data_ = true;
          return editor;
        }
      }
      return nullptr;
    }

    void ParamEditorDelegate::setEditorData(QWidget * editor, const QModelIndex & index) const
    {
      QString str = index.data(Qt::DisplayRole).toString();

      // only set editor data for first column (value column)
      if (index.column() == 1)
      {
        if (qobject_cast<QComboBox *>(editor))       //Drop-down list for enums
        {
          int index = static_cast<QComboBox *>(editor)->findText(str);
          if (index == -1)
          {
            index = 0;
          }
          static_cast<QComboBox *>(editor)->setCurrentIndex(index);
        }
        else if (qobject_cast<QLineEdit *>(editor))      // LineEdit for other values
        {
          QString dtype = index.sibling(index.row(), 2).data(Qt::DisplayRole).toString();
          if (dtype == "output file" || dtype == "input file")          /// for output/input file
          {
            if (!fileName_.isNull())
            {
              static_cast<QLineEdit *>(editor)->setText(fileName_);
            }
          }
          else
          {
            if (str == "" && (dtype == "int" || dtype == "float"))
            {
              if (dtype == "int")
                static_cast<QLineEdit *>(editor)->setText("0");
              else if (dtype == "float")
                static_cast<QLineEdit *>(editor)->setText("nan");
            }
            else
            {
              static_cast<QLineEdit *>(editor)->setText(str);
            }
          }
        }
        else                // ListEditor for lists
        {
          String list;
          list = str.mid(1, str.length() - 2);
          QString type = index.sibling(index.row(), 2).data(Qt::DisplayRole).toString();
          StringList rlist = ListUtils::create<String>(list);
          for (UInt i = 0; i < rlist.size(); ++i)
          {
            rlist[i]  = rlist[i].trim();
          }
          if (type == "int list")
          {
            static_cast<ListEditor *>(editor)->setList(rlist, ListEditor::INT);
          }
          else if (type == "double list")
          {
            static_cast<ListEditor *>(editor)->setList(rlist, ListEditor::FLOAT);
          }
          else if (type == "string list")
          {
            static_cast<ListEditor *>(editor)->setList(rlist, ListEditor::STRING);
          }
          else if (type == "input file list")
          {
            static_cast<ListEditor *>(editor)->setList(rlist, ListEditor::INPUT_FILE);
          }
          else if (type == "output file list")
          {
            static_cast<ListEditor *>(editor)->setList(rlist, ListEditor::OUTPUT_FILE);
          }
          static_cast<ListEditor *>(editor)->setListRestrictions(index.sibling(index.row(), 2).data(Qt::UserRole).toString());
        }
      }
    }

    void ParamEditorDelegate::setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const
    {
      QVariant present_value = index.data(Qt::DisplayRole);
      QVariant new_value;
      StringList list;
      bool new_list = false;
      // only set model data for first column (value column)
      if (index.column() == 1)
      {
        //extract new value
        if (qobject_cast<QComboBox *>(editor))       //Drop-down list for enums
        {
          new_value = QVariant(static_cast<QComboBox *>(editor)->currentText());
        }
        else if (qobject_cast<QLineEdit *>(editor))
        {
          QString dtype = index.sibling(index.row(), 2).data(Qt::DisplayRole).toString();
          if (dtype == "output file" || dtype == "input file")        // input/outut file
          {

            new_value = QVariant(static_cast<QLineEdit *>(editor)->text());
            fileName_ = "\0";
          }
          else if (static_cast<QLineEdit *>(editor)->text() == "" && ((dtype == "int") || (dtype == "float")))         //numeric
          {
            if (dtype == "int")
              new_value = QVariant("0");
            else if (dtype == "float")
              new_value = QVariant("nan");
          }
          else
          {
            new_value = QVariant(static_cast<QLineEdit *>(editor)->text());
          }
        }
        else
        {
          list = static_cast<ListEditor *>(editor)->getList();
          for (UInt i = 1; i < list.size(); ++i)
          {
            list[i] = "\n" + list[i];
          }
          new_list = true;
        }
        //check if it matches the restrictions or is empty
        if (new_value.toString() != "")
        {
          QString type = index.sibling(index.row(), 2).data(Qt::DisplayRole).toString();
          bool restrictions_met = true;
          String restrictions = index.sibling(index.row(), 2).data(Qt::UserRole).toString();
          if (type == "int")         //check if valid integer
          {
            bool ok(true);
            new_value.toString().toLong(&ok);
            if (!ok)
            {
              QMessageBox::warning(nullptr, "Invalid value", QString("Cannot convert '%1' to integer number!").arg(new_value.toString()));
              new_value = present_value;
            }
            //restrictions
            vector<String> parts;
            if (restrictions.split(' ', parts))
            {
              if (parts[0] != "" && new_value.toInt() < parts[0].toInt())
              {
                restrictions_met = false;
              }
              if (parts[1] != "" && new_value.toInt() > parts[1].toInt())
              {
                restrictions_met = false;
              }
            }
          }
          else if (type == "float")         //check if valid float
          {
            bool ok(true);
            new_value.toString().toDouble(&ok);
            if (!ok)
            {
              QMessageBox::warning(nullptr, "Invalid value", QString("Cannot convert '%1' to floating point number!").arg(new_value.toString()));
              new_value = present_value;
            }
            //restrictions
            vector<String> parts;
            if (restrictions.split(' ', parts))
            {
              if (parts[0] != "" && new_value.toDouble() < parts[0].toDouble())
              {
                restrictions_met = false;
              }
              if (parts[1] != "" && new_value.toDouble() > parts[1].toDouble())
              {
                restrictions_met = false;
              }
            }
          }
          if (!restrictions_met)
          {
            QMessageBox::warning(nullptr, "Invalid value", QString("Value restrictions not met: %1").arg(index.sibling(index.row(), 3).data(Qt::DisplayRole).toString()));
            new_value = present_value;
          }
        }
      }
      if (new_list)
      {
        stringstream ss;
        ss << list;
        QVariant new_value;
        new_value = QVariant(QString::fromStdString(ss.str()));
        model->setData(index, new_value);
        model->setData(index, QBrush(Qt::yellow), Qt::BackgroundRole);
        emit modified(true);
      }
      else
      {
        //check if modified
        if (new_value != present_value)
        {
          model->setData(index, new_value);
          model->setData(index, QBrush(Qt::yellow), Qt::BackgroundRole);
          emit modified(true);
        }
      }
    }

    void ParamEditorDelegate::updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex &) const
    {
      editor->setGeometry(option.rect);
    }

    bool ParamEditorDelegate::exists_(QString name, QModelIndex index) const
    {
      UInt current_index = 0;
      while (index.parent().child(current_index, 0).isValid())
      {
        if (
          current_index != (UInt)(index.row())
           &&
          index.parent().child(current_index, 0).data(Qt::DisplayRole).toString() == name
           &&
          (
            (index.data(Qt::UserRole).toInt() == 0 && index.parent().child(current_index, 0).data(Qt::UserRole).toInt() == 0)
          ||
            (index.data(Qt::UserRole).toInt() != 0 && index.parent().child(current_index, 0).data(Qt::UserRole).toInt() != 0)
          )
          )
        {
          return true;
        }
        ++current_index;
      }
      return false;
    }

    void ParamEditorDelegate::commitAndCloseListEditor_()
    {
      ListEditor * editor = qobject_cast<ListEditor *>(sender());
      emit commitData(editor);
      emit closeEditor(editor);
    }
    void ParamEditorDelegate::commitAndCloseLineEdit_()
    {
      has_uncommited_data_ = false;
      OpenMSLineEdit * editor = qobject_cast<OpenMSLineEdit *>(sender());
      emit commitData(editor);
      emit closeEditor(editor);
    }

    void ParamEditorDelegate::commitAndCloseComboBox_()
    {
      QComboBox * editor = qobject_cast<QComboBox *>(sender());
      emit commitData(editor);
      emit closeEditor(editor);
    }

    void ParamEditorDelegate::closeListEditor_()
    {
      ListEditor * editor = qobject_cast<ListEditor *>(sender());
      emit closeEditor(editor);
    }

    bool ParamEditorDelegate::hasUncommittedData() const
    {
      return has_uncommited_data_;
    }

    ///////////////////ParamTree/////////////////////////////////

    ParamTree::ParamTree(QWidget * parent) :
      QTreeWidget(parent)
    {
    }

    void ParamTree::selectionChanged(const QItemSelection & s, const QItemSelection &)
    {
      if (!s.empty())
      {
        emit selected(s.indexes().first());
      }
    }

    bool ParamTree::edit(const QModelIndex & index, EditTrigger trigger, QEvent * event)
    { // allow F2 or double click on any column in the current row
      if (trigger == QAbstractItemView::EditKeyPressed || trigger == QAbstractItemView::DoubleClicked)
      { // --> re-route to actual value column
        return QAbstractItemView::edit(index.sibling(index.row(), 1), trigger, event);
      }
      return QAbstractItemView::edit(index, trigger, event);
    }

  }

  ///////////////////ParamEditor/////////////////////////////////

  ParamEditor::ParamEditor(QWidget * parent) :
    QWidget(parent),
    param_(nullptr),
    modified_(false),
    advanced_mode_(false),
    ui_(new Ui::ParamEditorTemplate)
  {
    ui_->setupUi(this);
    tree_ = new Internal::ParamTree(this);
    //tree_->setMinimumSize(450, 200);
    tree_->setAllColumnsShowFocus(true);
    tree_->setColumnCount(4);
    tree_->setHeaderLabels(QStringList() << "parameter" << "value" << "type" << "restrictions");
    dynamic_cast<QVBoxLayout *>(layout())->insertWidget(0, tree_, 1);
    tree_->setItemDelegate(new Internal::ParamEditorDelegate(tree_));       // the delegate from above is set
    connect(tree_->itemDelegate(), SIGNAL(modified(bool)), this, SLOT(setModified(bool)));
    connect(ui_->advanced_, SIGNAL(toggled(bool)), this, SLOT(toggleAdvancedMode(bool)));
    connect(tree_, SIGNAL(selected(const QModelIndex &)), this, SLOT(showDocumentation(const QModelIndex &)));
  }


  ParamEditor::~ParamEditor()
  {
    delete ui_;
  }

  void ParamEditor::showDocumentation(const QModelIndex & index)
  {
    ui_->doc_->setText(index.sibling(index.row(), 1).data(Qt::UserRole).toString());
  }

  void ParamEditor::load(Param & param)
  {
    param_ = &param;

    tree_->clear();

    QTreeWidgetItem * parent = tree_->invisibleRootItem();
    QTreeWidgetItem * item = nullptr;

    bool has_advanced_item = false; // will be true if @p param has any advanced items; if still false, we disable the 'show advanced checkbox'

    for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
    {
      //********handle opened/closed nodes********
      const std::vector<Param::ParamIterator::TraceInfo> & trace = it.getTrace();
      for (std::vector<Param::ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
      {
        if (it2->opened)         //opened node
        {
          item = new QTreeWidgetItem(parent);
          //name
          item->setText(0, it2->name.toQString());
          item->setTextColor(0, Qt::darkGray);  // color of nodes with children

          //description
          item->setData(1, Qt::UserRole, it2->description.toQString());
          //role
          item->setData(0, Qt::UserRole, NODE);
          //flags
          if (param_ != nullptr)
          {
            item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
          }
          else
          {
            item->setFlags(Qt::ItemIsEnabled);
          }
          parent = item;
        }
        else         //closed node
        {
          parent = parent->parent();
          if (parent == nullptr)
            parent = tree_->invisibleRootItem();
        }
      }

      //********handle item********
      item = new QTreeWidgetItem(parent);

      // grey out non-editable columns (leaf nodes)
      bool is_required = it->tags.find("required") != it->tags.end();
      if (is_required)  // special color for required parameters
      {
        item->setTextColor(0, QColor(255, 140, 0, 255)); // orange
        item->setTextColor(2, QColor(255, 140, 0, 255));
        item->setTextColor(3, QColor(255, 140, 0, 255));
      }
      else
      {
        item->setTextColor(0, Qt::darkGray);
        item->setTextColor(2, Qt::darkGray);
        item->setTextColor(3, Qt::darkGray);
      }

      if (it->tags.count("advanced"))
      {
        item->setData(0, Qt::UserRole, ADVANCED_ITEM);
        has_advanced_item = true;
      }
      else       //advanced parameter
      {
        item->setData(0, Qt::UserRole, NORMAL_ITEM);
      }
      //name
      item->setText(0, it->name.toQString());
      //value
      if (it->value.valueType() == DataValue::STRING_LIST)
      {
        StringList string_list = it->value;
        String list_string = String("[") + ListUtils::concatenate(string_list, ",\n") + "]";
        item->setText(1, list_string.toQString());
      }
      else if (it->value.valueType() == DataValue::INT_LIST)
      {
        IntList list = it->value;
        StringList string_list;
        for (Size i = 0; i < list.size(); ++i)
        {
          string_list.push_back(list[i]);
        }
        String list_string = String("[") + ListUtils::concatenate(string_list, ",\n") + "]";
        item->setText(1, list_string.toQString());
      }
      else if (it->value.valueType() == DataValue::DOUBLE_LIST)
      {
        DoubleList list = it->value;
        StringList string_list;
        for (Size i = 0; i < list.size(); ++i)
        {
          string_list.push_back(list[i]);
        }
        String list_string = String("[") + ListUtils::concatenate(string_list, ",\n") + "]";
        item->setText(1, list_string.toQString());
      }
      else
      {
        item->setText(1, String(it->value).toQString());
      }
      //type
      switch (it->value.valueType())
      {
      case DataValue::INT_VALUE:
        item->setText(2, "int");
        break;

      case DataValue::DOUBLE_VALUE:
        item->setText(2, "float");
        break;

      case DataValue::STRING_VALUE:
        if (it->tags.count("input file"))
        {
          item->setText(2, "input file");
        }
        else if (it->tags.count("output file"))
        {
          item->setText(2, "output file");
        }
        else
        {
          item->setText(2, "string");
        }
        break;

      case DataValue::STRING_LIST:
        if (it->tags.count("input file"))
        {
          item->setText(2, "input file list");
        }
        else if (it->tags.count("output file"))
        {
          item->setText(2, "output file list");
        }
        else
        {
          item->setText(2, "string list");
        }
        break;

      case DataValue::INT_LIST:
        item->setText(2, "int list");
        break;

      case DataValue::DOUBLE_LIST:
        item->setText(2, "double list");
        break;

      default:
        break;
      }
      //restrictions (displayed and internal for easier parsing)
      switch (it->value.valueType())
      {
      case DataValue::INT_VALUE:
      case DataValue::INT_LIST:
      {
        String drest = "", irest = "";
        bool min_set = (it->min_int != -numeric_limits<Int>::max());
        bool max_set = (it->max_int != numeric_limits<Int>::max());
        if (max_set || min_set)
        {
          if (min_set)
          {
            drest += String("min: ") + it->min_int;
            irest += it->min_int;
          }
          irest += " ";
          if (max_set)
          {
            if (min_set && max_set)
              drest += " ";
            drest += String("max: ") + it->max_int;
            irest += it->max_int;
          }
          item->setText(3, drest.toQString());
        }
        item->setData(2, Qt::UserRole, irest.toQString());
      }
      break;

      case DataValue::DOUBLE_VALUE:
      case DataValue::DOUBLE_LIST:
      {
        String drest = "", irest = "";
        bool min_set = (it->min_float != -numeric_limits<double>::max());
        bool max_set = (it->max_float != numeric_limits<double>::max());
        if (max_set || min_set)
        {
          if (min_set)
          {
            drest += String("min: ") + it->min_float;
            irest += it->min_float;
          }
          irest += " ";
          if (max_set)
          {
            if (min_set && max_set)
              drest += " ";
            drest += String("max: ") + it->max_float;
            irest += it->max_float;
          }
          item->setText(3, drest.toQString());
        }
        item->setData(2, Qt::UserRole, irest.toQString());
      }
      break;

      case DataValue::STRING_VALUE:
      case DataValue::STRING_LIST:
      {
        String irest = "";
        if (it->valid_strings.size() != 0)
        {
          irest.concatenate(it->valid_strings.begin(), it->valid_strings.end(), ",");
        }
        if (irest != "")
        {
          String r_text = irest;
          if (r_text.size() > 255) // truncate restriction text, as some QT versions (4.6 & 4.7) will crash if text is too long
          {
            r_text = irest.prefix(251) + "...";
          }
          item->setText(3, r_text.toQString());
        }
        item->setData(2, Qt::UserRole, irest.toQString());
      }
      break;

      default:
        break;
      }

      //description
      item->setData(1, Qt::UserRole, it->description.toQString());
      //flags
      if (param_ != nullptr)
      {
        item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
      }
      else
      {
        item->setFlags(Qt::ItemIsEnabled);
      }
    }

    if (!has_advanced_item)
    {
      ui_->advanced_->setVisible(false);
    }

    tree_->expandAll();
    toggleAdvancedMode(advanced_mode_);

    tree_->resizeColumnToContents(0);
    tree_->resizeColumnToContents(1);
    tree_->resizeColumnToContents(2);
    tree_->resizeColumnToContents(3);
  }

  void ParamEditor::store()
  {
    //std::cerr << "store entered ...\n";

    // store only if no line-edit is opened (in which case data is uncommitted and will not be saved)
    // this applies only to INIFileEditor, where pressing Ctrl-s results in saving the current (but outdated) param
    if (param_ != nullptr &&
        !static_cast<Internal::ParamEditorDelegate*>(this->tree_->itemDelegate())->hasUncommittedData())
    {
      //std::cerr << "and done!...\n";
      QTreeWidgetItem * parent = tree_->invisibleRootItem();
      //param_->clear();

      for (Int i = 0; i < parent->childCount(); ++i)
      {
        map<String, String> section_descriptions;
        storeRecursive_(parent->child(i), "", section_descriptions);        //whole tree recursively
      }

      setModified(false);
    }
    //else std::cerr << "store aborted!\n";

  }

  void ParamEditor::clear()
  {
    tree_->clear();
  }

  void ParamEditor::storeRecursive_(QTreeWidgetItem * child, String path, map<String, String> & section_descriptions)
  {
    /**

        @todo: why would we "recreate" (setting restrictions etc) the the param object from scratch?
         updating everything that changed seems the better option, as
                     this is more robust against additions to Param

    */
    child->setData(1, Qt::BackgroundRole, QBrush(Qt::white));

    if (path == "")
    {
      path = child->text(0).toStdString();
    }
    else
    {
      path += String(":") + String(child->text(0).toStdString());
    }

    String description = child->data(1, Qt::UserRole).toString();

    StringList tag_list;
    try // might throw ElementNotFound
    {
      tag_list = param_->getTags(path);
    }
    catch (...)
    {
    }

    if (child->text(2) == "")  // node
    {
      if (description != "")
      {
        section_descriptions.insert(make_pair(path, description));
      }
    }
    else     //item + section descriptions
    {
      if (child->text(2) == "float")
      {
        param_->setValue(path, child->text(1).toDouble(), description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        vector<String> parts;
        if (restrictions.split(' ', parts))
        {
          if (parts[0] != "")
          {
            param_->setMinFloat(path, parts[0].toDouble());
          }
          if (parts[1] != "")
          {
            param_->setMaxFloat(path, parts[1].toDouble());
          }
        }
      }
      else if (child->text(2) == "string")
      {
        param_->setValue(path, child->text(1).toStdString(), description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        if (restrictions != "")
        {
          std::vector<String> parts;
          restrictions.split(',', parts);
          param_->setValidStrings(path, parts);
        }
      }
      else if (child->text(2) == "input file")
      {
        param_->setValue(path, child->text(1).toStdString(), description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        if (restrictions != "")
        {
          std::vector<String> parts;
          restrictions.split(',', parts);
          param_->setValidStrings(path, parts);
        }
      }
      else if (child->text(2) == "output file")
      {
        param_->setValue(path, child->text(1).toStdString(), description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        if (restrictions != "")
        {
          std::vector<String> parts;
          restrictions.split(',', parts);
          param_->setValidStrings(path, parts);
        }
      }
      else if (child->text(2) == "int")
      {
        param_->setValue(path, child->text(1).toInt(), description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        vector<String> parts;
        if (restrictions.split(' ', parts))
        {
          if (parts[0] != "")
          {
            param_->setMinInt(path, parts[0].toInt());
          }
          if (parts[1] != "")
          {
            param_->setMaxInt(path, parts[1].toInt());
          }
        }
      }
      String list;
      list = child->text(1).mid(1, child->text(1).length() - 2);
      StringList rlist = ListUtils::create<String>(list);
      for (UInt i = 0; i < rlist.size(); ++i)
      {
        rlist[i] = rlist[i].trim();
      }
      if (child->text(2) == "string list")
      {
        param_->setValue(path, rlist, description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        if (restrictions != "")
        {
          vector<String> parts;
          restrictions.split(',', parts);
          param_->setValidStrings(path, parts);
        }
      }
      else if (child->text(2) == "input file list")
      {
        param_->setValue(path, rlist, description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        if (restrictions != "")
        {
          vector<String> parts;
          restrictions.split(',', parts);
          param_->setValidStrings(path, parts);
        }
      }
      else if (child->text(2) == "output file list")
      {
        param_->setValue(path, rlist, description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        if (restrictions != "")
        {
          vector<String> parts;
          restrictions.split(',', parts);
          param_->setValidStrings(path, parts);
        }
      }
      else if (child->text(2) == "double list")
      {
        param_->setValue(path, ListUtils::create<double>(rlist), description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        vector<String> parts;
        if (restrictions.split(' ', parts))
        {
          if (parts[0] != "")
          {
            param_->setMinFloat(path, parts[0].toFloat());
          }
          if (parts[1] != "")
          {
            param_->setMaxFloat(path, parts[1].toFloat());
          }
        }
      }
      else if (child->text(2) == "int list")
      {
        param_->setValue(path, ListUtils::create<Int>(rlist), description, tag_list);
        String restrictions = child->data(2, Qt::UserRole).toString();
        vector<String> parts;
        if (restrictions.split(' ', parts))
        {
          if (parts[0] != "")
          {
            param_->setMinInt(path, parts[0].toInt());
          }
          if (parts[1] != "")
          {
            param_->setMaxInt(path, parts[1].toInt());
          }
        }
      }

      // set description node description if the prefix matches
      for (map<String, String>::const_iterator it = section_descriptions.begin(); it != section_descriptions.end(); ++it)
      {
        if (path.hasPrefix(it->first))
        {
          param_->setSectionDescription(it->first, it->second);
        }
      }
      section_descriptions.clear();
    }

    for (Int i = 0; i < child->childCount(); ++i)
    {
      storeRecursive_(child->child(i), path, section_descriptions);     //whole tree recursively
    }
  }

  void ParamEditor::setModified(bool is_modified)
  {
    if (is_modified != modified_)
    {
      modified_ = is_modified;
      emit modified(modified_);
    }
  }

  bool ParamEditor::isModified() const
  {
    return modified_;
  }

  void ParamEditor::toggleAdvancedMode(bool advanced)
  {
    advanced_mode_ = advanced;

    stack<QTreeWidgetItem *> stack, node_stack;

    //show/hide items
    stack.push(tree_->invisibleRootItem());
    while (!stack.empty())
    {
      QTreeWidgetItem * current = stack.top();
      stack.pop();

      Int type = current->data(0, Qt::UserRole).toInt();
      if (type != NODE)     //ITEM
      {
        if (advanced_mode_ && type == ADVANCED_ITEM)       //advanced mode
        {
          current->setHidden(false);
        }
        else if (!advanced_mode_ && type == ADVANCED_ITEM)       //Normal mode
        {
          current->setHidden(true);
        }
      }
      else       //NODE
      {
        for (Int i = 0; i < current->childCount(); ++i)
        {
          stack.push(current->child(i));
        }

        if (advanced_mode_)
        {
          current->setHidden(false);           //show all nodes in advanced mode
        }
        else
        {
          node_stack.push(current);           //store node pointers in normal mode
        }
      }
    }

    //hide sections that have no visible items in normal mode
    while (!node_stack.empty())
    {
      QTreeWidgetItem * current = node_stack.top();
      node_stack.pop();

      bool has_visible_children = false;
      for (Int i = 0; i < current->childCount(); ++i)
      {
        if (!current->child(i)->isHidden())
        {
          has_visible_children = true;
          break;
        }
      }
      if (!has_visible_children)
      {
        current->setHidden(true);
      }
    }

    //resize columns
    tree_->resizeColumnToContents(0);
    tree_->resizeColumnToContents(1);
    tree_->resizeColumnToContents(2);
    tree_->resizeColumnToContents(3);
  }

} // namespace OpenMS
