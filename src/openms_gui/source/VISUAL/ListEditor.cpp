// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ListEditor.h>
#include <OpenMS/DATASTRUCTURES/String.h>


// for DIALOG
#include <QtWidgets/QPushButton>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QFileDialog>

#include <utility>
#include <vector>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    //////////////////////////////////////////////////////////////
    //DELEGATE
    //////////////////////////////////////////////////////////////
    ListEditorDelegate::ListEditorDelegate(QObject * parent) :
      QItemDelegate(parent)
    {
    }

    QWidget * ListEditorDelegate::createEditor(QWidget * parent, const QStyleOptionViewItem &, const QModelIndex & index) const
    {
      if (type_ == ListEditor::INPUT_FILE)
      {
        QLineEdit * editor = new QLineEdit(parent);
        //editor->setReadOnly(true);
        QString str = index.data(Qt::DisplayRole).toString();

        editor->setFocusPolicy(Qt::StrongFocus);

        file_name_ = QFileDialog::getOpenFileName(editor, tr("Input File List"), str);

        return editor;
      }
      else if (type_ == ListEditor::OUTPUT_FILE)
      {
        QLineEdit * editor = new QLineEdit(parent);
        //editor->setReadOnly(true);
        QString str = index.data(Qt::DisplayRole).toString();
        file_name_ = QFileDialog::getSaveFileName(editor, tr("Output File List"), str);

        return editor;
      }
      else if (type_ == ListEditor::STRING && !restrictions_.empty())
      {
        QComboBox * editor = new QComboBox(parent);
        QStringList list;
        list.append("");
        list += restrictions_.toQString().split(",");
        editor->addItems(list);
        return editor;
      }
      else
      {
        QLineEdit * editor = new QLineEdit(parent);
        editor->setFocusPolicy(Qt::StrongFocus);

        return editor;
      }
    }

    void ListEditorDelegate::setEditorData(QWidget * editor, const QModelIndex & index) const
    {
      if (index.isValid())
      {
        QString str = index.data(Qt::DisplayRole).toString();

        if (type_ == ListEditor::INPUT_FILE || type_ == ListEditor::OUTPUT_FILE)
        {
          if (!file_name_.isNull())
          {
            static_cast<QLineEdit *>(editor)->setText(file_name_);
          }
        }
        else if (qobject_cast<QComboBox *>(editor))
        {
          int index = static_cast<QComboBox *>(editor)->findText(str);
          if (index == -1)
          {
            index = 0;
          }
          static_cast<QComboBox *>(editor)->setCurrentIndex(index);
        }
        else if (qobject_cast<QLineEdit *>(editor))
        {
          static_cast<QLineEdit *>(editor)->setText(str);
        }
      }
    }

    void ListEditorDelegate::setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const
    {
      QVariant present_value = index.data(Qt::DisplayRole);
      QVariant new_value;
      if (index.column() == 0)
      {
        if (qobject_cast<QComboBox *>(editor))        //Drop-down list for enums
        {
          new_value = QVariant(static_cast<QComboBox *>(editor)->currentText());
        }
        else if (type_ == ListEditor::INPUT_FILE || type_ == ListEditor::OUTPUT_FILE)
        {
          new_value = QVariant(static_cast<QLineEdit *>(editor)->text());         //file_name_;
          file_name_ = "\0";
        }
        else
        {
          if (type_ == ListEditor::FLOAT && static_cast<QLineEdit *>(editor)->text() == "")
          {
            new_value = QVariant("0.0");
          }
          else if (type_ == ListEditor::INT && static_cast<QLineEdit *>(editor)->text() == "")
          {
            new_value = QVariant("0");
          }
          else
          {
            new_value = QVariant(static_cast<QLineEdit *>(editor)->text());
          }
        }
        //check if it matches the restrictions or is empty
        if (new_value.toString() != "")
        {
          bool restrictions_met = true;
          switch (type_)
          {
            //check if valid integer
            case ListEditor::INT:
            {
              bool ok;
              new_value.toString().toLong(&ok);
              if (!ok)
              {
                QMessageBox::warning(nullptr, "Invalid value", QString("Cannot convert '%1' to integer number!").arg(new_value.toString()));
                new_value = present_value;
                if (new_value == "")
                  new_value = 0;
              }

              //restrictions
              vector<String> parts;
              if (restrictions_.split(' ', parts))
              {
                if (!parts[0].empty() && new_value.toInt() < parts[0].toInt())
                {
                  restrictions_met = false;
                }
                if (!parts[1].empty() && new_value.toInt() > parts[1].toInt())
                {
                  restrictions_met = false;
                }
              }
            }
            break;

            case ListEditor::FLOAT:               //check if valid float
            {
              bool ok;
              new_value.toString().toDouble(&ok);
              if (!ok)
              {
                QMessageBox::warning(nullptr, "Invalid value", QString("Cannot convert '%1' to floating point number!").arg(new_value.toString()));
                new_value = present_value;
                if (new_value == "")
                {
                  new_value = 0;
                }
              }

              //restrictions
              vector<String> parts;
              if (restrictions_.split(' ', parts))
              {
                if (!parts[0].empty() && new_value.toDouble() < parts[0].toDouble())
                {
                  restrictions_met = false;
                }
                if (!parts[1].empty() && new_value.toDouble() > parts[1].toDouble())
                {
                  restrictions_met = false;
                }
              }
            }
            break;

          default:
            break;
          }

          if (!restrictions_met)
          {
            QMessageBox::warning(nullptr, "Invalid value", QString("Value restrictions not met: %1").arg(index.sibling(index.row(), 3).data(Qt::DisplayRole).toString()));
            new_value = present_value;
          }
        }
      }

      //check if modified
      if (new_value != present_value)
      {
        model->setData(index, new_value);
        model->setData(index, QBrush(Qt::yellow), Qt::BackgroundRole);
      }
    }

    void ListEditorDelegate::updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex &) const
    {
      editor->setGeometry(option.rect);
    }

    void ListEditorDelegate::setType(const ListEditor::Type type)
    {
      type_ = type;
    }

    void ListEditorDelegate::setRestrictions(const String & restrictions)
    {
      restrictions_ = restrictions;
    }

    void ListEditorDelegate::setTypeName(QString name)
    {
      typeName_ = std::move(name);
    }

    void ListEditorDelegate::setFileName(QString name)
    {
      file_name_ = std::move(name);
    }

    //////////////////////////////////////////////////////////////
    //LISTTABLE
    //////////////////////////////////////////////////////////////
    ListTable::ListTable(QWidget * parent) :
      QListWidget(parent)
    {
    }

    StringList ListTable::getList()
    {
      String stringit;
      list_.clear();
      for (Int i = 0; i < count(); ++i)
      {
        stringit = item(i)->text();
        if (!stringit.empty())
        {
          stringit.trim();
        }
        list_.push_back(stringit);
      }
      return list_;
    }

    void ListTable::setList(const StringList & list, ListEditor::Type type)
    {
      type_ = type;

      for (UInt i = 0; i < list.size(); ++i)
      {
        QListWidgetItem * item = new QListWidgetItem(list[i].toQString());
        item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);

        insertItem(i, item);
      }
      list_ = list;
      adjustSize();
    }

    void ListTable::createNewRow()
    {

      QListWidgetItem * item = nullptr;
      switch (type_)
      {
      case ListEditor::INT:
        item = new QListWidgetItem("0");
        break;

      case ListEditor::FLOAT:
        item = new QListWidgetItem("0.0");
        break;

      default:
        item = new QListWidgetItem("");
      }
      item->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
      item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable);
      addItem(item);
      item->setSelected(true);
      setCurrentRow(row(item));
      itemActivated(item);
      edit(currentIndex());
    }

    void ListTable::removeCurrentRow()
    {
      takeItem(currentRow());
    }

  }  //namespace Internal

  ////////////////////////////////////////////////////////////
  //ListEditor
  ////////////////////////////////////////////////////////////
  ListEditor::ListEditor(QWidget * parent, const QString& title) :
    QDialog(parent)
  {
    listTable_ = new Internal::ListTable(this);
    listTable_->setRowHidden(-1, true);
    listDelegate_ = new Internal::ListEditorDelegate(listTable_);
    listTable_->setItemDelegate(listDelegate_);

    removeRowButton_ = new QPushButton(tr("&delete"));
    newRowButton_ = new QPushButton(tr("&new"));
    newRowButton_->setDefault(true);
    OkButton_ = new QPushButton(tr("&ok"));
    CancelButton_ = new QPushButton(tr("&cancel"));

    connect(newRowButton_, SIGNAL(clicked()), listTable_, SLOT(createNewRow()));
    connect(removeRowButton_, SIGNAL(clicked()), listTable_, SLOT(removeCurrentRow()));

    QDialogButtonBox * rightLayout = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Vertical);
    rightLayout->addButton(newRowButton_, QDialogButtonBox::ActionRole);
    rightLayout->addButton(removeRowButton_, QDialogButtonBox::ActionRole);

    connect(rightLayout, SIGNAL(accepted()), this, SLOT(accept()));
    connect(rightLayout, SIGNAL(rejected()), this, SLOT(reject()));
    QHBoxLayout * mainLayout = new QHBoxLayout;
    mainLayout->addWidget(listTable_);
    mainLayout->addWidget(rightLayout);
    setLayout(mainLayout);
    QString tit =  "List Editor" + title;
    setWindowTitle(tit);
    setMinimumSize(800, 500);
  }

  StringList ListEditor::getList() const
  {
    return listTable_->getList();
  }

  void ListEditor::setList(const StringList & list, ListEditor::Type type)
  {
    type_ = type;
    listTable_->setList(list, type_);
    listDelegate_->setType(type_);
  }

  void ListEditor::setListRestrictions(const String & restrictions)
  {
    listDelegate_->setRestrictions(restrictions);
  }

  void ListEditor::setTypeName(QString name)
  {
    listDelegate_->setTypeName(std::move(name));
  }

} //namespace OpenMS
