// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TableView.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Qt5Port.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QFile>
#include <QFileDialog>
#include <QHeaderView>
#include <QMenu>
#include <QTextStream>

#include <iostream>

using namespace std;

///@improvement write the visibility-status of the columns in toppview.ini and read at start


namespace OpenMS
{
  TableView::TableView(QWidget* parent) :
    QTableWidget(parent)
  {
    this->setObjectName("table_widget");

    this->setSortingEnabled(true);

    this->setEditTriggers(QAbstractItemView::NoEditTriggers);
    this->setSelectionBehavior(QAbstractItemView::SelectRows);
    this->setShowGrid(false);

    this->setSelectionMode(QAbstractItemView::SingleSelection);

    this->horizontalHeader()->setSectionsMovable(true);
    this->horizontalHeader()->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this->horizontalHeader(), &QHeaderView::customContextMenuRequested, this, &TableView::headerContextMenu_);

    this->verticalHeader()->setHidden(true); // hide vertical column
    {
      QTableWidgetItem* proto_item = new QTableWidgetItem();
      proto_item->setTextAlignment(Qt::AlignCenter);
      this->setItemPrototype(proto_item);
    }
  }

  void TableView::headerContextMenu_(const QPoint& pos)
  {
    // create menu
    QMenu context_menu(this);

    // add actions which show/hide columns
    for (int i = 0; i != columnCount(); ++i)
    {
      QTableWidgetItem* ti = horizontalHeaderItem(i);
      if (ti == nullptr)
      {
        continue;
      }
      QAction* action = context_menu.addAction(ti->text(), [=]() {
        // invert visibility upon clicking the item
        setColumnHidden(i, !isColumnHidden(i));
        });
      action->setCheckable(true);
      action->setChecked(!isColumnHidden(i));
    }
    context_menu.exec(mapToGlobal(pos));
  }

  void TableView::setMandatoryExportColumns(QStringList& cols)
  {
    mandatory_export_columns_ = cols;
  }

  void TableView::exportEntries()
  {
    QString filename = QFileDialog::getSaveFileName(this, "Save File", "", "tsv file (*.tsv)");
    QFile f(filename);

    if (!f.open(QIODevice::WriteOnly))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(filename));
    }
    QTextStream ts(&f);
    QStringList str_list;
    
    QStringList cols_to_export = (getHeaderNames(WidgetHeader::VISIBLE_ONLY, true) + mandatory_export_columns_);    
    cols_to_export.removeDuplicates();

    QStringList all_header_names = getHeaderNames(WidgetHeader::WITH_INVISIBLE, true);

    // write header
    bool first{true};
    for (int c = 0; c < columnCount(); ++c)
    {
      // columns marked for export
      if (cols_to_export.indexOf(all_header_names[c]) != -1)
      {
        if (!first) 
        { 
          ts << "\t"; 
        }
        else 
        {
          first = false;
        }
        ts << all_header_names[c];        
      }
    }
    ts << "\n";

    // write entries
    for (int r = 0; r < rowCount(); ++r)
    {
      for (int c = 0; c < columnCount(); ++c)
      {
        // only export columns we marked for export
        if (cols_to_export.indexOf(all_header_names[c]) == -1)
        {
          continue;
        }

        QTableWidgetItem* ti = this->item(r, c);
        if (ti == nullptr)
        {
          str_list << "";
          std::cerr << "Warning: Empty table cell found at position: ["<< r << ' ' << c << "]\n";
        }
        else
        {
          if (ti->data(Qt::UserRole).isValid())
          {
            str_list << ti->data(Qt::UserRole).toString();
          }
          else if (ti->data(Qt::CheckStateRole).isValid()) // Note: item with check box also has a display role, so this test needs to come first
          {
            str_list << ti->data(Qt::CheckStateRole).toString();
          }
          else if (ti->data(Qt::DisplayRole).isValid())
          {
            str_list << ti->data(Qt::DisplayRole).toString();
          }
          else
          {
            str_list << "";
            std::cerr << "Warning: table cell with unhandled role found at position: [" << r << ' ' << c << "]\n";
          }
        }
      }
      ts << str_list.join("\t") + "\n";
      str_list.clear();
    }
    f.close();
  }

  void TableView::setHeaders(const QStringList& headers)
  {
    setColumnCount(headers.size());
    setHorizontalHeaderLabels(headers);
  }

  void TableView::hideColumns(const QStringList& header_names)
  {
    auto hset = toQSet(header_names);
    // add actions which show/hide columns
    for (int i = 0; i != columnCount(); ++i)
    {
      QTableWidgetItem* ti = horizontalHeaderItem(i);
      if (ti == nullptr)
      {
        continue;
      }
      if (hset.contains(ti->text()))
      {
        setColumnHidden(i, true);
        hset.remove(ti->text());
      }
    }
    if (!hset.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "header_names contains a column name which is unknown: " + String(hset.values().join(", ")));
    }
  }

  void TableView::appendRow()
  {
    insertRow(rowCount());
  }

  QTableWidgetItem* TableView::setAtBottomRow(const QString& text, size_t column_index, const QColor& background, const QColor& foreground)
  {
    QTableWidgetItem* item = itemPrototype()->clone();
    item->setText(text);
    return setAtBottomRow(item, column_index, background, foreground);
  }

  QTableWidgetItem* TableView::setAtBottomRow(const char* text, size_t column_index, const QColor& background, const QColor& foreground)
  {
    QTableWidgetItem* item = itemPrototype()->clone();
    item->setText(text);
    return setAtBottomRow(item, column_index, background, foreground);
  }

  QTableWidgetItem* TableView::setAtBottomRow(const int i, size_t column_index, const QColor& background, const QColor& foreground)
  {
    QTableWidgetItem* item = itemPrototype()->clone();
    item->setData(Qt::DisplayRole, i);
    return setAtBottomRow(item, column_index, background, foreground);
  }

  QTableWidgetItem* TableView::setAtBottomRow(const double d, size_t column_index, const QColor& background, const QColor& foreground)
  {
    QTableWidgetItem* item = itemPrototype()->clone();
    item->setData(Qt::DisplayRole, d);
    return setAtBottomRow(item, column_index, background, foreground);
  }

  QTableWidgetItem* TableView::setAtBottomRow(const bool selected, size_t column_index, const QColor& background, const QColor& foreground)
  {
    QTableWidgetItem* item = itemPrototype()->clone();
    item->setCheckState(selected ? Qt::Checked : Qt::Unchecked);
    /// sorting of columns is done by the DisplayRole, not the checkstate. So we need different content.
    updateCheckBoxItem(item);
    return setAtBottomRow(item, column_index, background, foreground);
  }

  QTableWidgetItem* TableView::setAtBottomRow(QTableWidgetItem* item, size_t column_index, const QColor& background, const QColor& foreground)
  {
    item->setBackground(QBrush(background));
    if (foreground.isValid())
    {
      item->setForeground(QBrush(foreground));
    }
    setItem(rowCount() - 1, (int)column_index, item);
    return item;
  }

  void TableView::updateCheckBoxItem(QTableWidgetItem* item)
  {
    // check if this function is called on checkbox items only (either no DisplayRole set or the text is '' or ' ')
    if (!item->data(Qt::DisplayRole).isValid() || 
        (item->data(Qt::DisplayRole).type() == QVariant::Type::String
          && (item->data(Qt::DisplayRole).toString().isEmpty() || item->data(Qt::DisplayRole).toString() == " ")
        )
       )
    {
      item->setText(item->checkState() == Qt::Checked ? " " : "");
    }
    else
    { 
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Called on non-checkbox item");
    }
  }
 
  QStringList TableView::getHeaderNames(const WidgetHeader which, bool use_export_name)
  {
    QStringList header_labels;
    for (int i = 0; i != columnCount(); ++i)
    {
      // do not export hidden columns
      if (which == WidgetHeader::VISIBLE_ONLY && isColumnHidden(i))
      {
        continue;
      }
      if (use_export_name)
      {
        header_labels << getHeaderExportName(i);
      }
      else 
      {
        header_labels << getHeaderName(i);
      }
    }
    return header_labels;
  }

  void TableView::setHeaderExportName(const int header_column, const QString& export_name)
  {
    QTableWidgetItem* ti = horizontalHeaderItem(header_column);
    if (ti == nullptr)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " not found!");
    }
    ti->setData(Qt::UserRole, export_name);
  }


  QString TableView::getHeaderExportName(const int header_column)
  {
    QTableWidgetItem* ti = horizontalHeaderItem(header_column);
    if (ti == nullptr)
    {
      throw  Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " not found!");
    }
    // prefer user role over display role
    if (ti->data(Qt::UserRole).isValid())
    {
      return ti->data(Qt::UserRole).toString();
    }
    else if (ti->data(Qt::DisplayRole).isValid())
    {
      return ti->data(Qt::DisplayRole).toString();
    }

    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " has no data!");
  }

  QString TableView::getHeaderName(const int header_column)
  {
    QTableWidgetItem* ti = horizontalHeaderItem(header_column);
    if (ti == nullptr)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " not found!");
    }
    if (ti->data(Qt::DisplayRole).isValid())
    {
      return ti->data(Qt::DisplayRole).toString();
    }

    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " has no data!");
  }

  void TableView::resizeEvent(QResizeEvent* /*event*/)
  {
    this->resizeColumnsToContents();

    int widgetWidth = this->viewport()->size().width();
    int tableWidth = 0;

    for (int i = 0; i < this->columnCount(); ++i)
    {
      tableWidth += this->horizontalHeader()->sectionSize(i);
    } //sections already resized to fit all data

    double scale = (double) widgetWidth / tableWidth;
    if (scale > 1.)
    {
      for (int i = 0; i < this->columnCount(); ++i)
      {
        this->setColumnWidth(i, this->horizontalHeader()->sectionSize(i) * scale);
      }
    }

    emit resized();
  }

}
