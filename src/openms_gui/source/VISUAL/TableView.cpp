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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TableView.h>

#include <OpenMS/CONCEPT/Exception.h>
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
      if (ti == nullptr) continue;
      QAction* action = context_menu.addAction(ti->text(), [=]() {
        // invert visibility upon clicking the item
        setColumnHidden(i, !isColumnHidden(i));
        });
      action->setCheckable(true);
      action->setChecked(!isColumnHidden(i));
    }
    context_menu.exec(mapToGlobal(pos));
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

    // write header
    ts << getHeaderNames(WidgetHeader::VISIBLE_ONLY, true).join("\t") + "\n";

    // write entries
    for (int r = 0; r < rowCount(); ++r)
    {
      for (int c = 0; c < columnCount(); ++c)
      {
        // do not export hidden columns
        if (isColumnHidden(c))
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
          else if (ti->data(Qt::DisplayRole).isValid())
          {
            str_list << ti->data(Qt::DisplayRole).toString();
          }
          else if (ti->data(Qt::CheckStateRole).isValid())
          {
            str_list << ti->data(Qt::CheckStateRole).toString();
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
    auto hset = header_names.toSet();
    // add actions which show/hide columns
    for (int i = 0; i != columnCount(); ++i)
    {
      QTableWidgetItem* ti = horizontalHeaderItem(i);
      if (ti == nullptr) continue;
      if (hset.contains(ti->text()))
      {
        setColumnHidden(i, true);
        hset.remove(ti->text());
      }
    }
    if (!hset.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "header_names contains a column name which is unknown: " + String(hset.toList().join(", ")));
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
    item->setBackgroundColor(background);
    if (foreground.isValid()) item->setForeground(QBrush(foreground));
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
      if (use_export_name) header_labels << getHeaderExportName(i);
      else header_labels << getHeaderName(i);
    }
    return header_labels;
  }

  void TableView::setHeaderExportName(const int header_column, const QString& export_name)
  {
    QTableWidgetItem* ti = horizontalHeaderItem(header_column);
    if (ti == nullptr) throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " not found!");

    ti->setData(Qt::UserRole, export_name);
  }


  QString TableView::getHeaderExportName(const int header_column)
  {
    QTableWidgetItem* ti = horizontalHeaderItem(header_column);
    if (ti == nullptr) throw  Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " not found!");

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
    if (ti == nullptr) throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " not found!");

    if (ti->data(Qt::DisplayRole).isValid())
    {
      return ti->data(Qt::DisplayRole).toString();
    }

    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header item " + String(header_column) + " has no data!");
  }

  void TableView::resizeEvent(QResizeEvent* /*event*/)
  {
    emit resized();
  }

}
