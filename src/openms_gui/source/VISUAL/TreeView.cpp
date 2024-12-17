// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TreeView.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Qt5Port.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QHeaderView>
#include <QMenu>

using namespace std;

///@improvement write the visibility-status of the columns in toppview.ini and read at start


namespace OpenMS
{
  TreeView::TreeView(QWidget* parent) :
    QTreeWidget(parent)
  {
    this->setObjectName("tree_widget");

    this->header()->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this->header(), &QHeaderView::customContextMenuRequested, this, &TreeView::headerContextMenu_);
  }


  void TreeView::headerContextMenu_(const QPoint& pos)
  {
    // allows to hide/show columns
    QMenu context_menu(this->header());
    const auto& header = this->headerItem();

    for (int i = 0; i < header->columnCount(); ++i)
    {
      auto action = context_menu.addAction(header->text(i), [i, this]() {
        this->setColumnHidden(i, !this->isColumnHidden(i));
        });
      action->setCheckable(true);
      action->setChecked(!this->isColumnHidden(i));
    }

    // show and execute menu
    context_menu.exec(this->mapToGlobal(pos));
  }

  void TreeView::setHeaders(const QStringList& headers)
  {
    setColumnCount(headers.size());
    setHeaderLabels(headers);
  }

  void TreeView::hideColumns(const QStringList& header_names)
  {
    auto hset = toQSet(header_names);
    // add actions which show/hide columns
    const auto& header = this->headerItem();

    for (int i = 0; i < header->columnCount(); ++i)
    {
      if (hset.contains(header->text(i)))
      {
        setColumnHidden(i, true);
        hset.remove(header->text(i));
      }
    }
    if (!hset.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "header_names contains a column name which is unknown: " + String(hset.values().join(", ")));
    }
  }

  QStringList TreeView::getHeaderNames(const WidgetHeader which) const
  {
    QStringList header_labels;
    for (int i = 0; i != columnCount(); ++i)
    {
      // do not export hidden columns
      if (which == WidgetHeader::VISIBLE_ONLY && isColumnHidden(i))
      {
        continue;
      }
      header_labels << getHeaderName(i);
    }
    return header_labels;
  }

  /// get the displayed name of the header in column with index @p header_column
  /// @throws Exception::ElementNotFound if header at index @p header_column is not valid

  QString TreeView::getHeaderName(const int header_column) const
  {
    const auto& header = this->headerItem();
    if (header->columnCount() <= header_column)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header index " + String(header_column) + " is too large. There are only " + String(header->columnCount()) + " columns!");
    }
    return header->text(header_column);
  }

}
