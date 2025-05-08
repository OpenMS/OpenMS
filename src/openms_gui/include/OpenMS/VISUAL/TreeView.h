// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <QTreeWidget>

#include <OpenMS/VISUAL/MISC/CommonDefs.h>

namespace OpenMS
{
  /**
    @brief A better QTreeWidget for TOPPView, which supports header context menu and conveniently adding/getting headers names.
  */
  class TreeView :
    public QTreeWidget
  {
    Q_OBJECT
  public:
    /// Constructor
    TreeView(QWidget* parent = nullptr);
    /// Destructor
    ~TreeView() override = default;

    /// sets the visible headers (and the number of columns)
    void setHeaders(const QStringList& headers);

    /// hides columns with the given names
    /// @throws Exception::InvalidParameter if a name is not matching the current column names
    void hideColumns(const QStringList& header_names);

    /**
       @brief Obtain header names, either from all, or only the visible columns

       @param which With or without invisible columns?
       @return List of header names
    */
    QStringList getHeaderNames(const WidgetHeader which) const;

    /// get the displayed name of the header in column with index @p header_column
    /// @throws Exception::ElementNotFound if header at index @p header_column is not valid
    QString getHeaderName(const int header_column) const;

  private slots:
    /// Display header context menu; allows to show/hide columns
    void headerContextMenu_(const QPoint& pos);
  };
}
