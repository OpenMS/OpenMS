// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <QTableWidget>

#include <OpenMS/VISUAL/MISC/CommonDefs.h>

namespace OpenMS
{
  /**
    @brief A better QTable for TOPPView, which supports exporting to TSV and conveniently adding data to cells and headers.
  */
  class TableView :
    public QTableWidget
  {
    Q_OBJECT
  public:
    /// Constructor
    TableView(QWidget* parent = nullptr);
    /// Destructor
    ~TableView() override = default;

    /**
      @brief Export table entries as currently shown in the table in TSV format (only for visible data)
      
      A filename will be queried using a dialog, before exporting.
    
      Headers will be exported using their export name (if available, see @p setHeaderExportName()).
      
      All cells will be queried for their Qt::UserRole, then for Qt::DisplayRole and last for Qt::CheckStateRole.
      The first item to return data will be used!
      Thus, to export data which differs from the visible (==DisplayRole), use QTableWidgetItem::setData(Qt::UserRole, ...).
      
      Note: to force export of hidden columns use @p setMandatoryExportColumns()
    */
    virtual void exportEntries();

    /// adds a new row to the bottom
    void appendRow();

    QTableWidgetItem* setAtBottomRow(const QString& text, size_t column_index, const QColor& background, const QColor& foreground = QColor("SomeInvalidColor"));
    QTableWidgetItem* setAtBottomRow(const char* text, size_t column_index, const QColor& background, const QColor& foreground = QColor("SomeInvalidColor"));
    QTableWidgetItem* setAtBottomRow(const int i, size_t column_index, const QColor& background, const QColor& foreground = QColor("SomeInvalidColor"));
    QTableWidgetItem* setAtBottomRow(const double d, size_t column_index, const QColor& background, const QColor& foreground = QColor("SomeInvalidColor"));
    /// create a checkbox item (with no text)
    QTableWidgetItem* setAtBottomRow(const bool selected, size_t column_index, const QColor& background, const QColor& foreground = QColor("SomeInvalidColor"));
    /// create a custom item (if above methods are not sufficient)
    QTableWidgetItem* setAtBottomRow(QTableWidgetItem* item, size_t column_index, const QColor& background, const QColor& foreground);

    /// if the item is purely a checkbox (e.g. added with setAtBottomRow(const bool selected, ...)),
    /// we set its DisplayRole to either '' or ' ', depending on checked state, to allow for row sorting 
    /// This function should be called whenever the check-state of the item changes
    static void updateCheckBoxItem(QTableWidgetItem* item);

    /// sets the visible headers (and the number of columns)
    void setHeaders(const QStringList& headers);
    
    /// hides columns with the given names
    /// @throws Exception::InvalidParameter if a name is not matching the current column names
    void hideColumns(const QStringList& header_names);

    /**
       @brief Obtain header names, either from all, or only the visible columns

       Headers can be obtained as shown (@p use_export_name = false) or for exporting to CSV
       where the alternative export name is preferred (if exists). See setHeaderExportName().

       @param which With or without invisible columns?
       @param use_export_name If column has a hidden export name, use that instead of the displayed name 
       @return List of header names 
    */
    QStringList getHeaderNames(const WidgetHeader which, bool use_export_name = false);

    /**
      @brief Set the export-name of a column, which will be returned in getHeaderNames() when @p use_export_name it true

      Export names are useful when exporting the table to CSV (see @p exportEntries()), and the column header should be a bit more verbose.

      Internally, this uses the Qt::UserRole's data to store the value.

      @param header_column Index of column
      @param export_name New export name to set

      @throws Exception::ElementNotFound if header at index @p header_column is not valid
    */
    void setHeaderExportName(const int header_column, const QString& export_name);

    /**
     @brief Gets the export-name of a column.

     Export names are useful when exporting the table to CSV (see @p exportEntries()), and the column header should be a bit more verbose.

     Internally, this queries the Qt::UserRole's data to get the value.
     If the export name was not set (using @p setHeaderExportName()), it returns the display name.

     @param header_column Index of column

     @throws Exception::ElementNotFound if header at index @p header_column is not valid
    */
    QString getHeaderExportName(const int header_column);

    /// get the displayed name of the header in column with index @p header_column
    /// @throws Exception::ElementNotFound if header at index @p header_column is not valid
    QString getHeaderName(const int header_column);

    /// Set the mandatory export columns @p cols which get exported even if the user decided to hide them.
    void setMandatoryExportColumns(QStringList& cols);    
  signals:
    /// emitted when the widget is resized
    void resized();

  protected:
    /// emits the resized signal
    void resizeEvent(QResizeEvent* event) override;

    /// columns that are exported to tsv files even if they are hidden in the GUI
    QStringList mandatory_export_columns_;
  protected slots:
    /// Display header context menu; allows to show/hide columns
    void headerContextMenu_(const QPoint&);
  };
}
