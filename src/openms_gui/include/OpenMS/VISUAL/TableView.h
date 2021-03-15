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
    virtual ~TableView() = default;

    /**
      @brief Export table entries as currently shown in the table in TSV format (only for visible data)
      
      A filename will be queried using a dialog, before exporting.
    
      Headers will be exported using their export name (if available, see @p setHeaderExportName()).
      
      All cells will be queried for their Qt::UserRole, then for Qt::DisplayRole and last for Qt::CheckStateRole.
      The first item to return data will be used!
      Thus, to export data which differs from the visible (==DisplayRole), use QTableWidgetItem::setData(Qt::UserRole, ...).
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

  signals:
    /// emitted when the widget is resized
    void resized();

  protected:
    // emits the resized signal
    void resizeEvent(QResizeEvent* event) override;

  protected slots:
    /// Display header context menu; allows to show/hide columns
    void headerContextMenu_(const QPoint&);
  };
}
