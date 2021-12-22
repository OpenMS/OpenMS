// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/LayerStatisticsDialog.h>
#include <ui_LayerStatisticsDialog.h>

#include <OpenMS/VISUAL/PlotWidget.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>

#include <QtWidgets/QPushButton>

#include <array>
#include <variant>

using namespace std;

namespace OpenMS
{

  // helper for visitor pattern with std::visit
  template<class... Ts>
  struct overload : Ts... {
    using Ts::operator()...;
  };
  template<class... Ts>
  overload(Ts...) -> overload<Ts...>;

  /// stringify a number using thousand separator for better readability.
  /// The actual separator used depends on the system locale
  template<class T>
  QString toStringWithLocale(const T number)
  {
    std::stringstream iss;
    iss.imbue(std::locale("")); // use system locale, whatever it may be, e.g. "DE_de" 
    iss << number;
    return QString(iss.str().c_str());
    // custom locale is only valid for 'iss' and vanishes here
  }

  /// insert an intermediate row with a single spanning cell, which describes where the values came from (e.g. from DataArrays, or MetaData)
  void addHeaderRow(QTableWidget* table, int& row_i, const QString& row_name)
  {
    // clear row header name
    QTableWidgetItem* item = new QTableWidgetItem();
    item->setText("");
    table->setVerticalHeaderItem(row_i, item);
    
    item = new QTableWidgetItem();
    item->setText(row_name);
    //item->setBackgroundColor(QColor(Qt::darkGray));
    QFont font;
    font.setBold(true);
    item->setFont(font);
    item->setTextAlignment(Qt::AlignCenter);
    table->setItem(row_i, 0, item);
    table->setSpan(row_i, 0, 1, table->columnCount()); // extend row to span all columns
    ++row_i;
  }
  void addRow(LayerStatisticsDialog* lsd, QTableWidget* table, int& row_i, const QString& row_name, const StatsSummaryVariant& row_data, const bool enable_show_button = false)
  {
    // row name
    QTableWidgetItem* item = new QTableWidgetItem();
    item->setText(row_name);
    table->setVerticalHeaderItem(row_i, item);

    // column data
    constexpr int col_count = 4;
    // count, and min, max, avg (if present)
    std::array<QString, col_count> col_values = std::visit(overload{
      [&](const SSInt& d) -> std::array<QString, col_count> { return {toStringWithLocale(d.getCount()),
                                                                      QString::number(d.getMin()),
                                                                      QString::number(d.getMax()),
                                                                      QString::number(d.getAvg(), 'f', 2)}; },
      [&](const SSDouble& d) -> std::array<QString, col_count> { return {toStringWithLocale(d.getCount()),
                                                                         QString::number(d.getMin(), 'f', 2),
                                                                         QString::number(d.getMax(), 'f', 2),
                                                                         QString::number(d.getAvg(), 'f', 2)}; },
      [&](const StatsCounter& c) -> std::array<QString, col_count> { return {toStringWithLocale(c.counter),
                                                                             "-",
                                                                             "-",
                                                                             "-"}; }
      }, row_data);
    for (int col = 0; col < col_count; ++col)
    {
      item = new QTableWidgetItem();
      item->setText(col_values[col]);
      table->setItem(row_i, col, item);
    }
    bool show_button = std::visit(overload{
      [&](auto&& stats) { return stats.getCount() > 1 && stats.getMin() < stats.getMax(); },// for SSInt, SSDouble
      [&](const StatsCounter& /*c*/) { return false; }
      }, row_data);
    if (enable_show_button && show_button)
    {
      auto button = new QPushButton(row_name, table);
      table->setCellWidget(row_i, col_count, button);
      QObject::connect(button, SIGNAL(clicked()), lsd, SLOT(showDistribution_()));
    }
    
    // next row
    ++row_i;
  }

  LayerStatisticsDialog::LayerStatisticsDialog(PlotWidget* parent) :
    QDialog(parent),
    canvas_(parent->canvas()),
    ui_(new Ui::LayerStatisticsDialogTemplate)
  {
    ui_->setupUi(this);
    
    LayerStatistics ls;
    canvas_->getCurrentLayer().computeStats(ls);

    const auto& stats_meta = ls.getMetaStats(); // a bit expensive to get; so grab only once
    // compute total number of rows from all categories
    size_t total_rows = ls.getCoreStats().size() + stats_meta.size() + ls.getArrayStats().size();
    size_t header_rows = !stats_meta.empty() + !ls.getArrayStats().empty();
    ui_->table_->setRowCount((int) (total_rows + header_rows));
    // add each row
    int row_i = 0;
    for (const auto& item : ls.getCoreStats())
    {
      bool show_button = (item.first == "intensity");
      addRow(this, ui_->table_, row_i, QString(item.first.c_str()), item.second, show_button);
    }
    if (!stats_meta.empty())
    {
      addHeaderRow(ui_->table_, row_i, "Meta values");
      for (const auto& item : stats_meta)
      {
        addRow(this, ui_->table_, row_i, QString(item.first.c_str()), item.second, true);// with button to show distribution
      }
    }
    
    if (!ls.getArrayStats().empty())
    {
      addHeaderRow(ui_->table_, row_i, "DataArray values");
      for (const auto& item : ls.getArrayStats())
      {
        addRow(this, ui_->table_, row_i, QString(item.first.c_str()), item.second);
      }
    }
  }
  
  LayerStatisticsDialog::~LayerStatisticsDialog()
  {
    delete ui_;
  }

  void LayerStatisticsDialog::showDistribution_()
  {
    QPushButton* button = qobject_cast<QPushButton*>(sender());
    QString text = button->text();

    if (text == "intensity")
    {
      qobject_cast<PlotWidget*>(parent())->showIntensityDistribution();
    }
    else
    {
      qobject_cast<PlotWidget*>(parent())->showMetaDistribution(String(text));
    }
  }

} // namespace
