// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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


  void showDistribution(LayerStatisticsDialog* lsd, const QString& text,
                        const Math::Histogram<>& hist)
  {
    if (text == "intensity")
    {
      qobject_cast<PlotWidget*>(lsd->parent())->showIntensityDistribution(hist);
    }
    else
    {
      qobject_cast<PlotWidget*>(lsd->parent())->showMetaDistribution(String(text), hist);
    }
  }

  void addEmptyRow(QTableWidget* table, const int row_i, const QString& row_name)
  {
    table->setRowCount(row_i + 1);

    // set row header name
    QTableWidgetItem* item = new QTableWidgetItem();
    item->setText(row_name);
    table->setVerticalHeaderItem(row_i, item);
  }

  constexpr int col_count = 4; // columns: count, min, max, avg

  void populateRow(QTableWidget* table, const int row_i, const std::array<QString, col_count>& data)
  {
    for (int col = 0; col < col_count; ++col)
    {
      QTableWidgetItem* item = new QTableWidgetItem();
      item->setText(data[col]);
      table->setItem(row_i, col, item);
    }
  }

  /// insert an intermediate row with a single spanning cell, which describes where the values came from (e.g. from DataArrays, or MetaData)
  void addHeaderRow(QTableWidget* table, int& row_i, const QString& row_name)
  {
    addEmptyRow(table, row_i, "");
                                
    QTableWidgetItem* item = new QTableWidgetItem();
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

  void addRangeRow(LayerStatisticsDialog* lsd, QTableWidget* table, int& row_i,
                   const RangeStatsType& row_name, const RangeStatsVariant& row_data,
                   const bool enable_show_button, LayerStatistics* stats)
  {
    addEmptyRow(table, row_i, row_name.name.c_str());
    // get column data
    std::array<QString, col_count> col_values = std::visit(overload{
      [&](const RangeStatsInt& d) -> std::array<QString, col_count> { return {toStringWithLocale(d.getCount()),
                                                                      QString::number(d.getMin()),
                                                                      QString::number(d.getMax()),
                                                                      QString::number(d.getAvg(), 'f', 2)}; },
      [&](const RangeStatsDouble& d) -> std::array<QString, col_count> { return {toStringWithLocale(d.getCount()),
                                                                         QString::number(d.getMin(), 'f', 2),
                                                                         QString::number(d.getMax(), 'f', 2),
                                                                         QString::number(d.getAvg(), 'f', 2)}; }
      }, row_data);
    
    populateRow(table, row_i, col_values);

    if (enable_show_button)
    {
      auto button = new QPushButton(row_name.name.c_str(), table);
      table->setCellWidget(row_i, col_count, button);
      QObject::connect(button, &QPushButton::clicked, [=]() { showDistribution(lsd, row_name.name.c_str(), stats->getDistribution(row_name)); });
    }
    
    // next row
    ++row_i;
  }

  void addCountRow(QTableWidget* table, int& row_i, const QString& row_name, const StatsCounter& row_data)
  {
    addEmptyRow(table, row_i, row_name);

    // column data
    populateRow(table, row_i, {toStringWithLocale(row_data.counter), "-", "-", "-"});
    
    // next row
    ++row_i;
  } 

  LayerStatisticsDialog::LayerStatisticsDialog(PlotWidget* parent, std::unique_ptr<LayerStatistics>&& stats) :
    QDialog(parent),
    stats_(std::move(stats)),
    ui_(new Ui::LayerStatisticsDialogTemplate)
  {
    ui_->setupUi(this);

    ui_->table_->setColumnCount(col_count + 1); // +1 for button column

    const auto& stats_range = stats_->getRangeStatistics();
    const auto& stats_count = stats_->getCountStatistics();
    // add each row
    int row_i = 0;
    RangeStatsSource old_category = RangeStatsSource::SIZE_OF_STATSSOURCE;
    for (const auto& item : stats_range)
    {
      // add sections (relies on items being sorted!)
      if (old_category != item.first.src)
      {
        addHeaderRow(ui_->table_, row_i, StatsSourceNames[(size_t)item.first.src]);
        old_category = item.first.src;
      }
      bool show_button = (item.first == RangeStatsType{RangeStatsSource::CORE, "intensity"}) || item.first.src == RangeStatsSource::METAINFO;
      addRangeRow(this, ui_->table_, row_i, item.first, item.second, show_button, stats_.get());
    }

    if (!stats_count.empty())
    {
      addHeaderRow(ui_->table_, row_i, "Meta count values");
      for (const auto& item : stats_count)
      {
        addCountRow(ui_->table_, row_i, QString(item.first.c_str()), item.second);
      }
    }
  }
  
  LayerStatisticsDialog::~LayerStatisticsDialog()
  {
    delete ui_;
  }


} // namespace
