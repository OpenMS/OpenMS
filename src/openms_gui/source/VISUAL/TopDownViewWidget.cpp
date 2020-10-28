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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TopDownViewWidget.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMenu>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QFileDialog>
#include <QtCore/QTextStream>


#include <vector>

using namespace std;

///@improvement write the visibility-status of the columns in toppview.ini and read at start

//#define DEBUG_IDENTIFICATION_VIEW 1

namespace OpenMS
{
  TopDownViewWidget::TopDownViewWidget(const Param&, QWidget* parent) :
    QWidget(parent),
    ignore_update(false),
    layer_(nullptr)
  {
    // name can be displayed e.g., in a tab widget
    setObjectName("Top-Down Proteomics");

    QVBoxLayout* spectra_widget_layout = new QVBoxLayout(this);
    table_widget_ = new QTableWidget(this);
    table_widget_->setObjectName("table_widget");
    table_widget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to open it.");

    table_widget_->setSortingEnabled(true);

    table_widget_->setColumnWidth(0, 65); // MS Level
    table_widget_->setColumnWidth(1, 45); // index
    table_widget_->setColumnWidth(2, 70); // RT
    table_widget_->setColumnWidth(3, 70); // precursor m/z
    table_widget_->setColumnWidth(4, 55); // dissociation
    table_widget_->setColumnHidden(4, true);
    table_widget_->setColumnWidth(5, 45); // scan type
    table_widget_->setColumnHidden(5, true);
    table_widget_->setColumnWidth(6, 45);
    table_widget_->setColumnHidden(6, true);
    table_widget_->setColumnWidth(7, 45);
    table_widget_->setColumnWidth(8, 45);
    table_widget_->setColumnWidth(9, 45);
    table_widget_->setColumnWidth(10, 400);
    table_widget_->setColumnWidth(11, 45);
    table_widget_->setColumnWidth(12, 45);
    table_widget_->setColumnWidth(13, 45);

    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan type" << "zoom" << "score" << "rank" << "charge" << "sequence" << "accessions" << "#ID" << "#PH";
    table_widget_->setHorizontalHeaderLabels(header_labels);
    table_widget_->setColumnCount(header_labels.size());

    table_widget_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_widget_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_widget_->setShowGrid(false);

    spectra_widget_layout->addWidget(table_widget_);

    ////////////////////////////////////
    // additional checkboxes or buttons
    QHBoxLayout* tmp_hbox_layout = new QHBoxLayout();
    cb_deconvolute_ = new QCheckBox("Deconvolute", this);
    cb_deconvolute_->setChecked(false);

    tmp_hbox_layout->addWidget(cb_deconvolute_);

    spectra_widget_layout->addLayout(tmp_hbox_layout);
    table_widget_->sortByColumn(2, Qt::AscendingOrder);

    table_widget_->setEditTriggers(QAbstractItemView::NoEditTriggers);

    // select single rows
    table_widget_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_widget_->setSelectionMode(QAbstractItemView::SingleSelection);

    table_widget_->horizontalHeader()->setSectionsMovable(true);

    // header context menu
    table_widget_->horizontalHeader()->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(table_widget_->horizontalHeader(), &QHeaderView::customContextMenuRequested, this, &TopDownViewWidget::headerContextMenu_);
    connect(table_widget_, &QTableWidget::cellClicked, this, &TopDownViewWidget::cellClicked_);
    connect(table_widget_, &QTableWidget::currentItemChanged, this, &TopDownViewWidget::spectrumSelectionChange_);    
    connect(cb_deconvolute_, &QCheckBox::toggled, this, &TopDownViewWidget::updateEntries);
  }

  QTableWidget* TopDownViewWidget::getTableWidget()
  {
    return table_widget_;
  }

  void TopDownViewWidget::clear()
  {
    // remove all entries
    setLayer(nullptr);
  }

  void TopDownViewWidget::cellClicked_(int row, int column)
  {
    if (row >= table_widget_->rowCount()
    || column >= table_widget_->columnCount()
    || table_widget_->horizontalHeaderItem(column) == nullptr)
    {
      return;
    }

    int ms2_spectrum_index = table_widget_->item(row, 1)->data(Qt::DisplayRole).toInt();

    if (table_widget_->horizontalHeaderItem(column)->text() == "precursor m/z")
    {
      const auto& ms2_spectrum = (*layer_->getPeakData())[ms2_spectrum_index];
      if (ms2_spectrum.getPrecursors().empty()) { return; }  // no precursor
      
      // has precursor:
      // determine parent MS1 spectrum of current MS2 row  // TODO: use method to get precursor spectrum instead
      int ms1_spectrum_index = 0;
      for (ms1_spectrum_index = ms2_spectrum_index; ms1_spectrum_index >= 0; --ms1_spectrum_index)
      { 
        if ((*layer_->getPeakData())[ms1_spectrum_index].getMSLevel() == 1)
        {
          break;
        }
      }

      if (ms1_spectrum_index != -1)
      {
        double precursor_mz = ms2_spectrum.getPrecursors()[0].getMZ();
        // determine start and stop of isolation window
        double isolation_window_lower_mz = precursor_mz - ms2_spectrum.getPrecursors()[0].getIsolationWindowLowerOffset();
        double isolation_window_upper_mz = precursor_mz + ms2_spectrum.getPrecursors()[0].getIsolationWindowUpperOffset();

        emit spectrumDeselected(ms2_spectrum_index);
        emit spectrumSelected(ms1_spectrum_index);
      }
    }
  }

  void TopDownViewWidget::spectrumSelectionChange_(QTableWidgetItem* current, QTableWidgetItem* previous)
  {
    /*test for previous == 0 is important - without it,
      the wrong spectrum will be selected after finishing
      the execution of a TOPP tool on the whole data */
    if (current == nullptr || previous == nullptr) { return; }

    int previous_spectrum_index = table_widget_->item(previous->row(), 1)->data(Qt::DisplayRole).toInt();
    int current_spectrum_index = table_widget_->item(current->row(), 1)->data(Qt::DisplayRole).toInt();

#ifdef DEBUG_IDENTIFICATION_VIEW
    cout << "previous/now: " << previous_spectrum_index << "\t" << current_spectrum_index << endl;
#endif

    emit spectrumDeselected(previous_spectrum_index);

    if (current->column() == 3) // precursor mz column clicked
    {
      // handled by cell click event
    }
    else // !precursor mz column clicked
    {
      emit spectrumSelected(current_spectrum_index);
    }
  }

  void TopDownViewWidget::setLayer(LayerData* cl)
  {
    // do not try to be smart and check if layer_ == cl; to return early
    // since the layer content might have changed, e.g. pepIDs were added
    layer_ = cl;
    updateEntries();
  }

  LayerData* TopDownViewWidget::getLayer()
  {
    return layer_;
  }

  void TopDownViewWidget::updateEntries()
  {
    has_data_ = false;

    // no valid peak layer attached
    if (layer_ == nullptr
    || layer_->getPeakData()->size() == 0
    || layer_->type != LayerData::DT_PEAK)
    {
      table_widget_->clear();
      return;
    }

    if (ignore_update) { return; }

    if (!isVisible()) { return; }

    // create header labels (setting header labels must occur after fill)
    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "precursor charge"; 

    table_widget_->clear();
    table_widget_->setRowCount(0);

    table_widget_->verticalHeader()->setHidden(true); // hide vertical column
    table_widget_->setColumnCount(header_labels.size());
    table_widget_->setColumnWidth(0, 65);
    table_widget_->setColumnWidth(1, 45);
    table_widget_->setColumnWidth(2, 70);
    table_widget_->setColumnWidth(3, 70);
    table_widget_->setColumnWidth(4, 55);

    QTableWidgetItem* proto_item = new QTableWidgetItem();
    proto_item->setTextAlignment(Qt::AlignCenter);
    table_widget_->setItemPrototype(proto_item);

    table_widget_->setSortingEnabled(false);
    table_widget_->setUpdatesEnabled(false);
    table_widget_->blockSignals(true);

    // generate flat list
    int selected_row(-1);
    for (Size i = 0; i < layer_->getPeakData()->size(); ++i)
    {
      QTableWidgetItem* item = table_widget_->itemPrototype()->clone();
      QColor c{255, 255, 255};

      const MSSpectrum & spectrum = (*layer_->getPeakData())[i];
      const UInt ms_level = spectrum.getMSLevel();
      const vector<Precursor> & precursors = spectrum.getPrecursors();

      // add new row at the end of the table
      table_widget_->insertRow(table_widget_->rowCount());

      // ms level
      addTextItemToBottomRow_(QString::number(ms_level), 0, c);

      // index
      addIntItemToBottomRow_(static_cast<Int>(i), 1, c);

      // rt
      addDoubleItemToBottomRow_(spectrum.getRT(), 2, c);


      // fill precursor information in columns
      if (!precursors.empty())
      {
        const Precursor & first_precursor = precursors.front();

        // set precursor m/z
        item = table_widget_->itemPrototype()->clone();
        item->setData(Qt::DisplayRole, first_precursor.getMZ());
        item->setBackgroundColor(c);
        item->setTextColor(Qt::blue); // draw precursor information in blue
        table_widget_->setItem(table_widget_->rowCount() - 1, 3, item); // precursor m/z

        // set charge
        item = table_widget_->itemPrototype()->clone();
        item->setBackgroundColor(c);
        item->setData(Qt::DisplayRole, first_precursor.getCharge());
        table_widget_->setItem(table_widget_->rowCount() - 1, 4, item); // precursor charge
      }
      else
      { // has no precursor
        addTextItemToBottomRow_("-", 3, c); // precursor m/z
        addTextItemToBottomRow_("-", 4, c); // precursor charge
      }

      if (i == layer_->getCurrentSpectrumIndex())
      {
        selected_row = item->row(); // get model index of selected spectrum
      }
    }

    table_widget_->setSortingEnabled(true);
    table_widget_->setHorizontalHeaderLabels(header_labels);
    table_widget_->resizeColumnsToContents();

    if (selected_row != -1)  // select and scroll down to item
    {
      table_widget_->selectRow(selected_row);
      QTableWidgetItem* selected_item = table_widget_->item(selected_row, 0);
      selected_item->setSelected(true);
      table_widget_->setCurrentItem(selected_item);
      table_widget_->scrollToItem(selected_item);
    }

    table_widget_->blockSignals(false);
    table_widget_->setUpdatesEnabled(true);
    has_data_ = true;
  }

  void TopDownViewWidget::headerContextMenu_(const QPoint& pos)
  {
    // create menu
    QMenu context_menu(table_widget_);

    // extract header labels
    QStringList header_labels;
    for (int i = 0; i != table_widget_->columnCount(); ++i)
    {
      QTableWidgetItem* ti = table_widget_->horizontalHeaderItem(i);
      if (ti != nullptr)
      {
        header_labels.append(ti->text());
      }
    }

    // add actions which show/hide columns
    for (int i = 0; i != table_widget_->columnCount(); ++i)
    {
      QTableWidgetItem* ti = table_widget_->horizontalHeaderItem(i);
      if (ti == nullptr) continue;
      QAction* action = context_menu.addAction(ti->text(), [=]() {
              // invert visibility upon clicking the item
              table_widget_->setColumnHidden(i, !table_widget_->isColumnHidden(i));
        });
      action->setCheckable(true);
      action->setChecked(!table_widget_->isColumnHidden(i));
    }
    context_menu.exec(table_widget_->mapToGlobal(pos));
  }

  void TopDownViewWidget::addTextItemToBottomRow_(const QString& text, Size column_index, const QColor& c)
  {
    QTableWidgetItem * item = table_widget_->itemPrototype()->clone();
    item->setText(text);
    item->setBackgroundColor(c);
    table_widget_->setItem(table_widget_->rowCount() - 1, column_index, item);
  }

  void TopDownViewWidget::addIntItemToBottomRow_(const Int i, Size column_index, const QColor& c)
  {
    QTableWidgetItem * item = table_widget_->itemPrototype()->clone();
    item->setData(Qt::DisplayRole, i);
    item->setBackgroundColor(c);
    table_widget_->setItem(table_widget_->rowCount() - 1, column_index, item);
  }

  void TopDownViewWidget::addDoubleItemToBottomRow_(const double d, Size column_index, const QColor& c)
  {
    QTableWidgetItem * item = table_widget_->itemPrototype()->clone();
    item->setData(Qt::DisplayRole, d);
    item->setBackgroundColor(c);
    table_widget_->setItem(table_widget_->rowCount() - 1, column_index, item);
  }

  void TopDownViewWidget::addCheckboxItemToBottomRow_(bool selected,  Size column_index, const QColor& c)
  {
    QTableWidgetItem * item = table_widget_->itemPrototype()->clone();
    item->setCheckState(selected ? Qt::Checked : Qt::Unchecked);
    item->setBackgroundColor(c);
    table_widget_->setItem(table_widget_->rowCount() - 1, column_index, item);
  }

}
