// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/VISUAL/SpectraIdentificationViewWidget.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>

#include <QtGui/QVBoxLayout>
#include <QtGui/QTreeWidget>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtGui/QHeaderView>
#include <QtGui/QMenu>
#include <QtGui/QPushButton>
#include <QtGui/QFileDialog>
#include <QtCore/QTextStream>


#include <vector>

using namespace std;

///@improvement write the visibility-status of the columns in toppview.ini and read at start

namespace OpenMS
{
  SpectraIdentificationViewWidget::SpectraIdentificationViewWidget(const Param&, QWidget* parent) :
    QWidget(parent),
    DefaultParamHandler("SpectraIdentificationViewWidget"),
    ignore_update(false),
    layer_(0),
    is_ms1_shown_(false)
  {
    // set common defaults
    defaults_.setValue("default_path", ".", "Default path for loading/storing data.");

    // id view
    defaults_.setValue("a_intensity", 1.0, "Default intensity of a-ions");
    defaults_.setValue("b_intensity", 1.0, "Default intensity of b-ions");
    defaults_.setValue("c_intensity", 1.0, "Default intensity of c-ions");
    defaults_.setValue("x_intensity", 1.0, "Default intensity of x-ions");
    defaults_.setValue("y_intensity", 1.0, "Default intensity of y-ions");
    defaults_.setValue("z_intensity", 1.0, "Default intensity of z-ions");
    defaults_.setValue("relative_loss_intensity", 0.1, "Relative loss in percent");
    defaults_.setValue("max_isotope", 2, "Maximum number of isotopes");
    defaults_.setValue("charge", 1, "Charge state");
    defaults_.setValue("show_a_ions", "false", "Show a-ions");
    defaults_.setValue("show_b_ions", "true", "Show b-ions");
    defaults_.setValue("show_c_ions", "false", "Show c-ions");
    defaults_.setValue("show_x_ions", "false", "Show x-ions");
    defaults_.setValue("show_y_ions", "true", "Show y-ions");
    defaults_.setValue("show_z_ions", "false", "Show z-ions");
    defaults_.setValue("show_precursor", "false", "Show precursor");
    defaults_.setValue("add_losses", "false", "Show neutral losses");
    defaults_.setValue("add_isotopes", "false", "Show isotopes");
    defaults_.setValue("add_abundant_immonium_ions", "false", "Show abundant immonium ions");
    defaults_.setValue("tolerance", 0.5, "Mass tolerance in Th used in the automatic alignment."); // unfortunately we don't support alignment with ppm error

    QVBoxLayout* spectra_widget_layout = new QVBoxLayout(this);
    table_widget_ = new QTableWidget(this);
    table_widget_->setObjectName("table_widget");
    table_widget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to open it.");

    table_widget_->setSortingEnabled(true);

    table_widget_->setColumnWidth(0, 65); //MS Level
    table_widget_->setColumnWidth(1, 45); //index
    table_widget_->setColumnWidth(2, 70);
    table_widget_->setColumnWidth(3, 70);
    table_widget_->setColumnWidth(4, 55);
    table_widget_->setColumnHidden(4, true);
    table_widget_->setColumnWidth(5, 45);
    table_widget_->setColumnHidden(5, true);
    table_widget_->setColumnWidth(6, 45);
    table_widget_->setColumnHidden(6, true);
    table_widget_->setColumnWidth(7, 45);
    table_widget_->setColumnWidth(8, 45);
    table_widget_->setColumnWidth(9, 45);
    table_widget_->setColumnWidth(10, 400);
    table_widget_->setColumnWidth(11, 45);

    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan type" << "zoom" << "score" << "rank" << "charge" << "sequence" << "accessions";
    table_widget_->setHorizontalHeaderLabels(header_labels);
    table_widget_->setColumnCount(header_labels.size());

    table_widget_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_widget_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_widget_->setShowGrid(false);

    connect(table_widget_, SIGNAL(currentItemChanged(QTableWidgetItem*, QTableWidgetItem*)), this, SLOT(spectrumSelectionChange_(QTableWidgetItem*, QTableWidgetItem*)));

    spectra_widget_layout->addWidget(table_widget_);

    ////////////////////////////////////
    // additional checkboxes and buttons
    QHBoxLayout* tmp_hbox_layout = new QHBoxLayout();

    hide_no_identification_ = new QCheckBox("Only hits", this);
    hide_no_identification_->setChecked(true);
    connect(hide_no_identification_, SIGNAL(toggled(bool)), this, SLOT(updateEntries()));

    create_rows_for_commmon_metavalue_ = new QCheckBox("Show advanced\nannotations", this);
    connect(create_rows_for_commmon_metavalue_, SIGNAL(toggled(bool)), this, SLOT(updateEntries()));

    QPushButton* save_idXML = new QPushButton("save idXML", this);
    connect(save_idXML, SIGNAL(clicked()), this, SLOT(saveIdXML_()));

    QPushButton* export_table = new QPushButton("export table", this);
    connect(export_table, SIGNAL(clicked()), this, SLOT(exportEntries_()));

    tmp_hbox_layout->addWidget(hide_no_identification_);
    tmp_hbox_layout->addWidget(create_rows_for_commmon_metavalue_);
    tmp_hbox_layout->addWidget(save_idXML);
    tmp_hbox_layout->addWidget(export_table);

    spectra_widget_layout->addLayout(tmp_hbox_layout);
    table_widget_->sortByColumn(2, Qt::AscendingOrder);

    table_widget_->setEditTriggers(QAbstractItemView::NoEditTriggers);

    // select single rows
    table_widget_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_widget_->setSelectionMode(QAbstractItemView::SingleSelection);

    table_widget_->horizontalHeader()->setMovable(true);

    // header context menu
    table_widget_->horizontalHeader()->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(table_widget_->horizontalHeader(), SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(headerContextMenu_(const QPoint &)));

    connect(table_widget_, SIGNAL(cellClicked(int, int)), this, SLOT(cellClicked_(int, int)));
  }

  QTableWidget* SpectraIdentificationViewWidget::getTableWidget()
  {
    return table_widget_;
  }

  void SpectraIdentificationViewWidget::cellClicked_(int row, int column)
  {
    if (row >= table_widget_->rowCount() || column >= table_widget_->columnCount())
    {
      return;
    }

    int ms2_spectrum_index = table_widget_->item(row, 1)->data(Qt::DisplayRole).toInt();

    if (column == 3) // precursor mz column
    {
      if (!(*layer_->getPeakData())[ms2_spectrum_index].getPrecursors().empty()) // has precursor
      {
        // determine parent MS1 spectrum of current MS2 row
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
          double precursor_mz = (*layer_->getPeakData())[ms2_spectrum_index].getPrecursors()[0].getMZ();
          // determine start and stop of isolation window
          double isolation_window_lower_mz = precursor_mz - (*layer_->getPeakData())[ms2_spectrum_index].getPrecursors()[0].getIsolationWindowLowerOffset();
          double isolation_window_upper_mz = precursor_mz + (*layer_->getPeakData())[ms2_spectrum_index].getPrecursors()[0].getIsolationWindowUpperOffset();

          if (!is_ms1_shown_)
          {
#ifdef DEBUG_IDENTIFICATION_VIEW
            cout << "cellClicked_ deselect MS2: " << ms2_spectrum_index << endl;
#endif
            emit spectrumDeselected(ms2_spectrum_index);
          }

#ifdef DEBUG_IDENTIFICATION_VIEW
          cout << "cellClicked_ select MS1: " << ms1_spectrum_index << endl;
#endif
          emit spectrumSelected(ms1_spectrum_index);
          is_ms1_shown_ = true;
          emit requestVisibleArea1D(isolation_window_lower_mz - 50.0, isolation_window_upper_mz +  50.0);
        }
      }
    }
  }

  void SpectraIdentificationViewWidget::spectrumSelectionChange_(QTableWidgetItem* current, QTableWidgetItem* previous)
  {
    /*test for previous == 0 is important - without it,
      the wrong spectrum will be selected after finishing
      the execution of a TOPP tool on the whole data */
    if (current == 0 || previous == 0)
    {
      return;
    }

    int previous_spectrum_index = table_widget_->item(previous->row(), 1)->data(Qt::DisplayRole).toInt();
    int current_spectrum_index = table_widget_->item(current->row(), 1)->data(Qt::DisplayRole).toInt();

    if (is_ms1_shown_)
    {
#ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "selection Change MS1 deselect: " << layer_->current_spectrum << endl;
#endif
      emit spectrumDeselected(int(layer_->getCurrentSpectrumIndex()));
    }
    else
    {
#ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "selection Change MS2 deselect: " << previous_spectrum_index << endl;
#endif
      emit spectrumDeselected(previous_spectrum_index);
    }

    if (current->column() == 3) // precursor mz column clicked
    {
      // handled by cell click event
    }
    else // !precursor mz column clicked
    {
#ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "selection Change MS2 select " << current_spectrum_index << endl;
#endif
      emit spectrumSelected(current_spectrum_index);
    }
  }

  void SpectraIdentificationViewWidget::attachLayer(LayerData* cl)
  {
    layer_ = cl;
  }

  void SpectraIdentificationViewWidget::updateEntries()
  {
    // no valid peak layer attached
    if (layer_ == 0 || layer_->getPeakData()->size() == 0 || layer_->type != LayerData::DT_PEAK)
    {
      table_widget_->clear();
      return;
    }

    if (ignore_update)
    {
      return;
    }

    if (!this->isVisible())
    {
      return;
    }

    set<String> common_keys;
    // determine meta values common to all hits
    if (create_rows_for_commmon_metavalue_->isChecked())
    {
      for (Size i = 0; i < layer_->getPeakData()->size(); ++i)
      {
        UInt ms_level = (*layer_->getPeakData())[i].getMSLevel();
        const vector<PeptideIdentification>& peptide_ids = (*layer_->getPeakData())[i].getPeptideIdentifications();

        if (ms_level != 2 || peptide_ids.size() == 0) // skip non ms2 spectra and spectra with no identification
        {
          continue;
        }

        for (vector<PeptideIdentification>::const_iterator pids_it = peptide_ids.begin(); pids_it != peptide_ids.end(); ++pids_it)
        {
          const vector<PeptideHit>& phits = pids_it->getHits();
          set<String> current_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<vector<PeptideHit>, set<String> >(phits.begin(), phits.end(), 100.0);
          if (common_keys.empty()) // first MS2 peptide hit found. Now insert keys.
          {
            swap(current_keys, common_keys);
          }
          else // calculate intersection between current keys and common keys -> set as common_keys
          {
            set<String> new_common_keys;
            set_intersection(current_keys.begin(), current_keys.end(), common_keys.begin(), common_keys.end(), inserter(new_common_keys, new_common_keys.begin()));
            swap(common_keys, new_common_keys);
          }
        }
      }
    }

    // create header labels (setting header labels must occur after fill)
    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan type" << "zoom" << "score" << "rank" << "charge" << "sequence" << "accession";
    for (set<String>::iterator sit = common_keys.begin(); sit != common_keys.end(); ++sit)
    {
      header_labels << sit->toQString();
    }

    table_widget_->clear();
    table_widget_->setRowCount(0);

    table_widget_->verticalHeader()->setHidden(true); // hide vertical column
    table_widget_->setColumnCount(header_labels.size());
    table_widget_->setColumnWidth(0, 65);
    table_widget_->setColumnWidth(1, 45);
    table_widget_->setColumnWidth(2, 70);
    table_widget_->setColumnWidth(3, 70);
    table_widget_->setColumnWidth(4, 55);
    table_widget_->setColumnHidden(4, true);
    table_widget_->setColumnWidth(5, 45);
    table_widget_->setColumnHidden(5, true);
    table_widget_->setColumnWidth(6, 45);
    table_widget_->setColumnHidden(6, true);
    table_widget_->setColumnWidth(7, 45);
    table_widget_->setColumnWidth(8, 45);
    table_widget_->setColumnHidden(8, true);
    table_widget_->setColumnWidth(9, 45);
    table_widget_->setColumnWidth(10, 400);
    table_widget_->setColumnWidth(11, 45);

    QTableWidgetItem* proto_item = new QTableWidgetItem();
    proto_item->setTextAlignment(Qt::AlignCenter);

    table_widget_->setItemPrototype(proto_item);

    table_widget_->setSortingEnabled(false);
    table_widget_->setUpdatesEnabled(false);
    table_widget_->blockSignals(true);

    if (layer_ == 0)
    {
      return;
    }

    QTableWidgetItem* item = 0;
    QTableWidgetItem* selected_item = 0;
    Size selected_row = 0;

    // generate flat list
    for (Size i = 0; i < layer_->getPeakData()->size(); ++i)
    {
      UInt ms_level = (*layer_->getPeakData())[i].getMSLevel();
      const vector<PeptideIdentification>& pi = (*layer_->getPeakData())[i].getPeptideIdentifications();
      Size id_count = pi.size();

      // allow only MS2    OR  MS1 with peptideIDs (from Mass Fingerprinting)
      if (ms_level != 2 && id_count == 0)
      {
        continue;
      }

      // skip
      if (hide_no_identification_->isChecked() && id_count == 0)
      {
        continue;
      }

      // coloring
      QColor c;

      // get peptide identifications of current spectrum
      if (id_count != 0)
      {
        c = Qt::green; // with identification
      }
      else
      {
        c = Qt::white; // without identification
      }

      // add new row at the end of the table
      table_widget_->insertRow(table_widget_->rowCount());

      // ms level
      item = table_widget_->itemPrototype()->clone();
      item->setText(QString::number(ms_level));
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount() - 1, 0, item);

      // index
      item = table_widget_->itemPrototype()->clone();
      item->setData(Qt::DisplayRole, Int(i));
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount() - 1, 1, item);

      // rt
      item = table_widget_->itemPrototype()->clone();
      item->setData(Qt::DisplayRole, (*layer_->getPeakData())[i].getRT());
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount() - 1, 2, item);

#ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "peptide identifications found:  " << id_count << endl;
#endif

      if (id_count != 0)
      {
#ifdef DEBUG_IDENTIFICATION_VIEW
        cout << "  peptide hits found: " << ph.size() << endl;
#endif

        PeptideHit best_ph;
        if (IDFilter().getBestHit(pi, false, best_ph))
        {
          // score
          item = table_widget_->itemPrototype()->clone();
          item->setData(Qt::DisplayRole, best_ph.getScore());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount() - 1, 7, item);

          // rank
          item = table_widget_->itemPrototype()->clone();
          item->setData(Qt::DisplayRole, best_ph.getRank());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount() - 1, 8, item);

          // charge
          item = table_widget_->itemPrototype()->clone();
          item->setData(Qt::DisplayRole, best_ph.getCharge());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount() - 1, 9, item);

          //sequence
          item = table_widget_->itemPrototype()->clone();
          item->setText(best_ph.getSequence().toString().toQString());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount() - 1, 10, item);

          //Accession
          item = table_widget_->itemPrototype()->clone();
          item->setTextAlignment(Qt::AlignLeft);

          set<String> protein_accessions = best_ph.extractProteinAccessions();
          String accessions = ListUtils::concatenate(vector<String>(protein_accessions.begin(), protein_accessions.end()), ", ");
          item->setText(accessions.toQString());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount() - 1, 11, item);

          // add additional meta value columns
          if (create_rows_for_commmon_metavalue_->isChecked())
          {
            Int current_col = 12;
            for (set<String>::iterator sit = common_keys.begin(); sit != common_keys.end(); ++sit)
            {
              DataValue dv = best_ph.getMetaValue(*sit);
              item = table_widget_->itemPrototype()->clone();
              item->setTextAlignment(Qt::AlignLeft);
              if (dv.valueType() == DataValue::DOUBLE_VALUE)
              {
                item->setData(Qt::DisplayRole, (double)dv);
              }
              else
              {
                item->setText(dv.toQString());
              }
              item->setBackgroundColor(c);
              table_widget_->setItem(table_widget_->rowCount() - 1, current_col, item);
              ++current_col;
            }
          }
        }
      }
      else // no identification
      {
        // score
        item = table_widget_->itemPrototype()->clone();
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount() - 1, 7, item);

        // rank
        item = table_widget_->itemPrototype()->clone();
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount() - 1, 8, item);

        // charge
        item = table_widget_->itemPrototype()->clone();
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount() - 1, 9, item);

        //sequence
        item = table_widget_->itemPrototype()->clone();
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount() - 1, 10, item);

        //accession
        item = table_widget_->itemPrototype()->clone();
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount() - 1, 11, item);

        // add additional meta value columns
        if (create_rows_for_commmon_metavalue_->isChecked())
        {
          Int current_col = 12;
          for (set<String>::iterator sit = common_keys.begin(); sit != common_keys.end(); ++sit)
          {
            item = table_widget_->itemPrototype()->clone();
            item->setTextAlignment(Qt::AlignLeft);
            item->setText("-");
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount() - 1, current_col, item);
            ++current_col;
          }
        }
      }

      if (!(*layer_->getPeakData())[i].getPrecursors().empty()) // has precursor
      {
        item = table_widget_->itemPrototype()->clone();
        item->setData(Qt::DisplayRole, (*layer_->getPeakData())[i].getPrecursors()[0].getMZ());
        item->setBackgroundColor(c);
        item->setTextColor(Qt::blue);
        table_widget_->setItem(table_widget_->rowCount() - 1, 3, item);

        item = table_widget_->itemPrototype()->clone();
        if (!(*layer_->getPeakData())[i].getPrecursors().front().getActivationMethods().empty())
        {
          QString t;
          for (std::set<Precursor::ActivationMethod>::const_iterator it = (*layer_->getPeakData())[i].getPrecursors().front().getActivationMethods().begin(); it != (*layer_->getPeakData())[i].getPrecursors().front().getActivationMethods().end(); ++it)
          {
            if (!t.isEmpty())
            {
              t.append(",");
            }
            t.append(QString::fromStdString((*layer_->getPeakData())[i].getPrecursors().front().NamesOfActivationMethod[*((*layer_->getPeakData())[i].getPrecursors().front().getActivationMethods().begin())]));
          }
          item->setText(t);
        }
        else
        {
          item->setText("-");
        }
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount() - 1, 4, item);
      }
      else // has no precursor (leave fields 3 and 4 empty)
      {
        item = table_widget_->itemPrototype()->clone();
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount() - 1, 3, item);

        item = table_widget_->itemPrototype()->clone();
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount() - 1, 4, item);
      }

      // scan mode
      item = table_widget_->itemPrototype()->clone();
      if ((*layer_->getPeakData())[i].getInstrumentSettings().getScanMode() > 0)
      {
        item->setText(QString::fromStdString((*layer_->getPeakData())[i].getInstrumentSettings().NamesOfScanMode[(*layer_->getPeakData())[i].getInstrumentSettings().getScanMode()]));
      }
      else
      {
        item->setText("-");
      }
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount() - 1, 5, item);

      // zoom scan
      item = table_widget_->itemPrototype()->clone();
      if ((*layer_->getPeakData())[i].getInstrumentSettings().getZoomScan())
      {
        item->setText("yes");
      }
      else
      {
        item->setText("no");
      }
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount() - 1, 6, item);

      if (i == layer_->getCurrentSpectrumIndex())
      {
        // just remember it, select later
        selected_item = item;
        selected_row = i;
      }
    }


    table_widget_->setSortingEnabled(true);
    table_widget_->setHorizontalHeaderLabels(header_labels);
    table_widget_->resizeColumnsToContents();

    if (selected_item)
    {
      // now, select and scroll down to item
      table_widget_->selectRow(int(selected_row));
      selected_item->setSelected(true);
      table_widget_->setCurrentItem(selected_item);
      table_widget_->scrollToItem(selected_item);
    }

    table_widget_->blockSignals(false);
    table_widget_->setUpdatesEnabled(true);
  }

  void SpectraIdentificationViewWidget::headerContextMenu_(const QPoint& pos)
  {
    // create menu
    QMenu* context_menu = new QMenu(table_widget_);

    // extract header labels
    QStringList header_labels;
    for (int i = 0; i != table_widget_->columnCount(); ++i)
    {
      QTableWidgetItem* ti = table_widget_->horizontalHeaderItem(i);
      if (ti != 0)
      {
        header_labels.append(ti->text());
      }
    }

    // add actions
    for (int i = 0; i < header_labels.size(); ++i)
    {
      QAction* tmp = new QAction(header_labels[i], context_menu);
      tmp->setCheckable(true);
      tmp->setChecked(!table_widget_->isColumnHidden(i));
      context_menu->addAction(tmp);
    }

    // show menu and hide selected columns
    QAction* selected = context_menu->exec(table_widget_->mapToGlobal(pos));
    if (selected != 0)
    {
      for (int i = 0; i < header_labels.size(); ++i)
      {
        if (selected->text() == header_labels[i])
        {
          selected->isChecked() ? table_widget_->setColumnHidden(i, false) : table_widget_->setColumnHidden(i, true);
        }
      }
    }
    delete (context_menu);
  }

  void SpectraIdentificationViewWidget::exportEntries_()
  {
    if (layer_ == 0 || layer_->getPeakData()->size() == 0 || layer_->type != LayerData::DT_PEAK)
    {
      return;
    }

    QString filename = QFileDialog::getSaveFileName(this, "Save File", "", "csv file (*.csv)");
    QFile f(filename);

    // extract header labels
    QStringList header_labels;
    for (int i = 0; i != table_widget_->columnCount(); ++i)
    {
      QTableWidgetItem* ti = table_widget_->horizontalHeaderItem(i);
      if (ti != 0)
      {
        header_labels.append(ti->text());
      }
    }

    if (f.open(QIODevice::WriteOnly))
    {
      QTextStream ts(&f);
      QStringList strList;

      // write header
      ts << header_labels.join("\t") + "\n";

      // write entries
      for (int r = 0; r < table_widget_->rowCount(); ++r)
      {
        strList.clear();
        for (int c = 0; c < table_widget_->columnCount(); ++c)
        {
          strList << table_widget_->item(r, c)->text();
        }
        ts << strList.join("\t") + "\n";
      }
      f.close();
    }
  }

  void SpectraIdentificationViewWidget::saveIdXML_()
  {
    // no valid peak layer attached
    if (layer_ == 0 || layer_->getPeakData()->size() == 0 || layer_->type != LayerData::DT_PEAK)
    {
      return;
    }

    QString filename = QFileDialog::getSaveFileName(this, "Save File", "", "idXML file (*.idXML)");
    vector<ProteinIdentification> prot_id = (*layer_->getPeakData()).getProteinIdentifications();
    vector<PeptideIdentification> all_pep_ids;
    for (int r = 0; r < table_widget_->rowCount(); ++r)
    {
      int spectrum_index = table_widget_->item(r, 1)->data(Qt::DisplayRole).toInt();
      vector<PeptideIdentification> pep_id = (*layer_->getPeakData())[spectrum_index].getPeptideIdentifications();
      copy(pep_id.begin(), pep_id.end(), back_inserter(all_pep_ids));
    }
    IdXMLFile().store(filename, prot_id, all_pep_ids);
  }

  SpectraIdentificationViewWidget::~SpectraIdentificationViewWidget()
  {
  }

}
