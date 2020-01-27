// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/FORMAT/MzIdentMLFile.h>
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
  SpectraIdentificationViewWidget::SpectraIdentificationViewWidget(const Param&, QWidget* parent) :
    QWidget(parent),
    DefaultParamHandler("SpectraIdentificationViewWidget"),
    ignore_update(false),
    layer_(nullptr),
    is_ms1_shown_(false),
    fragment_window_(nullptr)
  {
    setObjectName("Identifications");

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
    // additional checkboxes and buttons
    QHBoxLayout* tmp_hbox_layout = new QHBoxLayout();

    hide_no_identification_ = new QCheckBox("Only hits", this);
    hide_no_identification_->setChecked(true);


    create_rows_for_commmon_metavalue_ = new QCheckBox("Show advanced\nannotations", this);


    QPushButton* save_IDs = new QPushButton("Save IDs", this);
    connect(save_IDs, SIGNAL(clicked()), this, SLOT(saveIDs_()));

    QPushButton* export_table = new QPushButton("Export table", this);

    tmp_hbox_layout->addWidget(hide_no_identification_);
    tmp_hbox_layout->addWidget(create_rows_for_commmon_metavalue_);
    tmp_hbox_layout->addWidget(save_IDs);
    tmp_hbox_layout->addWidget(export_table);

    spectra_widget_layout->addLayout(tmp_hbox_layout);
    table_widget_->sortByColumn(2, Qt::AscendingOrder);

    table_widget_->setEditTriggers(QAbstractItemView::NoEditTriggers);

    // select single rows
    table_widget_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_widget_->setSelectionMode(QAbstractItemView::SingleSelection);

    table_widget_->horizontalHeader()->setSectionsMovable(true);

    // header context menu
    table_widget_->horizontalHeader()->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(table_widget_->horizontalHeader(), SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(headerContextMenu_(const QPoint &)));
    connect(table_widget_, SIGNAL(cellClicked(int, int)), this, SLOT(cellClicked_(int, int)));
    connect(table_widget_, SIGNAL(currentItemChanged(QTableWidgetItem*, QTableWidgetItem*)), this, SLOT(spectrumSelectionChange_(QTableWidgetItem*, QTableWidgetItem*)));
    connect(table_widget_, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(updateData_(QTableWidgetItem*)));
    connect(hide_no_identification_, SIGNAL(toggled(bool)), this, SLOT(updateEntries()));
    connect(create_rows_for_commmon_metavalue_, SIGNAL(toggled(bool)), this, SLOT(updateEntries()));
    connect(export_table, SIGNAL(clicked()), this, SLOT(exportEntries_()));
  }

  QTableWidget* SpectraIdentificationViewWidget::getTableWidget()
  {
    return table_widget_;
  }

  void SpectraIdentificationViewWidget::cellClicked_(int row, int column)
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
          emit spectrumSelected(ms1_spectrum_index, -1, -1); // no identification or hit selected (-1)
          is_ms1_shown_ = true;
          emit requestVisibleArea1D(isolation_window_lower_mz - 50.0, isolation_window_upper_mz +  50.0);
        }
      }
    }
    else if (table_widget_->horizontalHeaderItem(column)->text() == "PeakAnnotations")
    {
      QTableWidgetItem* current = table_widget_->item(row, column);

      int current_identification_index = table_widget_->item(current->row(), 12)->data(Qt::DisplayRole).toInt();  // peptide id. index

      if (current_identification_index < 0
       || current_identification_index >= static_cast<int>((*layer_->getPeakData())[ms2_spectrum_index].getPeptideIdentifications().size()))
       {
         return;
       }

      int current_peptide_hit_index = table_widget_->item(current->row(), 13)->data(Qt::DisplayRole).toInt();  // peptide hit index

      const vector<PeptideIdentification>& peptide_ids = (*layer_->getPeakData())[ms2_spectrum_index].getPeptideIdentifications();
      const vector<PeptideHit>& phits = peptide_ids[current_identification_index].getHits();

      if (current_peptide_hit_index < 0
       || current_peptide_hit_index >= static_cast<int>(phits.size()))
      {
        return;
      }
      const PeptideHit& hit = phits[current_peptide_hit_index];

      // initialize window, when the table is requested for the first time
      // afterwards the size will stay at the manually resized window size
      if (fragment_window_ == nullptr)
      {
        fragment_window_ = new QTableWidget();
        fragment_window_->resize(320, 500);

        fragment_window_->verticalHeader()->setHidden(true); // hide vertical column

        QStringList header_labels;
        header_labels << "m/z" << "name" << "intensity" << "charge";
        fragment_window_->setColumnCount(header_labels.size());
        fragment_window_->setHorizontalHeaderLabels(header_labels);

        QTableWidgetItem* proto_item = new QTableWidgetItem();
        proto_item->setTextAlignment(Qt::AlignCenter);
        fragment_window_->setItemPrototype(proto_item);
        fragment_window_->setSortingEnabled(true);
        fragment_window_->setWindowTitle(QApplication::translate("tr_fragment_annotation", "Peak Annotations"));
      }

      // reset table, if a new ID is chosen
      fragment_window_->setRowCount(0);

      for (const PeptideHit::PeakAnnotation & pa : hit.getPeakAnnotations())
      {
        fragment_window_->insertRow(fragment_window_->rowCount());
        QTableWidgetItem * item = fragment_window_->itemPrototype()->clone();
        item->setData(Qt::DisplayRole, pa.mz);
        fragment_window_->setItem(fragment_window_->rowCount() - 1, 0, item);
        item = fragment_window_->itemPrototype()->clone();
        item->setData(Qt::DisplayRole, pa.annotation.toQString());
        fragment_window_->setItem(fragment_window_->rowCount() - 1, 1, item);
        item = fragment_window_->itemPrototype()->clone();
        item->setData(Qt::DisplayRole, pa.intensity);
        fragment_window_->setItem(fragment_window_->rowCount() - 1, 2, item);
        item = fragment_window_->itemPrototype()->clone();
        item->setData(Qt::DisplayRole, pa.charge);
        fragment_window_->setItem(fragment_window_->rowCount() - 1, 3, item);
      }

      fragment_window_->resizeColumnsToContents();
      fragment_window_->resizeRowsToContents();
      fragment_window_->show();
      fragment_window_->setFocus(Qt::ActiveWindowFocusReason);
      QApplication::setActiveWindow(fragment_window_);
    }
  }

  void SpectraIdentificationViewWidget::spectrumSelectionChange_(QTableWidgetItem* current, QTableWidgetItem* previous)
  {
    /*test for previous == 0 is important - without it,
      the wrong spectrum will be selected after finishing
      the execution of a TOPP tool on the whole data */
    if (current == nullptr
    || previous == nullptr)
    {
      return;
    }

    int previous_spectrum_index = table_widget_->item(previous->row(), 1)->data(Qt::DisplayRole).toInt();
    int current_spectrum_index = table_widget_->item(current->row(), 1)->data(Qt::DisplayRole).toInt();
    int current_identification_index = table_widget_->item(current->row(), 12)->data(Qt::DisplayRole).toInt();  // peptide id. index
    int current_peptide_hit_index = table_widget_->item(current->row(), 13)->data(Qt::DisplayRole).toInt();  // peptide hit index

    if (is_ms1_shown_)
    {
#ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "selection Change MS1 deselect: " << layer_->getCurrentSpectrumIndex() << endl;
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
      emit spectrumSelected(current_spectrum_index, current_identification_index, current_peptide_hit_index);
    }
  }

  void SpectraIdentificationViewWidget::setLayer(LayerData* cl)
  {
    layer_ = cl;
    updateEntries();
  }

  LayerData* SpectraIdentificationViewWidget::getLayer()
  {
    return layer_;
  }

  void SpectraIdentificationViewWidget::updateEntries()
  {
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

    set<String> common_keys;
    bool has_peak_annotations(false);
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
          if (!has_peak_annotations && !phits[0].getPeakAnnotations().empty())
          {
            has_peak_annotations = true;
          }
        }
      }
    }

    // create header labels (setting header labels must occur after fill)
    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan type" << "zoom" << "score" << "rank" << "charge" << "sequence" << "accessions" << "#ID" << "#PH" << "Curated" << "precursor error (|ppm|)" << "precursor intensity";
    for (set<String>::iterator sit = common_keys.begin(); sit != common_keys.end(); ++sit)
    {
      header_labels << sit->toQString();
    }

    // add fragment annotation column, if fragment annotations were found in IDs
    if (has_peak_annotations)
    {
      header_labels << "PeakAnnotations";
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
    table_widget_->setColumnWidth(12, 45);
    table_widget_->setColumnWidth(13, 45);
    table_widget_->setColumnWidth(14, 45);
    table_widget_->setColumnWidth(15, 70);
    table_widget_->setColumnWidth(16, 70);

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
      QTableWidgetItem* item = nullptr;

      const MSSpectrum & spectrum = (*layer_->getPeakData())[i];
      const UInt ms_level = spectrum.getMSLevel();
      const vector<PeptideIdentification>& pi = spectrum.getPeptideIdentifications();
      const Size id_count = pi.size();
      const vector<Precursor> & precursors = spectrum.getPrecursors();


      // allow only MS2 OR MS1 with peptideIDs (from Mass Fingerprinting)
      if (ms_level != 2 && id_count == 0) { continue; }

      // skip
      if (hide_no_identification_->isChecked() && id_count == 0)  { continue; }

      // set row color
      QColor c;
      if (id_count == 0)
      {
        c = Qt::white; // no identification
      }
      else
      {
        c = Qt::green; // identification
      }

#ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "peptide identifications found:  " << id_count << endl;
#endif
      // get peptide identifications of current spectrum
      if (id_count == 0)
      {
        // add new row at the end of the table
        table_widget_->insertRow(table_widget_->rowCount());

        // ms level
        addTextItemToBottomRow_(QString::number(ms_level), 0, c);

        // index
        addIntItemToBottomRow_(static_cast<Int>(i), 1, c);

        // rt
        addDoubleItemToBottomRow_(spectrum.getRT(), 2, c);

        // score
        addTextItemToBottomRow_("-", 7, c);

        // rank
        addTextItemToBottomRow_("-", 8, c);

        // charge
        addTextItemToBottomRow_("-", 9, c);

        // sequence
        addTextItemToBottomRow_("-", 10, c);

        // accession
        addTextItemToBottomRow_("-", 11, c);

        // peptide identification index
        addTextItemToBottomRow_("-", 12, c);

        // peptide identification index
        addTextItemToBottomRow_("-", 13, c);

        // peptide identification index
        addTextItemToBottomRow_("-", 14, c);

        // ppm error
        addTextItemToBottomRow_("-", 15, c);

        // fill precursor information in columns
        if (!precursors.empty())
        {
          const Precursor & first_precursor = precursors.front();

          // set precursor m/z
          item = table_widget_->itemPrototype()->clone();
          item->setData(Qt::DisplayRole, first_precursor.getMZ());
          item->setBackgroundColor(c);
          item->setTextColor(Qt::blue); // draw precursor information in blue
          table_widget_->setItem(table_widget_->rowCount() - 1, 3, item);

          // set activation method
          item = table_widget_->itemPrototype()->clone();
          if (!first_precursor.getActivationMethods().empty())
          {
            QString t;
            for (auto it = first_precursor.getActivationMethods().begin();
              it != first_precursor.getActivationMethods().end();
              ++it)
            {
              if (!t.isEmpty()) { t.append(","); }
              t.append(QString::fromStdString(first_precursor.NamesOfActivationMethod[*first_precursor.getActivationMethods().begin()]));
            }
            item->setText(t);
          }
          else
          {
            item->setText("-");
          }
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount() - 1, 4, item);

          // set precursor intensity
          item = table_widget_->itemPrototype()->clone();
          item->setData(Qt::DisplayRole, first_precursor.getIntensity());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount() - 1, 17, item);
        }
        else
        { // has no precursor
          addTextItemToBottomRow_("-", 3, c);
          addTextItemToBottomRow_("-", 4, c);
          addTextItemToBottomRow_("-", 16, c); // precursor intensity
        }

        // scan mode
        QString scan_mode;
        if (spectrum.getInstrumentSettings().getScanMode() > 0)
        {
          scan_mode = QString::fromStdString(spectrum.getInstrumentSettings().NamesOfScanMode[spectrum.getInstrumentSettings().getScanMode()]);
        }
        else
        {
          scan_mode = "-";
        }
        addTextItemToBottomRow_(scan_mode, 5, c);

        // zoom scan
        QString is_zoom;
        if (spectrum.getInstrumentSettings().getZoomScan())
        {
          is_zoom = "yes";
        }
        else
        {
          is_zoom = "no";
        }
        addTextItemToBottomRow_(is_zoom, 6, c);
      }
      else
      {
        c = QColor(175, 255, 175); // with identification: light green

        for (Size pi_idx = 0; pi_idx != id_count; ++pi_idx)
        {
          for (Size ph_idx = 0; ph_idx != pi[pi_idx].getHits().size(); ++ph_idx)
          {
            const PeptideHit & ph = pi[pi_idx].getHits()[ph_idx];

            // add new row at the end of the table
            table_widget_->insertRow(table_widget_->rowCount());

            // ms level
            addTextItemToBottomRow_(QString::number(ms_level), 0, c);

            // index
            addIntItemToBottomRow_(static_cast<Int>(i), 1, c);

            // rt
            addDoubleItemToBottomRow_(spectrum.getRT(), 2, c);

            // score
            addDoubleItemToBottomRow_(ph.getScore(), 7, c);

            // rank
            addDoubleItemToBottomRow_(ph.getRank(), 8, c);

            // charge
            addDoubleItemToBottomRow_(ph.getCharge(), 9, c);

            //sequence
            String seq = ph.getSequence().toString();
            if (seq.empty()) seq = ph.getMetaValue("label");
            addTextItemToBottomRow_(seq.toQString(), 10, c);

            //Accession
            item = table_widget_->itemPrototype()->clone();
            item->setTextAlignment(Qt::AlignLeft);

            set<String> protein_accessions = ph.extractProteinAccessionsSet();
            String accessions = ListUtils::concatenate(vector<String>(protein_accessions.begin(), protein_accessions.end()), ", ");
            addTextItemToBottomRow_(accessions.toQString(), 11, c);

            // peptide identification index
            addIntItemToBottomRow_(static_cast<Int>(pi_idx), 12, c);

            // peptide identification index
            addIntItemToBottomRow_(static_cast<Int>(ph_idx), 13, c);

            bool selected(false);
            if (ph.metaValueExists("selected"))
            {
               selected = ph.getMetaValue("selected").toString() == "true";
            }
            addCheckboxItemToBottomRow_(selected, 14, c);

            // ppm error
            if (!precursors.empty()) // has precursor
            {
              const Precursor & first_precursor = precursors.front();
              // TODO compute theoretical precursor
              // following code is from precursor
              item = table_widget_->itemPrototype()->clone();
              item->setData(Qt::DisplayRole, first_precursor.getMZ());
              item->setBackgroundColor(c);
              item->setTextColor(Qt::blue);

              double ppm_error(0);

              // Protein:RNA cross-link, Protein-Protein cross-link, or other data with a precomputed precursor error
              if (pi[pi_idx].getHits()[0].metaValueExists(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM))
              {
                ppm_error = fabs((double)pi[pi_idx].getHits()[0].getMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM));
              }
              else if (!ph.getSequence().empty()) // works for normal linear fragments with the correct modifications included in the AASequence
              {
                double exp_precursor = first_precursor.getMZ();
                int charge = first_precursor.getCharge();
                double theo_mass = ph.getSequence().getMonoWeight();
                double theo_precursor= (theo_mass + (static_cast<double>(charge) * Constants::PROTON_MASS_U)) / static_cast<double>(charge);
                ppm_error = fabs((exp_precursor - theo_precursor) / exp_precursor / 1e-6);
              }
              addDoubleItemToBottomRow_(ppm_error, 15, c);
            }

            // fill precursor information in columns
            if (!precursors.empty()) // has precursor
            {
              const Precursor & first_precursor = precursors.front();
              item = table_widget_->itemPrototype()->clone();
              item->setData(Qt::DisplayRole, precursors.front().getMZ());
              item->setBackgroundColor(c);
              item->setTextColor(Qt::blue);
              table_widget_->setItem(table_widget_->rowCount() - 1, 3, item);

              item = table_widget_->itemPrototype()->clone();
              if (!first_precursor.getActivationMethods().empty())
              {
                QString t;
                for (auto it = first_precursor.getActivationMethods().begin();
                  it != first_precursor.getActivationMethods().end();
                  ++it)
                {
                  if (!t.isEmpty()) { t.append(","); }
                  t.append(QString::fromStdString(first_precursor.NamesOfActivationMethod[*first_precursor.getActivationMethods().begin()]));
                }
                item->setText(t);
              }
              else
              {
                item->setText("-");
              }
              item->setBackgroundColor(c);
              table_widget_->setItem(table_widget_->rowCount() - 1, 4, item);

              // set precursor intensity
              item = table_widget_->itemPrototype()->clone();
              item->setData(Qt::DisplayRole, first_precursor.getIntensity());
              item->setBackgroundColor(c);
              table_widget_->setItem(table_widget_->rowCount() - 1, 16, item);
            }
            else // has no precursor (leave fields 3 and 4 empty)
            {
              addTextItemToBottomRow_("-", 3, c);
              addTextItemToBottomRow_("-", 4, c);
              addTextItemToBottomRow_("-", 16, c); // precursor intensity
            }

            // add additional meta value columns
            if (create_rows_for_commmon_metavalue_->isChecked())
            {
              Int current_col = 17;
              for (set<String>::iterator sit = common_keys.begin(); sit != common_keys.end(); ++sit)
              {
                DataValue dv = ph.getMetaValue(*sit);
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
              // add peak annotation column
              if (has_peak_annotations)
              {
                addTextItemToBottomRow_("show", current_col, c);
              }
            }

            // scan mode
            QString scan_mode;
            if (spectrum.getInstrumentSettings().getScanMode() > 0)
            {
              scan_mode = QString::fromStdString(spectrum.getInstrumentSettings().NamesOfScanMode[spectrum.getInstrumentSettings().getScanMode()]);
            }
            else
            {
              scan_mode = "-";
            }
            addTextItemToBottomRow_(scan_mode, 5, c);

            // zoom scan
            QString is_zoom;
            if (spectrum.getInstrumentSettings().getZoomScan())
            {
              is_zoom = "yes";
            }
            else
            {
              is_zoom = "no";
            }
            addTextItemToBottomRow_(is_zoom, 6, c);
          }
        }
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
      if (ti != nullptr)
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
    if (selected != nullptr)
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
    if (layer_ == nullptr
      || layer_->getPeakData()->size() == 0
      || layer_->type != LayerData::DT_PEAK)
    {
      return;
    }

    QString filename = QFileDialog::getSaveFileName(this, "Save File", "", "csv file (*.csv)");
    QFile f(filename);

    // extract header labels
    QStringList header_labels;
    for (int i = 0; i != table_widget_->columnCount(); ++i)
    {
      // do not export hidden columns
      if (table_widget_->isColumnHidden(i))
      {
        continue;
      }

      QTableWidgetItem* ti = table_widget_->horizontalHeaderItem(i);
      if (ti != nullptr)
      {
        // add the format of this complex column to the header
        if (ti->text() == "PeakAnnotations")
        {
          header_labels.append("PeakAnnotations(mz|intensity|charge|annotation;)");
        }
        else
        {
          header_labels.append(ti->text());
        }
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
          // do not export hidden columns
          if (table_widget_->isColumnHidden(c))
          {
            continue;
          }

          QTableWidgetItem* ti = table_widget_->item(r, c);
          if (ti != nullptr)
          {
            // column 14 is the "Curated" column of checkboxes
            if (c == 14)
            {
              QString sel("0");
              if (ti->checkState() == Qt::Checked) // if the box is checked
              {
                sel = "1";
              }
              strList << sel;
            }
            else if (table_widget_->horizontalHeaderItem(c)->text() == "PeakAnnotations") // write out peak annotations instead of the "show" string in the table
            {
              int ms2_spectrum_index = table_widget_->item(r, 1)->data(Qt::DisplayRole).toInt();
              int current_identification_index = table_widget_->item(r, 12)->data(Qt::DisplayRole).toInt();  // peptide id. index

              if (current_identification_index < 0
               || current_identification_index >= static_cast<int>((*layer_->getPeakData())[ms2_spectrum_index].getPeptideIdentifications().size()))
               {
                 continue;
               }

              int current_peptide_hit_index = table_widget_->item(r, 13)->data(Qt::DisplayRole).toInt();  // peptide hit index

              const vector<PeptideIdentification>& peptide_ids = (*layer_->getPeakData())[ms2_spectrum_index].getPeptideIdentifications();
              const vector<PeptideHit>& phits = peptide_ids[current_identification_index].getHits();

              if (current_peptide_hit_index < 0
               || current_peptide_hit_index >= static_cast<int>(phits.size()))
              {
                continue;
              }
              const PeptideHit& hit = phits[current_peptide_hit_index];
              QString annotation = "";

              for (const PeptideHit::PeakAnnotation & pa : hit.getPeakAnnotations())
              {
                annotation += String(pa.mz).toQString() + "|" +
                                String(pa.intensity).toQString() + "|" +
                                String(pa.charge).toQString() + "|" +
                                pa.annotation.toQString() + ";";
              }
              strList << annotation;
            }
            else
            {
              strList << table_widget_->item(r, c)->text();
            }
          }
        }
        ts << strList.join("\t") + "\n";
      }
      f.close();
    }
  }
  void SpectraIdentificationViewWidget::addTextItemToBottomRow_(const QString& text, Size column_index, const QColor& c)
  {
    QTableWidgetItem * item = table_widget_->itemPrototype()->clone();
    item->setText(text);
    item->setBackgroundColor(c);
    table_widget_->setItem(table_widget_->rowCount() - 1, column_index, item);
  }

  void SpectraIdentificationViewWidget::addIntItemToBottomRow_(const Int i, Size column_index, const QColor& c)
  {
    QTableWidgetItem * item = table_widget_->itemPrototype()->clone();
    item->setData(Qt::DisplayRole, i);
    item->setBackgroundColor(c);
    table_widget_->setItem(table_widget_->rowCount() - 1, column_index, item);
  }

  void SpectraIdentificationViewWidget::addDoubleItemToBottomRow_(const double d, Size column_index, const QColor& c)
  {
    QTableWidgetItem * item = table_widget_->itemPrototype()->clone();
    item->setData(Qt::DisplayRole, d);
    item->setBackgroundColor(c);
    table_widget_->setItem(table_widget_->rowCount() - 1, column_index, item);
  }

  void SpectraIdentificationViewWidget::addCheckboxItemToBottomRow_(bool selected,  Size column_index, const QColor& c)
  {
    QTableWidgetItem * item = table_widget_->itemPrototype()->clone();
    if (selected)
    {
      item->setCheckState(Qt::Checked);
    }
    else
    {
      item->setCheckState(Qt::Unchecked);
    }
    item->setBackgroundColor(c);
    table_widget_->setItem(table_widget_->rowCount() - 1, column_index, item);
  }

  void SpectraIdentificationViewWidget::saveIDs_()
  {
    // no valid peak layer attached
    if (layer_ == nullptr || layer_->getPeakData()->size() == 0 || layer_->type != LayerData::DT_PEAK)
    {
      return;
    }

    // synchronize PeptideHits with the annotations in the spectrum
    layer_->synchronizePeakAnnotations();

    QString selectedFilter;
    QString filename = QFileDialog::getSaveFileName(this, "Save File", "", "idXML file (*.idXML);;mzIdentML file (*.mzid)", &selectedFilter);
    vector<ProteinIdentification> prot_id = (*layer_->getPeakData()).getProteinIdentifications();
    vector<PeptideIdentification> all_pep_ids;

    // collect PeptideIdentifications from each spectrum, while making sure each spectrum is only considered once
    // otherwise duplicates will be stored, if more than one PeptideHit is contained in a PeptideIdentification
    vector<int> added_spectra;
    for (int r = 0; r < table_widget_->rowCount(); ++r)
    {
      // get spectrum index of current table line
      int spectrum_index = table_widget_->item(r, 1)->data(Qt::DisplayRole).toInt();

      // skip this row, if this spectrum was already processed
      if (std::find(added_spectra.begin(), added_spectra.end(), spectrum_index) != added_spectra.end())
      {
        continue;
      }
      added_spectra.push_back(spectrum_index);

      // collect all PeptideIdentifications from this spectrum
      vector<PeptideIdentification> pep_id = (*layer_->getPeakData())[spectrum_index].getPeptideIdentifications();
      copy(pep_id.begin(), pep_id.end(), back_inserter(all_pep_ids));
    }

    if (String(filename).hasSuffix(String(".mzid")))
    {
      MzIdentMLFile().store(filename, prot_id, all_pep_ids);
    }
    else if (String(filename).hasSuffix(String(".idXML")))
    {
      IdXMLFile().store(filename, prot_id, all_pep_ids);
    }
    else if (String(selectedFilter).hasSubstring(String(".mzid")))
    {
      filename = filename + ".mzid";
      MzIdentMLFile().store(filename, prot_id, all_pep_ids);
    }
    else
    {
      filename = filename + ".idXML";
      IdXMLFile().store(filename, prot_id, all_pep_ids);
    }
  }

  // Upon changes in the table data (only possible by checking or unchecking a checkbox right now),
  // update the corresponding PeptideIdentification / PeptideHits
  void SpectraIdentificationViewWidget::updateData_(QTableWidgetItem* item)
  {
    // no valid peak layer attached
    if (layer_ == nullptr || layer_->getPeakData()->size() == 0 || layer_->type != LayerData::DT_PEAK)
    {
      return;
    }

    // Find column indices by names
    Size n_col = table_widget_->columnCount();
    Size id_col = 0;
    Size ph_col = 0;

    for (Size c = 0; c < n_col; ++c)
    {
      String col_head = table_widget_->horizontalHeaderItem(c)->text();
      if (col_head == "#ID")
      {
        id_col = c;
      }
      if (col_head == "#PH")
      {
        ph_col = c;
      }
    }

    // extract position of the correct Spectrum, PeptideIdentification and PeptideHit from the table
    int r = item->row();
    bool selected = item->checkState() == Qt::Checked;
    int spectrum_index = table_widget_->item(r, 1)->data(Qt::DisplayRole).toInt();
    int num_id = table_widget_->item(r, id_col)->text().toInt();
    int num_ph = table_widget_->item(r, ph_col)->text().toInt();

    vector<PeptideIdentification> pep_id = (*layer_->getPeakData())[spectrum_index].getPeptideIdentifications();

    // update "selected" value in the correct PeptideHits
    vector<PeptideHit> hits = pep_id[num_id].getHits();
    String sel = selected ? "true" : "false";

    if (hits[0].metaValueExists("xl_chain")) // XL-MS specific case, both PeptideHits belong to the same cross-link
    {
      hits[0].setMetaValue("selected", sel);
      if (hits.size() >= 2)
      {
        hits[1].setMetaValue("selected", sel);
      }
    }
    else // general case, update only the selected PepideHit
    {
      hits[num_ph].setMetaValue("selected", sel);
    }
    pep_id[num_id].setHits(hits);
    (*layer_->getPeakDataMuteable())[spectrum_index].setPeptideIdentifications(pep_id);

  }

  SpectraIdentificationViewWidget::~SpectraIdentificationViewWidget()
  {
  }

}
