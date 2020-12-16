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

#include <OpenMS/VISUAL/SpectraIDViewTab.h>

#include <OpenMS/VISUAL/TableView.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QFileDialog>

#include <vector>

using namespace std;

///@improvement write the visibility-status of the columns in toppview.ini and read at start

// Use a namespace to encapsulate names, yet use c-style 'enum' for fast conversion to int.
// So we can write: 'Clmn::MS_LEVEL', but get implicit conversion to int
namespace Clmn
{
  enum HeaderNames
  { // indices into QTableWidget's columns (which start at index 0)
    MS_LEVEL, SPEC_INDEX, RT, PRECURSOR_MZ, DISSOCIATION, SCANTYPE, ZOOM, SCORE, RANK, 
    CHARGE, SEQUENCE, ACCESSIONS, ID_NR, PEPHIT_NR, CURATED, PREC_PPM, PREC_INT, PEAK_ANNOTATIONS, /* last entry --> */ SIZE_OF_HEADERNAMES
  };
  // keep in SYNC with enum HeaderNames
  const QStringList HEADER_NAMES = QStringList()
                                    << "MS" << "index" << "RT"
                                    << "precursor m/z" << "dissociation" << "scan type" << "zoom" << "score"
                                    << "rank" << "charge" << "sequence" << "accessions" << "#ID" << "#PH"
                                    << "Curated" << "precursor error (|ppm|)" << "precursor intensity" << "peak annotations";
}

namespace OpenMS
{
  SpectraIDViewTab::SpectraIDViewTab(const Param&, QWidget* parent) :
    QWidget(parent),
    DefaultParamHandler("SpectraIDViewTab")
  {
    setObjectName("Identifications");

    // make sure they are in sync
    assert(Clmn::HEADER_NAMES.size() == Clmn::HeaderNames::SIZE_OF_HEADERNAMES);

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
    table_widget_ = new TableView(this);
    table_widget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to open it.");

    spectra_widget_layout->addWidget(table_widget_);

    ////////////////////////////////////
    // additional checkboxes and buttons
    QHBoxLayout* tmp_hbox_layout = new QHBoxLayout();

    hide_no_identification_ = new QCheckBox("Only hits", this);
    hide_no_identification_->setChecked(true);

    create_rows_for_commmon_metavalue_ = new QCheckBox("Show advanced\nannotations", this);

    QPushButton* save_IDs = new QPushButton("Save IDs", this);
    connect(save_IDs, &QPushButton::clicked, this, &SpectraIDViewTab::saveIDs_);

    QPushButton* export_table = new QPushButton("Export table", this);

    tmp_hbox_layout->addWidget(hide_no_identification_);
    tmp_hbox_layout->addWidget(create_rows_for_commmon_metavalue_);
    tmp_hbox_layout->addWidget(save_IDs);
    tmp_hbox_layout->addWidget(export_table);
    spectra_widget_layout->addLayout(tmp_hbox_layout);

    connect(table_widget_, &QTableWidget::currentCellChanged, this, &SpectraIDViewTab::currentCellChanged_);
    connect(table_widget_, &QTableWidget::itemChanged, this, &SpectraIDViewTab::updatedSingleCell_);
    connect(hide_no_identification_, &QCheckBox::toggled, this, &SpectraIDViewTab::updateEntries_);
    connect(create_rows_for_commmon_metavalue_, &QCheckBox::toggled, this, &SpectraIDViewTab::updateEntries_);
    connect(export_table, &QPushButton::clicked, table_widget_, &TableView::exportEntries);
  }

  void SpectraIDViewTab::clear()
  {
    table_widget_->clear();
    layer_ = nullptr;
  }

  void SpectraIDViewTab::currentCellChanged_(int row, int column, int /*old_row*/, int /*old_column*/)
  {
    // sometimes Qt calls this function when table empty during refreshing
    if (row < 0 || column < 0) return;

    if (row >= table_widget_->rowCount()
        ||  column >= table_widget_->columnCount())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid cell clicked.", String(row) + " " + column);
    }

    // deselect whatever is currently shown
    int last_spectrum_index = int(layer_->getCurrentSpectrumIndex());
    emit spectrumDeselected(last_spectrum_index);


    int current_spectrum_index = table_widget_->item(row, Clmn::SPEC_INDEX)->data(Qt::DisplayRole).toInt();
    const auto& exp = *layer_->getPeakData();
    const auto& spec2 = exp[current_spectrum_index];

    //
    // Signal for a new spectrum to be shown
    //
    // show precursor spectrum (usually MS1)
    if (column == Clmn::PRECURSOR_MZ)
    {
      const auto prec_it = exp.getPrecursorSpectrum(exp.begin() + current_spectrum_index);
      
      if (prec_it != exp.end() && !spec2.getPrecursors().empty())
      {
        double precursor_mz = spec2.getPrecursors()[0].getMZ();
        // determine start and stop of isolation window
        double isolation_window_lower_mz = precursor_mz - spec2.getPrecursors()[0].getIsolationWindowLowerOffset();
        double isolation_window_upper_mz = precursor_mz + spec2.getPrecursors()[0].getIsolationWindowUpperOffset();

        emit spectrumSelected(std::distance(exp.begin(), prec_it), -1, -1); // no identification or hit selected (-1)
        // zoom into precursor area
        emit requestVisibleArea1D(isolation_window_lower_mz - 50.0, isolation_window_upper_mz +  50.0);
      }
    }
    else
    { // if spectrum with no PepIDs is selected, there is nothing to show...
      auto item_pepid = table_widget_->item(row, Clmn::ID_NR);
      if (item_pepid == nullptr   // null for MS1 spectra
          || (!(item_pepid->data(Qt::DisplayRole).isValid())))
      {
        return;
      }
      int current_identification_index = item_pepid->data(Qt::DisplayRole).toInt();
      int current_peptide_hit_index = table_widget_->item(row, Clmn::PEPHIT_NR)->data(Qt::DisplayRole).toInt();
      emit spectrumSelected(current_spectrum_index, current_identification_index, current_peptide_hit_index);
    }

    //
    // show extra peak-fragment window
    //
    if (column == Clmn::PEAK_ANNOTATIONS 
        // column might not be present. Check the header name to make sure
        && table_widget_->horizontalHeaderItem(Clmn::PEAK_ANNOTATIONS)->text() == Clmn::HEADER_NAMES[Clmn::PEAK_ANNOTATIONS])
    {         

      auto item_pepid = table_widget_->item(row, Clmn::ID_NR);
      if (item_pepid)  // might be null for MS1 spectra
      {
        int current_identification_index = item_pepid->data(Qt::DisplayRole).toInt();
        int current_peptide_hit_index = table_widget_->item(row, Clmn::PEPHIT_NR)->data(Qt::DisplayRole).toInt();

        const vector<PeptideIdentification>& peptide_ids = spec2.getPeptideIdentifications();
        const vector<PeptideHit>& phits = peptide_ids[current_identification_index].getHits();
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
    } // PeakAnnotation cell clicked
  }

  bool SpectraIDViewTab::hasData(const LayerData* layer)
  {
    // this is a very easy check.
    // We do not check for PeptideIdentifications attached to Spectra, because the user could just 
    // want the list of unidentified MS2 spectra (obtained by unchecking the 'just hits' button).
    bool no_data = (layer == nullptr
                || (layer->type == LayerData::DT_PEAK && layer->getPeakData()->empty())
                || (layer->type == LayerData::DT_CHROMATOGRAM && layer->getChromatogramData()->empty()));
    return !no_data;
  }

  void SpectraIDViewTab::updateEntries(LayerData* cl)
  {
    // do not try to be smart and check if layer_ == cl; to return early
    // since the layer content might have changed, e.g. pepIDs were added
    layer_ = cl;
    updateEntries_(); // we need this extra function since its an internal slot
  }

  LayerData* SpectraIDViewTab::getLayer()
  {
    return layer_;
  }

  namespace Detail
  {
    template<>
    struct MetaKeyGetter<std::reference_wrapper<const PeptideHit>>
    {
      static void getKeys(const std::reference_wrapper<const PeptideHit>& object, std::vector<String>& keys)
      {
        object.get().getKeys(keys);
      };
    };
  }

  void SpectraIDViewTab::updateEntries_()
  {
    // no valid peak layer attached
    if (!hasData(layer_))
    {
      clear();
      return;
    }

    if (ignore_update) { return; }

    if (!isVisible()) { return; }

    int restore_spec_index = layer_->getCurrentSpectrumIndex();

    set<String> common_keys;
    bool has_peak_annotations(false);
    // determine meta values common to all hits
    Detail::MetaKeyGetter<std::reference_wrapper<const PeptideHit>> getter;
    if (create_rows_for_commmon_metavalue_->isChecked())
    {
      std::vector<std::reference_wrapper<const PeptideHit>> all_hits;

      for (const auto& spec : layer_->getPeakData()->getSpectra())
      {
        UInt ms_level = spec.getMSLevel();
        const vector<PeptideIdentification>& peptide_ids = spec.getPeptideIdentifications();

        if (ms_level != 2 || peptide_ids.size() == 0) // skip non ms2 spectra and spectra with no identification
        {
          continue;
        }

        for (const auto& pep_id : peptide_ids)
        {
          const vector<PeptideHit>& phits = pep_id.getHits();
          all_hits.insert(all_hits.end(), phits.begin(), phits.end());
          if (!has_peak_annotations && !phits[0].getPeakAnnotations().empty())
          {
            has_peak_annotations = true;
          }
        }
      }

      common_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<
                      std::vector<std::reference_wrapper<const PeptideHit>>,
                      set<String> >(all_hits.begin(), all_hits.end(), 100.0, getter);
    }

    // create header labels (setting header labels must occur after fill)
    QStringList headers = Clmn::HEADER_NAMES;
    if (!has_peak_annotations)
    { // remove peak annotations column                   
      headers.pop_back();
    }
    // add common meta columns (not indexed anymore, but we don't need them to be)
    for (const auto& ck : common_keys)
    {
      headers << ck.toQString();
    }

    table_widget_->clear();
    table_widget_->setRowCount(0);
    table_widget_->setColumnCount(headers.size());
    table_widget_->setSortingEnabled(false);
    table_widget_->setUpdatesEnabled(false);
    table_widget_->blockSignals(true);

    // generate flat list
    int selected_row(-1);
    // index i is needed, so iterate the old way...
    for (Size i = 0; i < layer_->getPeakData()->size(); ++i)
    {
      const MSSpectrum& spectrum = (*layer_->getPeakData())[i];
      const UInt ms_level = spectrum.getMSLevel();
      const vector<PeptideIdentification>& pi = spectrum.getPeptideIdentifications();
      const Size id_count = pi.size();
      const vector<Precursor> & precursors = spectrum.getPrecursors();

      // allow only MS2 OR MS1 with peptideIDs (from Mass Fingerprinting)
      if (ms_level != 2 && id_count == 0) { continue; }

      // skip
      if (hide_no_identification_->isChecked() && id_count == 0)  { continue; }

      // set row background color
      QColor bg_color = (id_count == 0 ? Qt::white : Qt::green);

      // get peptide identifications of current spectrum
      if (id_count == 0)
      {
        // add new row at the end of the table
        table_widget_->insertRow(table_widget_->rowCount());

        fillRow_(spectrum, i, bg_color);
      }
      else
      {
        for (Size pi_idx = 0; pi_idx != id_count; ++pi_idx)
        {
          for (Size ph_idx = 0; ph_idx != pi[pi_idx].getHits().size(); ++ph_idx)
          {
            const PeptideHit& ph = pi[pi_idx].getHits()[ph_idx];

            // add new row at the end of the table
            table_widget_->insertRow(table_widget_->rowCount());

            fillRow_(spectrum, i, bg_color);

            table_widget_->setAtBottomRow(ph.getScore(), Clmn::SCORE, bg_color);
            table_widget_->setAtBottomRow((int)ph.getRank(), Clmn::RANK, bg_color);
            table_widget_->setAtBottomRow(ph.getCharge(), Clmn::CHARGE, bg_color);

            // sequence
            String seq = ph.getSequence().toString();
            if (seq.empty()) seq = ph.getMetaValue("label");
            table_widget_->setAtBottomRow(seq.toQString(), Clmn::SEQUENCE, bg_color);

            // accession
            set<String> protein_accessions = ph.extractProteinAccessionsSet();
            String accessions = ListUtils::concatenate(vector<String>(protein_accessions.begin(), protein_accessions.end()), ", ");
            table_widget_->setAtBottomRow(accessions.toQString(), Clmn::ACCESSIONS, bg_color);
            table_widget_->setAtBottomRow((int)(pi_idx), Clmn::ID_NR, bg_color);
            table_widget_->setAtBottomRow((int)(ph_idx), Clmn::PEPHIT_NR, bg_color);

            bool selected(false);
            if (ph.metaValueExists("selected"))
            {
               selected = ph.getMetaValue("selected").toString() == "true";
            }
            table_widget_->setAtBottomRow(selected, Clmn::CURATED, bg_color);

            // additional precursor infos, e.g. ppm error
            if (!precursors.empty())
            {
              const Precursor& first_precursor = precursors.front();
              double ppm_error(0);
              // Protein:RNA cross-link, Protein-Protein cross-link, or other data with a precomputed precursor error
              if (ph.metaValueExists(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM))
              {
                ppm_error = fabs((double)ph.getMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM));
              }
              else if (ph.metaValueExists("OMS:precursor_mz_error_ppm")) // for legacy reasons added in OpenMS 2.5
              {
                ppm_error = fabs((double)ph.getMetaValue("OMS:precursor_mz_error_ppm"));
              }
              else if (!ph.getSequence().empty()) // works for normal linear fragments with the correct modifications included in the AASequence
              {
                double exp_precursor = first_precursor.getMZ();
                int charge = first_precursor.getCharge();
                double theo_precursor= ph.getSequence().getMZ(charge);
                ppm_error = fabs((exp_precursor - theo_precursor) / exp_precursor / 1e-6);
              }
              table_widget_->setAtBottomRow(ppm_error, Clmn::PREC_PPM, bg_color);
            }

            // add additional meta value columns
            if (create_rows_for_commmon_metavalue_->isChecked())
            {
              Int current_col = Clmn::PEAK_ANNOTATIONS;
              // add peak annotation column (part of meta-value assessment above)
              if (has_peak_annotations)
              {
                // set hidden data for export to TSV
                QString annotation;
                for (const PeptideHit::PeakAnnotation& pa : ph.getPeakAnnotations())
                {
                  annotation += String(pa.mz).toQString() + "|" +
                    String(pa.intensity).toQString() + "|" +
                    String(pa.charge).toQString() + "|" +
                    pa.annotation.toQString() + ";";
                }
                QTableWidgetItem* item = table_widget_->setAtBottomRow("show", current_col, bg_color);
                item->setData(Qt::UserRole, annotation);
                ++current_col;
              }
              for (const auto& ck : common_keys)
              {
                DataValue dv = ph.getMetaValue(ck);
                if (dv.valueType() == DataValue::DOUBLE_VALUE)
                {
                  table_widget_->setAtBottomRow(double(dv), current_col, bg_color);
                }
                else
                {
                  table_widget_->setAtBottomRow(dv.toQString(), current_col, bg_color);
                }
                
                ++current_col;
              }
            }
          }
        }
      }

      if ((int)i == restore_spec_index)
      {
        selected_row = table_widget_->rowCount(); // get model index of selected spectrum
      }
    }

    table_widget_->setHeaders(headers);
    String s = headers.join(';');
    table_widget_->hideColumns(QStringList() << "dissociation" << "scan type" << "zoom" << "rank");
    if (has_peak_annotations) table_widget_->setHeaderExportName(Clmn::PEAK_ANNOTATIONS, "PeakAnnotations(mz|intensity|charge|annotation");

    table_widget_->resizeColumnsToContents();
    table_widget_->setSortingEnabled(true);
    table_widget_->sortByColumn(Clmn::SPEC_INDEX, Qt::AscendingOrder);

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
 
  void SpectraIDViewTab::saveIDs_()
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
    set<int> added_spectra;
    for (int r = 0; r < table_widget_->rowCount(); ++r)
    {
      // get spectrum index of current table line
      int spectrum_index = table_widget_->item(r, Clmn::SPEC_INDEX)->data(Qt::DisplayRole).toInt();

      // skip this row, if this spectrum was already processed
      if (std::find(added_spectra.begin(), added_spectra.end(), spectrum_index) != added_spectra.end())
      {
        continue;
      }
      added_spectra.insert(spectrum_index);

      // collect all PeptideIdentifications from this spectrum
      const vector<PeptideIdentification>& pep_id = (*layer_->getPeakData())[spectrum_index].getPeptideIdentifications();
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
  // update the corresponding PeptideIdentification / PeptideHits by adding a metavalue: 'selected'
  void SpectraIDViewTab::updatedSingleCell_(QTableWidgetItem* item)
  {
    // extract position of the correct Spectrum, PeptideIdentification and PeptideHit from the table
    int row = item->row();
    String selected = item->checkState() == Qt::Checked ? "true" : "false";
    int spectrum_index = table_widget_->item(row, Clmn::SPEC_INDEX)->data(Qt::DisplayRole).toInt();
    int num_id = table_widget_->item(row, Clmn::ID_NR)->data(Qt::DisplayRole).toInt();
    int num_ph = table_widget_->item(row, Clmn::PEPHIT_NR)->data(Qt::DisplayRole).toInt();

    // maintain sortability of our checkbox column
    TableView::updateCheckBoxItem(item);

    vector<PeptideIdentification>& pep_id = (*layer_->getPeakDataMuteable())[spectrum_index].getPeptideIdentifications();

    // update "selected" value in the correct PeptideHits
    vector<PeptideHit>& hits = pep_id[num_id].getHits();
    // XL-MS specific case, both PeptideHits belong to the same cross-link
    if (hits[0].metaValueExists("xl_chain")) 
    {
      hits[0].setMetaValue("selected", selected);
      if (hits.size() >= 2)
      {
        hits[1].setMetaValue("selected", selected);
      }
    }
    else // general case, update only the selected PepideHit
    {
      hits[num_ph].setMetaValue("selected", selected);
    }
  }

  void SpectraIDViewTab::fillRow_(const MSSpectrum& spectrum, const int spec_index, const QColor background_color)
  {
    const vector<Precursor>& precursors = spectrum.getPrecursors();

    table_widget_->setAtBottomRow(QString::number(spectrum.getMSLevel()), Clmn::MS_LEVEL, background_color);
    table_widget_->setAtBottomRow(spec_index, Clmn::SPEC_INDEX, background_color);
    table_widget_->setAtBottomRow(spectrum.getRT(), Clmn::RT, background_color);

    // scan mode
    table_widget_->setAtBottomRow(QString::fromStdString(spectrum.getInstrumentSettings().NamesOfScanMode[spectrum.getInstrumentSettings().getScanMode()]), Clmn::SCANTYPE, background_color);

    // zoom scan
    table_widget_->setAtBottomRow(spectrum.getInstrumentSettings().getZoomScan() ? "yes" : "no", Clmn::ZOOM, background_color);

    // fill precursor information in columns
    if (!precursors.empty())
    {
      const Precursor& first_precursor = precursors.front();

      // draw precursor information in blue
      table_widget_->setAtBottomRow(first_precursor.getMZ(), Clmn::PRECURSOR_MZ, background_color, Qt::blue);

      // set activation method
      table_widget_->setAtBottomRow(ListUtils::concatenate(first_precursor.getActivationMethodsAsString(), ",").toQString(), Clmn::DISSOCIATION, background_color);

      // set precursor intensity
      table_widget_->setAtBottomRow(first_precursor.getIntensity(), Clmn::PREC_INT, background_color);
    }
  }
}
