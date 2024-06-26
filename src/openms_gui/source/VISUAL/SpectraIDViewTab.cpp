// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "OpenMS/VISUAL/LayerDataPeak.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/SYSTEM/NetworkGetRequest.h>
#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/SequenceVisualizer.h>
#include <OpenMS/VISUAL/SpectraIDViewTab.h>
#include <OpenMS/VISUAL/TableView.h>

#include <QJsonArray>
#include <QJsonObject>
#include <QJsonValue>
#include <QRegularExpression>
#include <QString>
#include <QStringList>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>

#include <vector>
#include <string>

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

// Use a namespace to encapsulate names, yet use c-style 'enum' for fast conversion to int.
// So we can write: 'Clmn::MS_LEVEL', but get implicit conversion to int
namespace ProteinClmn
{
  enum HeaderNames
      { // indices into QTableWidget's columns (which start at index 0)
    ACCESSION,
    FULL_PROTEIN_SEQUENCE,
    SEQUENCE,
    DESCRIPTION,
    SCORE,
    COVERAGE,
    NR_PSM,
    /* last entry --> */ SIZE_OF_HEADERNAMES
      };
  // keep in SYNC with enum HeaderNames
  const QStringList HEADER_NAMES = QStringList()
      << "accession" << "full sequence" << "sequence" << "description" << "score" << "coverage" << "#PSMs";
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

    // id view parameters (warning: must be matched in TOPPViewPrefDialog)
    defaults_.insert("tsg:", TheoreticalSpectrumGenerator().getParameters());
    defaults_.insert("align:", SpectrumAlignment().getParameters());

    QVBoxLayout* all = new QVBoxLayout(this);
    tables_splitter_ = new QSplitter(Qt::Horizontal);

    table_widget_ = new TableView(tables_splitter_);

    // exported protein accessions and PSM rank even if hidden
    table_widget_->setMandatoryExportColumns(QStringList() << "accessions" << "rank");
    
    table_widget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to open it.");
    tables_splitter_->addWidget(table_widget_);

    protein_table_widget_ = new TableView(tables_splitter_);
    protein_table_widget_->setWhatsThis("Protein selection bar<BR><BR>Here all proteins of the current experiment are shown. TODO what can you do with it");

    tables_splitter_->addWidget(protein_table_widget_);

    all->addWidget(tables_splitter_);
    
    ////////////////////////////////////
    // additional checkboxes and buttons
    QHBoxLayout* buttons_hbox_layout = new QHBoxLayout();

    hide_no_identification_ = new QCheckBox("Only hits", this);
    hide_no_identification_->setChecked(true);

    create_rows_for_commmon_metavalue_ = new QCheckBox("Show advanced\nannotations", this);

    QPushButton* save_IDs = new QPushButton("Save IDs", this);
    connect(save_IDs, &QPushButton::clicked, this, &SpectraIDViewTab::saveIDs_);

    QPushButton* export_table = new QPushButton("Export table", this);

    QPushButton* switch_orientation_ = new QPushButton("Switch orientation", this);
    connect(switch_orientation_, &QPushButton::clicked, this, &SpectraIDViewTab::switchOrientation_);

    buttons_hbox_layout->addWidget(hide_no_identification_);
    buttons_hbox_layout->addWidget(create_rows_for_commmon_metavalue_);
    buttons_hbox_layout->addWidget(save_IDs);
    buttons_hbox_layout->addWidget(export_table);
    buttons_hbox_layout->addWidget(switch_orientation_);
    all->addLayout(buttons_hbox_layout);
    //TODO add search boxes like in spectrum list view
    //TODO add score filter box or Apply Tool to identifications

    connect(table_widget_, &QTableWidget::currentCellChanged, this, &SpectraIDViewTab::currentCellChanged_);
    connect(table_widget_, &QTableWidget::itemChanged, this, &SpectraIDViewTab::updatedSingleCell_);
    connect(table_widget_->selectionModel(), &QItemSelectionModel::selectionChanged, this, &SpectraIDViewTab::currentSpectraSelectionChanged_);
    connect(protein_table_widget_, &QTableWidget::cellClicked, this, &SpectraIDViewTab::proteinCellClicked_);
    connect(protein_table_widget_, &QTableWidget::itemChanged, this, &SpectraIDViewTab::updatedSingleProteinCell_);
    connect(hide_no_identification_, &QCheckBox::toggled, this, &SpectraIDViewTab::updateEntries_);
    connect(create_rows_for_commmon_metavalue_, &QCheckBox::toggled, this, &SpectraIDViewTab::updateEntries_);
    connect(export_table, &QPushButton::clicked, table_widget_, &TableView::exportEntries);
  }

  void SpectraIDViewTab::clear()
  {
    table_widget_->clear();
    protein_table_widget_->clear();
    layer_ = nullptr;
  }

  // Create the protein accession to peptide identification map using C++ STL unordered_map
  void SpectraIDViewTab::createProteinToPeptideIDMap_()
  {
    //clear the map each time entries are updated with updateEntries()
    protein_to_peptide_id_map.clear();

    if (is_first_time_loading_ && layer_)
    {
      for (const auto& spec : *layer_->getPeakData())
      {
        if (!spec.getPeptideIdentifications().empty())
        {
          const vector<PeptideIdentification>& peptide_ids = spec.getPeptideIdentifications();

          for (const auto& pepid : peptide_ids)
          {
            const vector<PeptideHit>& pep_hits = pepid.getHits();
            //add id_accession as the key of the map and push the peptideID to the vector value-
            for (const auto & pep_hit : pep_hits)
            {
              const vector<PeptideEvidence>& evidences = pep_hit.getPeptideEvidences();

              for (const auto & evidence : evidences)
              {
                const String& id_accession = evidence.getProteinAccession();
                protein_to_peptide_id_map[id_accession].push_back(&pepid);
              }
            }
          }
        }
      }
      // set is_first_time_loading to false so that the map gets created only the first time!
      is_first_time_loading_ = false;
    }
  }

  //extract required part of accession and open browser
  QString SpectraIDViewTab::extractNumFromAccession_(const QString& full_accession)
  {
    // anchored (^...$) regex for matching accession
    QRegularExpression reg_pre_accession("^(tr|sp)$", QRegularExpression::PatternOption::CaseInsensitiveOption);
    QRegularExpression reg_uniprot_accession("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$");

    // The full accession is in the form "tr|A9GID7|A9GID7_SORC5" or "P02769|ALBU_BOVIN", 
    // so split it with | and get the individual parts
    QStringList acsn = full_accession.split("|");

    for (const QString& substr : acsn)
    {
      //eg, substr2 = tr, substr2 = p02769 etc
      // if substr = tr/sp then skip
      if (reg_pre_accession.match(substr.simplified()).hasMatch())
      {
        continue;
      }
      else
      {
        if (reg_uniprot_accession.match(substr.simplified()).hasMatch())
        {
          return substr.simplified();
        }
        else
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid accession found!", 
              String(full_accession));
        }
      }
    }
    return {};
  }

  void SpectraIDViewTab::openUniProtSiteWithAccession_(const QString& accession)
  {
    QString accession_num;
    try
    {
      accession_num = extractNumFromAccession_(accession);
    }
    catch (Exception::InvalidValue&)
    {
      // TODO: print in status(?) that accession format is not supported
    }

    if (!accession_num.isEmpty()) 
    {
      QString base_url = "https://www.uniprot.org/uniprot/";
      QString url = base_url + accession_num;
      GUIHelpers::openURL(url);
    }
  }

  void SpectraIDViewTab::proteinCellClicked_(int row, int column)
  {
    //TODO maybe highlight/filter all PepHits that may provide evidence for this protein (or at least that are top scorer)
    if (row < 0 || column < 0)
      return;

    if (row >= protein_table_widget_->rowCount() || column >= protein_table_widget_->columnCount())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid cell clicked.", String(row) + " " + column);
    }

    // Open browser with accession when clicked on the accession column on a row
    if (column == ProteinClmn::ACCESSION)
    {
      // This stores the complete accession, eg, "tr|A9GID7|A9GID7_SORC5"
      QString accession = protein_table_widget_->item(row, ProteinClmn::ACCESSION)->data(Qt::DisplayRole).toString();
      // As with the current logic, we have only one accession per row, we can directly use that accession 
      // while opening the window instead of showing another widget that lists all accessions
      openUniProtSiteWithAccession_(accession);
    }

    //
    // Check if Qt WebEngineWidgets is installed on user's machine and if so,
    // open a new window to visualize protein sequence
    #ifdef QT_WEBENGINEWIDGETS_LIB
    if (column == ProteinClmn::SEQUENCE)
    {
      // store the current sequence clicked from the FULL_PROTEIN_SEQUENCE column. This column(hidden by default) 
      // stores the full protein sequence
      QString protein_sequence = protein_table_widget_->item(row, ProteinClmn::FULL_PROTEIN_SEQUENCE)->data(Qt::DisplayRole).toString();
      // store the accession as string, eg: tr|P02769|ALBU_BOVIN
      QString current_accession = protein_table_widget_->item(row, ProteinClmn::ACCESSION)->data(Qt::DisplayRole).toString();

      // extract the part of accession , eg: P02769
      QString accession_num;
      try
      {
        accession_num = extractNumFromAccession_(current_accession);
      }
      catch (Exception::InvalidValue&)
      {
        // TODO: print in status(?) that accession format is not supported
      }    

      auto item_pepid = table_widget_->item(row, Clmn::ID_NR);

      if (item_pepid)
      {

        //array to store object of start-end positions, sequence and mod data of peptides;
        QJsonArray peptides_data;
       
        //use data from the protein_to_peptide_id_map map and store the start/end position to the QJsonArray
        for (auto pep_id_ptr : protein_to_peptide_id_map[current_accession])
        {
          const vector<PeptideHit>& pep_hits = pep_id_ptr->getHits();

          //store start and end positions
          //TODO maybe we could store the index of the hit that belongs to that specific protein in the map as well
          // or we generally should only look at the first hit
          for (const auto & pep_hit : pep_hits)
          {
            const vector<PeptideEvidence>& evidences = pep_hit.getPeptideEvidences();
            const AASequence& aaseq = pep_hit.getSequence();
            const auto qstrseq = aaseq.toString().toQString();

            for (const auto & evidence : evidences)
            {
              const String& id_accession = evidence.getProteinAccession();
              QJsonObject pep_data_obj;
              int pep_start = evidence.getStart();
              int pep_end = evidence.getEnd();
              if (id_accession.toQString() == current_accession)
              {
                // contains key-value of modName and vector of indices
                QJsonObject mod_data;

                for (int i = 0; i < (int)aaseq.size(); ++i)
                {
                  if (aaseq[i].isModified())
                  {
                    const String& mod_name = aaseq[i].getModificationName();

                    if (!mod_data.contains(mod_name.toQString()))
                    {
                      mod_data[mod_name.toQString()] = QJsonArray{i + pep_start}; // add pep_start to get the correct location in the whole sequence
                    }
                    else
                    {
                      QJsonArray values = mod_data.value(mod_name.toQString()).toArray();
                      // add pep_start to get the correct location in the whole sequence
                      values.push_back(i + pep_start); 
                      mod_data[mod_name.toQString()] = values;
                    }
                  }
                }
                pep_data_obj["start"] = pep_start;
                pep_data_obj["end"] = pep_end;
                pep_data_obj["seq"] = qstrseq;
                pep_data_obj["mod_data"] = mod_data;
                //Push objects to array that will be passed to html
                peptides_data.push_back(pep_data_obj);
              }
            }
          }
        }

        auto* widget = new SequenceVisualizer(this); // no parent since we want a new window
        widget->setWindowFlags(Qt::Window);
        widget->resize(1500,500); // make a bit bigger
        widget->setProteinPeptideDataToJsonObj(accession_num, protein_sequence, peptides_data);
        widget->show();
      }
    }
    #endif
  }

  void SpectraIDViewTab::currentSpectraSelectionChanged_()
  {
    if (table_widget_->selectionModel()->selectedRows().empty())
    {
      // deselect whatever is currently shown
      //layer_->getCurrentSpectrumIndex();
      // Deselecting spectrum does not do what you think it does. It still paints stuff. Without annotations..
      // so just leave it for now.
      //
      // PARTLY SOLVED: The problem was, that if you defocus the TOPPView window, somehow
      // selectionChange is called, with EMPTY selection. Maybe this is a feature and we have to store the
      // selected spectrum indices as well. I want to support multi-selection in the future to see shared peptides
      // Actually this might be solved by the removal of the unnecessary updates in activateSubWindow.
      // I think updateEntries resets selections as well.. not sure how we could avoid that. We really have to avoid
      // calling this crazy function when only small updates are needed.
      //emit spectrumDeselected(last_spectrum_index);
      // TODO also currently, the current active spectrum can be restored after deselection by clicking on
      //  the Scans tab and then switching back to ID tab. (Scans will get the current scan in the 1D View, which
      //  is still there. I guess I have to deselect in the 1D view, too, after all.
      updateProteinEntries_(-1);
    }
    //TODO if you deselected the current spectrum, you currently cannot click on/navigate to the same spectrum
    // because currentCellChanged_ will not trigger. We would need to do it here.
  }

  void SpectraIDViewTab::currentCellChanged_(int row, int column, int /*old_row*/, int /*old_column*/)
  {
    // TODO you actually only have to do repainting if the row changes..
    // sometimes Qt calls this function when table empty during refreshing
    if (row < 0 || column < 0)
    {
      return;
    }

    if (row >= table_widget_->rowCount()
        ||  column >= table_widget_->columnCount())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid cell clicked.", String(row) + " " + column);
    }
    
    // deselect whatever is currently shown (if we are in 1D view)
    auto* layer_1d = dynamic_cast<LayerData1DPeak*>(layer_);
    if (layer_1d)
    {
      emit spectrumDeselected(int(layer_1d->getCurrentIndex()));
    }

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

        emit spectrumSelected(std::distance(exp.begin(), prec_it), -1, -1);// no identification or hit selected (-1)
        // zoom into precursor area
        emit requestVisibleArea1D(isolation_window_lower_mz - 50.0, isolation_window_upper_mz + 50.0);
      }
    }
    else
    {// if spectrum with no PepIDs is selected, there is nothing to show...
      auto item_pepid = table_widget_->item(row, Clmn::ID_NR);
      if (item_pepid == nullptr// null for MS1 spectra
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
      if (item_pepid)// might be null for MS1 spectra
      {
        int current_identification_index = item_pepid->data(Qt::DisplayRole).toInt();
        int current_peptide_hit_index = table_widget_->item(row, Clmn::PEPHIT_NR)->data(Qt::DisplayRole).toInt();

        const vector<PeptideIdentification>& peptide_ids = spec2.getPeptideIdentifications();
        const vector<PeptideHit>& pep_hits = peptide_ids[current_identification_index].getHits();
        const PeptideHit& hit = pep_hits[current_peptide_hit_index];

        // initialize window, when the table is requested for the first time
        // afterwards the size will stay at the manually resized window size
        if (fragment_window_ == nullptr)
        {
          fragment_window_ = new QTableWidget();
          fragment_window_->resize(320, 500);

          fragment_window_->verticalHeader()->setHidden(true);// hide vertical column

          QStringList header_labels;
          header_labels << "m/z"
                        << "name"
                        << "intensity"
                        << "charge";
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

        for (const PeptideHit::PeakAnnotation& pa : hit.getPeakAnnotations())
        {
          fragment_window_->insertRow(fragment_window_->rowCount());
          QTableWidgetItem* item = fragment_window_->itemPrototype()->clone();
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

    // Update the protein table with data of the id row that was clicked
    updateProteinEntries_(row);
  }

  bool SpectraIDViewTab::hasData(const LayerDataBase* layer)
  {
    // this is a very easy check.
    // We do not check for PeptideIdentifications attached to Spectra, because the user could just
    // want the list of unidentified MS2 spectra (obtained by unchecking the 'just hits' button).
    auto* ptr_peak = dynamic_cast<const LayerDataPeak*>(layer);
    bool no_data = (ptr_peak == nullptr
                    || (ptr_peak && ptr_peak->getPeakData()->empty()));
    return !no_data;
  }

  void SpectraIDViewTab::updateEntries(LayerDataBase* cl)
  {
    // do not try to be smart and check if layer_ == cl; to return early
    // since the layer content might have changed, e.g. pepIDs were added
    auto* ptr_peak = dynamic_cast<LayerDataPeak*>(cl);
    layer_ = ptr_peak; // might be nullptr

    // setting "is_first_time_loading_ = true;" here currently negates the logic of creating the map only the first time
    // the data loads, but in future, after fixing the issue of calling updateEntries() multiple times, we can use it to only
    // create the map when the table data loads completely new data from idXML file. Currently the map gets created each time 
    // the updateEntries() is called.
    is_first_time_loading_ = true;
    createProteinToPeptideIDMap_();
    updateEntries_(); // we need this extra function since it's an internal slot
  }

  LayerDataBase* SpectraIDViewTab::getLayer()
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
  }// namespace Detail

  void SpectraIDViewTab::updateProteinEntries_(int selected_spec_row_idx)
  {
    //TODO Currently when switching to 2D view of the same dataset and then switching back to the fragment spectrum,
    // the spectrum table (almost; annotations gone) correctly restores the row, while the proteins do not get newly
    // refreshed. Check why and fix. It is not too bad though.
    // no valid peak layer attached
    if (!hasData(layer_) || layer_->getPeakData()->getProteinIdentifications().empty())
    {
      //clear(); this was done in updateEntries_() already.
      return;
    }

    if (ignore_update)
    {
      return;
    }

    if (!isVisible())
    {
      return;
    }

    set<String> accs;
    if(selected_spec_row_idx >= 0)
      //TODO another option would be a "Filter proteins" checkbox that filters for proteins for this Hit
      // only when checked, otherwise only highlights
    {
      int row = selected_spec_row_idx;
      int spectrum_index = table_widget_->item(row, Clmn::SPEC_INDEX)->data(Qt::DisplayRole).toInt();
      int num_id = table_widget_->item(row, Clmn::ID_NR)->data(Qt::DisplayRole).toInt();
      int num_ph = table_widget_->item(row, Clmn::PEPHIT_NR)->data(Qt::DisplayRole).toInt();
      const auto& spec = layer_->getPeakData()->operator[](spectrum_index);
      const vector<PeptideIdentification>& pep_id = spec.getPeptideIdentifications();

      if(!spec.getPeptideIdentifications().empty())
      {
        const vector<PeptideHit>& hits = pep_id[num_id].getHits();
        if (!hits.empty()) accs = hits[num_ph].extractProteinAccessionsSet();
      }
    }

    // create header labels (setting header labels must occur after fill)
    QStringList headers = ProteinClmn::HEADER_NAMES;

    protein_table_widget_->clear();
    protein_table_widget_->setRowCount(0);
    protein_table_widget_->setColumnCount(headers.size());
    protein_table_widget_->setSortingEnabled(false);
    protein_table_widget_->setUpdatesEnabled(false);
    protein_table_widget_->blockSignals(true);

    // generate flat list
    int selected_row(-1);
    // index i is needed, so iterate the old way...
    for (Size i = 0; i < layer_->getPeakData()->getProteinIdentifications()[0].getHits().size(); ++i)
    {
      const auto& protein = layer_->getPeakData()->getProteinIdentifications()[0].getHits()[i];
      if (accs.empty() || accs.find(protein.getAccession()) != accs.end())
      {
        // set row background color
        QColor bg_color = accs.empty() ? Qt::white : Qt::lightGray;

        int total_pepids = protein_to_peptide_id_map[protein.getAccession()].size();
        
        // add new row at the end of the table
        protein_table_widget_->insertRow(protein_table_widget_->rowCount());

        protein_table_widget_->setAtBottomRow(protein.getAccession().toQString(), ProteinClmn::ACCESSION, bg_color, Qt::blue);
        protein_table_widget_->setAtBottomRow(protein.getSequence().toQString(), ProteinClmn::FULL_PROTEIN_SEQUENCE, bg_color);
        protein_table_widget_->setAtBottomRow("show", ProteinClmn::SEQUENCE, bg_color, Qt::blue);
        protein_table_widget_->setAtBottomRow(protein.getDescription().toQString(), ProteinClmn::DESCRIPTION, bg_color);
        protein_table_widget_->setAtBottomRow(protein.getScore(), ProteinClmn::SCORE, bg_color);
        protein_table_widget_->setAtBottomRow(protein.getCoverage(), ProteinClmn::COVERAGE, bg_color);
        protein_table_widget_->setAtBottomRow(total_pepids, ProteinClmn::NR_PSM, bg_color);

        /*if ((int)i == restore_spec_index) //TODO actually extract the accessions for the selected spectrum and compare
        {
          selected_row = protein_table_widget_->rowCount() - 1; // get model index of selected spectrum
        }*/
      }
    }

    protein_table_widget_->setHeaders(headers);
    protein_table_widget_->setColumnHidden(ProteinClmn::FULL_PROTEIN_SEQUENCE, true);
    #ifndef QT_WEBENGINEWIDGETS_LIB
      protein_table_widget_->setColumnHidden(ProteinClmn::SEQUENCE, true); // no web engine? hide sequence column used to do the JS query
    #endif
    protein_table_widget_->resizeColumnsToContents();
    protein_table_widget_->setSortingEnabled(true);
    protein_table_widget_->sortByColumn(ProteinClmn::SCORE, Qt::AscendingOrder); //TODO figure out higher_score_better

    if (selected_row != -1)  // select and scroll down to item
    {
      protein_table_widget_->selectRow(selected_row);
      QTableWidgetItem* selected_item = protein_table_widget_->item(selected_row, 0);
      selected_item->setSelected(true);
      protein_table_widget_->setCurrentItem(selected_item);
      protein_table_widget_->scrollToItem(selected_item);
    }

    protein_table_widget_->blockSignals(false);
    protein_table_widget_->setUpdatesEnabled(true);
    protein_table_widget_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    protein_table_widget_->verticalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  }

  void SpectraIDViewTab::updateEntries_()
  {

    // no valid peak layer attached
    if (!hasData(layer_))
    {
      clear();
      return;
    }

    if (ignore_update)
    { 
      return; 
    }

    if (!isVisible())
    { 
      return;
    }

    auto layer_peak = dynamic_cast<LayerData1DPeak*>(layer_);
    if (!layer_peak) return;

    int restore_spec_index = int(layer_peak->getCurrentIndex());

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

        if (ms_level != 2 || peptide_ids.empty()) // skip non ms2 spectra and spectra with no identification
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

    table_widget_->blockSignals(true); // to be safe, that clear does not trigger anything.
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
      if (ms_level != 2 && id_count == 0)
      { 
        continue;
      }

      // skip
      if (hide_no_identification_->isChecked() && id_count == 0) 
      { 
        continue;
      }
      // set row background color
      QColor bg_color = (id_count == 0 ? Qt::white : QColor::fromRgb(127,255,148));

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
            if (seq.empty())
            {
              seq = ph.getMetaValue("label");
            }
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
                QTableWidgetItem* item = table_widget_->setAtBottomRow("show", current_col, bg_color, Qt::blue);
                item->setData(Qt::UserRole, annotation);
                ++current_col;
              }
              for (const auto& ck : common_keys)
              {
                const DataValue& dv = ph.getMetaValue(ck);
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
        // get model index of selected spectrum, 
        // as table_widget_->rowCount() returns rows starting from 1, selected row is 1 less than the returned row
        selected_row = table_widget_->rowCount() - 1; 
      }
    }

    table_widget_->setHeaders(headers);
    String s = headers.join(';');
    table_widget_->hideColumns(QStringList() << "accessions"
                                             << "dissociation"
                                             << "scan type"
                                             << "zoom"
                                             << "rank"
                                             << "#ID"
                                             << "#PH");
    if (has_peak_annotations) table_widget_->setHeaderExportName(Clmn::PEAK_ANNOTATIONS, "PeakAnnotations(mz|intensity|charge|annotation");

    table_widget_->setSortingEnabled(true);
    table_widget_->sortByColumn(Clmn::SPEC_INDEX, Qt::AscendingOrder);

    if (selected_row != -1)  // select and scroll down to item
    {
      table_widget_->selectRow(selected_row);
      QTableWidgetItem* selected_item = table_widget_->item(selected_row, 0);
      selected_item->setSelected(true);
      table_widget_->setCurrentItem(selected_item);
      table_widget_->scrollToItem(selected_item);
      currentCellChanged_(selected_row, 0, 0, 0); // simulate cell change to trigger repaint and reannotation of spectrum 1D view
    }

    table_widget_->verticalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    table_widget_->blockSignals(false);
    table_widget_->setUpdatesEnabled(true);

    // call this updateProteinEntries_(-1) function after the table_widget data is filled, 
    // otherwise table_widget_->item(row, clm) returns nullptr;
    updateProteinEntries_(selected_row);
  }

  void SpectraIDViewTab::switchOrientation_()
  {
    if (tables_splitter_->orientation() == Qt::Vertical) 
    {
      tables_splitter_->setOrientation(Qt::Horizontal);
    }
    else
    {
      tables_splitter_->setOrientation(Qt::Vertical);
    }

  }
 
  void SpectraIDViewTab::saveIDs_()
  {
    // no valid peak layer attached
    if (layer_ == nullptr || layer_->getPeakData()->empty() || layer_->type != LayerDataBase::DT_PEAK)
    {
      return;
    }

    // synchronize PeptideHits with the annotations in the spectrum
    dynamic_cast<LayerData1DPeak*>(layer_)->synchronizePeakAnnotations();

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
      if (added_spectra.find(spectrum_index) != added_spectra.end())
      {
        continue;
      }
      added_spectra.insert(spectrum_index);

      // collect all PeptideIdentifications from this spectrum
      const vector<PeptideIdentification>& pep_id = (*layer_->getPeakData())[spectrum_index].getPeptideIdentifications();
      copy(pep_id.begin(), pep_id.end(), back_inserter(all_pep_ids));
    }

    QString filename = GUIHelpers::getSaveFilename(this, "Save file", "", FileTypeList({FileTypes::IDXML, FileTypes::MZIDENTML}), true, FileTypes::IDXML);
    if (filename.isEmpty())
    {
      return;
    }      
    FileHandler().storeIdentifications(filename, prot_id, all_pep_ids, {FileTypes::IDXML, FileTypes::MZIDENTML});
  }

  void SpectraIDViewTab::updatedSingleProteinCell_(QTableWidgetItem* /*item*/)
  {    
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
    else // general case, update only the selected PeptideHit
    {
      hits[num_ph].setMetaValue("selected", selected);
    }
  }

  void SpectraIDViewTab::fillRow_(const MSSpectrum& spectrum, const int spec_index, const QColor& background_color)
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

  void SpectraIDViewTab::SelfResizingTableView_::resizeEvent(QResizeEvent * /*event*/)
  {

  }
}
