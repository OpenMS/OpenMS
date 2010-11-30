// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/SpectraIdentificationViewWidget.h>
#include <QtGui/QVBoxLayout>
#include <QtGui/QTreeWidget>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtGui/QHeaderView>
#include <QtGui/QMenu>
#include <QtGui/QPushButton>

#include <vector>

using namespace std;

namespace OpenMS
{
  SpectraIdentificationViewWidget::SpectraIdentificationViewWidget(const Param&, QWidget* parent)
    : QWidget(parent),
      DefaultParamHandler("SpectraIdentificationViewWidget"),
      ignore_update(false),
      layer_(0)
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
    defaults_.setValue("relative_loss_intensity", 0.1, "Relativ loss in percent");
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

    QVBoxLayout* spectra_widget_layout = new QVBoxLayout(this);
    table_widget_ = new QTableWidget(this);
    table_widget_->setObjectName("table_widget");
    table_widget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to open it.");

    table_widget_->setSortingEnabled(true);

    table_widget_->setColumnWidth(0,65);  //MS Level
    table_widget_->setColumnWidth(1,45);  //index
    table_widget_->setColumnWidth(2,70);
    table_widget_->setColumnWidth(3,70);
    table_widget_->setColumnWidth(4,55);
    table_widget_->setColumnHidden(4, true);
    table_widget_->setColumnWidth(5,45);
    table_widget_->setColumnHidden(5, true);
    table_widget_->setColumnWidth(6,45);
    table_widget_->setColumnHidden(6, true);
    table_widget_->setColumnWidth(7,45);
    table_widget_->setColumnWidth(8,45);
    table_widget_->setColumnWidth(9,45);
    table_widget_->setColumnWidth(10,400);
    table_widget_->setColumnWidth(11,45);

    ///@improvement write the visibility-status of the columns in toppview.ini and read at start

    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan type" << "zoom" << "score" << "rank" << "charge" << "          sequence          " << "accessions";
    table_widget_->setHorizontalHeaderLabels(header_labels);
    table_widget_->setColumnCount(header_labels.size());

    table_widget_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_widget_->setShowGrid(false);

    connect(table_widget_, SIGNAL(currentItemChanged(QTableWidgetItem*, QTableWidgetItem*)), this, SLOT(spectrumSelectionChange_(QTableWidgetItem*, QTableWidgetItem*)));

    spectra_widget_layout->addWidget(table_widget_);

    // additional checkboxes and buttons
    QHBoxLayout* tmp_hbox_layout = new QHBoxLayout();
    hide_ms1_ = new QCheckBox("Hide MS1", this);
    connect(hide_ms1_, SIGNAL(toggled(bool)), this, SLOT(updateEntries()));

    hide_no_identification_ = new QCheckBox("Only hits", this);
    connect(hide_no_identification_, SIGNAL(toggled(bool)), this, SLOT(updateEntries()));

    create_rows_for_commmon_metavalue_ = new QCheckBox("Extra columns for additional annotations (slow)", this);
    connect(create_rows_for_commmon_metavalue_, SIGNAL(toggled(bool)), this, SLOT(updateEntries()));

    tmp_hbox_layout->addWidget(hide_ms1_);
    tmp_hbox_layout->addWidget(hide_no_identification_);
    tmp_hbox_layout->addWidget(create_rows_for_commmon_metavalue_);
    spectra_widget_layout->addLayout(tmp_hbox_layout);
    table_widget_->sortByColumn ( 2, Qt::AscendingOrder);
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

    int spectrum_index = table_widget_->item(current->row(), 1)->data(Qt::DisplayRole).toInt() ;
    emit spectrumSelected(spectrum_index);
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
    // determine metavalues common to all hits
    if (create_rows_for_commmon_metavalue_->isChecked())
    {

      for (Size i = 0; i < layer_->getPeakData()->size(); ++i)
      {
        UInt ms_level = (*layer_->getPeakData())[i].getMSLevel();
        const vector<PeptideIdentification> peptide_ids = (*layer_->getPeakData())[i].getPeptideIdentifications();
        Size peptide_ids_count = peptide_ids.size();

        if (ms_level != 2 || peptide_ids_count == 0)  // skip non ms2 spectra and spectra with no identification
        {
          continue;    table_widget_->setRowCount(0);
        }

        for (vector<PeptideIdentification>::const_iterator pids_it = peptide_ids.begin(); pids_it != peptide_ids.end(); ++pids_it)
        {
          const vector<PeptideHit> phits = pids_it->getHits();
          for (vector<PeptideHit>::const_iterator phits_it = phits.begin(); phits_it != phits.end(); ++phits_it)
          {
            // get meta value keys
            vector<String> keys;
            phits_it->getKeys(keys);
            if (common_keys.size() == 0) // first MS2 peptide hit found. Now insert keys.
            {
              for(vector<String>::iterator sit = keys.begin(); sit != keys.end(); ++sit)
              {
                common_keys.insert(*sit);
              }
            } else // calculate intersection between current keys and common keys -> set as common_keys
            {
              set<String> current_keys;
              current_keys.insert(keys.begin(), keys.end());
              set<String> new_common_keys;
              set_intersection(current_keys.begin(), current_keys.end(), common_keys.begin(), common_keys.end(), inserter(new_common_keys, new_common_keys.begin()));
              swap(new_common_keys, common_keys);
            }
          }
        }
      }
    }

    // create header labels (setting header labels must occur after fill)
    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "dissociation"
        << "scan type" << "zoom" << "score" << "rank" << "charge" << "          sequence          " << "accession";

    for(set<String>::iterator sit = common_keys.begin(); sit != common_keys.end(); ++sit)
    {
      header_labels << sit->toQString();
    }

    table_widget_->clear();
    table_widget_->setRowCount(0);

    table_widget_->verticalHeader()->setHidden(true); // hide vertical column
    table_widget_->setColumnCount(header_labels.size());
    table_widget_->setColumnWidth(0,65);
    table_widget_->setColumnWidth(1,45);
    table_widget_->setColumnWidth(2,70);
    table_widget_->setColumnWidth(3,70);
    table_widget_->setColumnWidth(4,55);
    table_widget_->setColumnWidth(5,45);
    table_widget_->setColumnWidth(6,45);
    table_widget_->setColumnWidth(7,45);
    table_widget_->setColumnWidth(8,45);
    table_widget_->setColumnWidth(9,45);
    table_widget_->setColumnWidth(10,400);
    table_widget_->setColumnWidth(11,45);

    table_widget_->setSortingEnabled(false);
    table_widget_->setUpdatesEnabled(false);
    table_widget_->blockSignals(true);

    if (layer_ == 0)
    {
      return;
    }

    QTableWidgetItem* item = 0;
    QTableWidgetItem* selected_item = 0;

    // generate flat list
    selected_item = 0;
    for (Size i = 0; i < layer_->getPeakData()->size(); ++i)
    {
      UInt ms_level = (*layer_->getPeakData())[i].getMSLevel();

      // check whether to skip MS1 level scans
      if (ms_level == 1 && hide_ms1_->isChecked())
      {
        continue;
      }

      // hiding of entries with no identifications
      // for convinience MS1 are shown if their MS2 contain identification (if not hide MS1 is selected)
      Size id_count = (*layer_->getPeakData())[i].getPeptideIdentifications().size();
      if (hide_no_identification_->isChecked())
      {
        if (id_count == 0 && ms_level == 2)  // hide MS2 level with no identification
        {
          continue;
        } else if (ms_level == 1)
        {
          if (hide_ms1_->isChecked())
          {
            continue;
          } else // hide MS1 level if there is no MS2 with Identification
          {
            Size j = i + 1;
            bool has_ms2_with_id = false;
            while(j < layer_->getPeakData()->size())
            {
              int level = (*layer_->getPeakData())[j].getMSLevel();
              if (level == 1)
              {
                break;  // reached next MS1
              } else if (level ==2)
              {
                if ((*layer_->getPeakData())[j].getPeptideIdentifications().size() >= 1)
                {
                  has_ms2_with_id = true;
                  break;
                }
              }
              ++j;
            }

            if (!has_ms2_with_id)
            {
              continue;
            }
          }
        }
      }

      // coloring
      QColor c;
      if (ms_level == 1)
      {
        c = Qt::lightGray;  // default color for MS1
      } else if (ms_level == 2)
      {
        // get peptide identifications of current spectrum
        vector<PeptideIdentification> pi = (*layer_->getPeakData())[i].getPeptideIdentifications();
        if (pi.size() != 0)
        {
          c = Qt::green; // default color for MS1 with identification
        } else
        {
          c = Qt::white;
        }
      }

      // add new row at the end of the table
      table_widget_->insertRow(table_widget_->rowCount());

      // ms level
      item = new QTableWidgetItem(QString::number(ms_level));
      item->setBackgroundColor(c);
      item->setTextAlignment(Qt::AlignCenter);
      table_widget_->setItem(table_widget_->rowCount()-1 , 0, item);

      // index
      item = new QTableWidgetItem();
      item->setTextAlignment(Qt::AlignCenter);
      item->setData(Qt::DisplayRole, Int(i));
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount()-1 , 1, item);

      // rt
      item = new QTableWidgetItem();
      item->setTextAlignment(Qt::AlignCenter);
      item->setData(Qt::DisplayRole, (*layer_->getPeakData())[i].getRT());
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount()-1 , 2, item);

      vector<PeptideIdentification> pi = (*layer_->getPeakData())[i].getPeptideIdentifications();
      //cout << "peptide identifications: " << pi.size() << endl;
      if (pi.size()!=0)
      {
        vector<PeptideHit> ph = pi[0].getHits();

        Size best_j_index = 0;
        bool is_higher_score_better = false;
        Size best_score = pi[0].getHits()[0].getScore();
        is_higher_score_better = pi[0].isHigherScoreBetter(); // TODO: check whether its ok to assume this holds for all

        if (ph.size()!=0)
        {
          for(Size j=0; j!=pi[0].getHits().size(); ++j)
          {
            PeptideHit ph = pi[0].getHits()[j];
            // better score?
            if ((ph.getScore() < best_score && !is_higher_score_better)
              || (ph.getScore() > best_score && is_higher_score_better))
              {
              best_score = ph.getScore();
              best_j_index = j;
            }
          }
          PeptideHit best_ph = pi[0].getHits()[best_j_index];
          // score
          item = new QTableWidgetItem();
          item->setTextAlignment(Qt::AlignCenter);
          item->setData(Qt::DisplayRole, best_ph.getScore());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 7, item);

          // rank
          item = new QTableWidgetItem();
          item->setTextAlignment(Qt::AlignCenter);
          item->setData(Qt::DisplayRole, best_ph.getRank());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 8, item);

          // charge
          item = new QTableWidgetItem();
          item->setTextAlignment(Qt::AlignCenter);
          item->setData(Qt::DisplayRole, best_ph.getCharge());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 9, item);

          //sequence
          item = new QTableWidgetItem();
          item->setTextAlignment(Qt::AlignLeft);
          item->setText(best_ph.getSequence().toString().toQString());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 10, item);

          //Accession
          item = new QTableWidgetItem();
          item->setTextAlignment(Qt::AlignLeft);
          String accessions = "";

          // add protein accessions:
          for(vector<String>::const_iterator it = best_ph.getProteinAccessions().begin(); it != best_ph.getProteinAccessions().end(); ++it)
          {
            if (accessions != "")
            {
              accessions += ", " + *it;
            } else
            {
              accessions = *it;
            }
          }
          item->setText(accessions.toQString());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 11, item);

          // add additional meta value columns
          Size current_col = 12;
          if (create_rows_for_commmon_metavalue_->isChecked())
          {
            for(set<String>::iterator sit = common_keys.begin(); sit != common_keys.end(); ++sit)
            {
              DataValue dv = best_ph.getMetaValue(*sit);
              item = new QTableWidgetItem();
              item->setTextAlignment(Qt::AlignLeft);
              if (dv.valueType() == DataValue::DOUBLE_VALUE)
              {
                item->setData(Qt::DisplayRole, (DoubleReal)dv);
              } else
              {
                item->setText(dv.toQString());
              }
              item->setBackgroundColor(c);
              table_widget_->setItem(table_widget_->rowCount()-1 , current_col, item);
              ++current_col;
            }
          }
        }
      } else // no identification
      {
        // score
        item = new QTableWidgetItem();
        item->setText("-");
        item->setTextAlignment(Qt::AlignCenter);
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount()-1 , 7, item);

        // rank
        item = new QTableWidgetItem();
        item->setText("-");
        item->setTextAlignment(Qt::AlignCenter);
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount()-1 , 8, item);

        // charge
        item = new QTableWidgetItem();
        item->setText("-");
        item->setTextAlignment(Qt::AlignCenter);
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount()-1 , 9, item);

        //sequence
        item = new QTableWidgetItem();
        item->setText("-");
        item->setTextAlignment(Qt::AlignCenter);
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount()-1 , 10, item);

        //accession
        item = new QTableWidgetItem();
        item->setText("-");
        item->setTextAlignment(Qt::AlignCenter);
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount()-1 , 11, item);

        // add additional meta value columns
        Size current_col = 12;
        if (create_rows_for_commmon_metavalue_->isChecked())
        {
          for(set<String>::iterator sit = common_keys.begin(); sit != common_keys.end(); ++sit)
          {
            item = new QTableWidgetItem();
            item->setTextAlignment(Qt::AlignLeft);
            item->setText("-");
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount()-1 , current_col, item);
            ++current_col;
          }
        }
      }

      if (!(*layer_->getPeakData())[i].getPrecursors().empty())  // has precursor
      {
        item = new QTableWidgetItem();
        item->setTextAlignment(Qt::AlignCenter);
        item->setData(Qt::DisplayRole, (*layer_->getPeakData())[i].getPrecursors()[0].getMZ());
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount()-1 , 3, item);

        item = new QTableWidgetItem();
        item->setTextAlignment(Qt::AlignCenter);
        if (!(*layer_->getPeakData())[i].getPrecursors().front().getActivationMethods().empty())
        {
          QString t;
          for(std::set<Precursor::ActivationMethod>::const_iterator it = (*layer_->getPeakData())[i].getPrecursors().front().getActivationMethods().begin(); it != (*layer_->getPeakData())[i].getPrecursors().front().getActivationMethods().end(); ++it)
          {
            if(!t.isEmpty()){
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
        table_widget_->setItem(table_widget_->rowCount()-1 , 4, item);
      }
      else  // has no precursor (leave fields 3 and 4 empty)
      {
        item = new QTableWidgetItem();
        item->setTextAlignment(Qt::AlignCenter);
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount()-1 , 3, item);
        item = new QTableWidgetItem();
        item->setTextAlignment(Qt::AlignCenter);
        item->setText("-");
        item->setBackgroundColor(c);
        table_widget_->setItem(table_widget_->rowCount()-1 , 4, item);
      }

      // scan mode
      item = new QTableWidgetItem();
      item->setTextAlignment(Qt::AlignCenter);
      if ((*layer_->getPeakData())[i].getInstrumentSettings().getScanMode()>0)
      {
        item->setText(QString::fromStdString((*layer_->getPeakData())[i].getInstrumentSettings().NamesOfScanMode[(*layer_->getPeakData())[i].getInstrumentSettings().getScanMode()]));
      }
      else
      {
        item->setText("-");
      }
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount()-1 , 5, item);

      // zoom scan
      item = new QTableWidgetItem();
      item->setTextAlignment(Qt::AlignCenter);
      if ((*layer_->getPeakData())[i].getInstrumentSettings().getZoomScan())
      {
        item->setText("yes");
      }
      else
      {
        item->setText("no");
      }
      item->setBackgroundColor(c);
      table_widget_->setItem(table_widget_->rowCount()-1 , 6, item);

      if (i == layer_->current_spectrum)
      {
        // just remember it, select later
        selected_item = item;
      }
    }

    if (selected_item)
    {
      // now, select and scroll down to item
      selected_item->setSelected(true);
      table_widget_->scrollToItem(selected_item);
    }

    table_widget_->blockSignals(false);
    table_widget_->setUpdatesEnabled(true);
    table_widget_->setSortingEnabled(true);
    table_widget_->setHorizontalHeaderLabels(header_labels);
    table_widget_->resizeColumnsToContents();
    table_widget_->resizeRowsToContents();
  }

  SpectraIdentificationViewWidget::~SpectraIdentificationViewWidget()
  {
  }

}
