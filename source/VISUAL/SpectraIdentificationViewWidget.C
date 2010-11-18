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

#include <vector>

using namespace std;

namespace OpenMS
{
  SpectraIdentificationViewWidget::SpectraIdentificationViewWidget(const Param&, QWidget* parent)
    : QWidget(parent),
      DefaultParamHandler("SpectraIdentificationViewWidget"),
      ignore_update(false)
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

    ///@improvement write the visibility-status of the columns in toppview.ini and read at start

    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan type" << "zoom" << "score" << "rank" << "charge" << "sequence";
    table_widget_->setHorizontalHeaderLabels(header_labels);
    table_widget_->setColumnCount(header_labels.size());

    table_widget_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_widget_->setShowGrid(false);

    connect(table_widget_, SIGNAL(currentItemChanged(QTableWidgetItem*, QTableWidgetItem*)), this, SLOT(spectrumSelectionChange_(QTableWidgetItem*, QTableWidgetItem*)));

    spectra_widget_layout->addWidget(table_widget_);

    QHBoxLayout* tmp_hbox_layout = new QHBoxLayout();

    spectra_search_box_ = new QLineEdit("", this);

    QStringList qsl;
    qsl.push_back("index");
    qsl.push_back("RT");
    qsl.push_back("MZ");
    qsl.push_back("dissociation");
    qsl.push_back("scan");
    qsl.push_back("zoom");
    spectra_combo_box_ = new QComboBox(this);
    spectra_combo_box_->addItems(qsl);

    connect(spectra_search_box_, SIGNAL(textEdited ( const QString &)), this, SLOT(spectrumSelected_(const QString&)));

    tmp_hbox_layout->addWidget(spectra_search_box_);
    tmp_hbox_layout->addWidget(spectra_combo_box_);
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

  void SpectraIdentificationViewWidget::updateEntries(const LayerData& cl)
  {
    if (ignore_update)
    {
      return;
    }
    if (!table_widget_->isVisible())
    {
      return;
    }

    table_widget_->setSortingEnabled(false);
    table_widget_->blockSignals(true);
    table_widget_->clear();
    table_widget_->setRowCount(0);

    QStringList header_labels;
    header_labels << "MS" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan type" << "zoom" << "score" << "rank" << "charge" << "sequence";
    table_widget_->setHorizontalHeaderLabels(header_labels);
    table_widget_->verticalHeader()->setHidden(true); // hide vertical column
    table_widget_->setColumnCount(header_labels.size());
    table_widget_->setColumnWidth(0,65);  //MS Level
    table_widget_->setColumnWidth(1,45);  //index
    table_widget_->setColumnWidth(2,70);
    table_widget_->setColumnWidth(3,70);
    table_widget_->setColumnWidth(4,55);
    table_widget_->setColumnWidth(5,45);
    table_widget_->setColumnWidth(6,45);
    table_widget_->setColumnWidth(7,45);

    QTableWidgetItem* item = 0;
    QTableWidgetItem* selected_item = 0;

    if(cl.type == LayerData::DT_PEAK)
    {
        // generate flat list
        selected_item = 0;
        for (Size i = 0; i < cl.getPeakData()->size(); ++i)
        {
          UInt ms_level = (*cl.getPeakData())[i].getMSLevel();

          // coloring
          QColor c;
          if (ms_level == 1)
          {
            c = Qt::lightGray;  // default color for MS1
          } else if (ms_level == 2)
          {
            // get peptide identifications of current spectrum
            vector<PeptideIdentification> pi = (*cl.getPeakData())[i].getPeptideIdentifications();
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
          table_widget_->setItem(table_widget_->rowCount()-1 , 0, item);

          // index
          item = new QTableWidgetItem();          
          item->setData(Qt::DisplayRole, i);
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 1, item);

          // rt
          item = new QTableWidgetItem();
          item->setData(Qt::DisplayRole, (*cl.getPeakData())[i].getRT());
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 2, item);

          vector<PeptideIdentification> pi = (*cl.getPeakData())[i].getPeptideIdentifications();
          //cout << "peptide identifications: " << pi.size() << endl;
          if (pi.size()!=0)
          {
            vector<PeptideHit> ph = pi[0].getHits();
            if (ph.size()!=0)  // @TODO: select best scoring hit
            {
              // score
              item = new QTableWidgetItem();
              item->setData(Qt::DisplayRole, ph[0].getScore());
              item->setBackgroundColor(c);
              table_widget_->setItem(table_widget_->rowCount()-1 , 7, item);

              // rank
              item = new QTableWidgetItem();
              item->setData(Qt::DisplayRole, ph[0].getRank());
              item->setBackgroundColor(c);
              table_widget_->setItem(table_widget_->rowCount()-1 , 8, item);

              // charge
              item = new QTableWidgetItem();
              item->setData(Qt::DisplayRole, ph[0].getCharge());
              item->setBackgroundColor(c);
              table_widget_->setItem(table_widget_->rowCount()-1 , 9, item);

              //sequence
              item = new QTableWidgetItem();
              item->setText(ph[0].getSequence().toString().toQString());
              item->setBackgroundColor(c);
              table_widget_->setItem(table_widget_->rowCount()-1 , 10, item);
            }
          } else
          {
            // score
            item = new QTableWidgetItem();
            item->setText("-");
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount()-1 , 7, item);

            // rank
            item = new QTableWidgetItem();
            item->setText("-");
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount()-1 , 8, item);

            // charge
            item = new QTableWidgetItem();
            item->setText("-");
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount()-1 , 9, item);

            //sequence
            item = new QTableWidgetItem();
            item->setText("-");
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount()-1 , 10, item);
          }

          if (!(*cl.getPeakData())[i].getPrecursors().empty())  // has precursor
          {            
            item = new QTableWidgetItem();
            item->setData(Qt::DisplayRole, (*cl.getPeakData())[i].getPrecursors()[0].getMZ());
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount()-1 , 3, item);

            item = new QTableWidgetItem();
            if (!(*cl.getPeakData())[i].getPrecursors().front().getActivationMethods().empty())
            {
              QString t;
              for(std::set<Precursor::ActivationMethod>::const_iterator it = (*cl.getPeakData())[i].getPrecursors().front().getActivationMethods().begin(); it != (*cl.getPeakData())[i].getPrecursors().front().getActivationMethods().end(); ++it)
              {
                if(!t.isEmpty()){
                  t.append(",");
                }
                t.append(QString::fromStdString((*cl.getPeakData())[i].getPrecursors().front().NamesOfActivationMethod[*((*cl.getPeakData())[i].getPrecursors().front().getActivationMethods().begin())]));
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
            item->setText("-");
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount()-1 , 3, item);
            item = new QTableWidgetItem();
            item->setText("-");
            item->setBackgroundColor(c);
            table_widget_->setItem(table_widget_->rowCount()-1 , 4, item);
          }

          // scan mode
          item = new QTableWidgetItem();
          if ((*cl.getPeakData())[i].getInstrumentSettings().getScanMode()>0)
          {
            item->setText(QString::fromStdString((*cl.getPeakData())[i].getInstrumentSettings().NamesOfScanMode[(*cl.getPeakData())[i].getInstrumentSettings().getScanMode()]));
          }
          else
          {
            item->setText("-");
          }
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 5, item);

          // zoom scan
          item = new QTableWidgetItem();
          if ((*cl.getPeakData())[i].getInstrumentSettings().getZoomScan())
          {
            item->setText("yes");
          }
          else
          {
            item->setText("no");
          }
          item->setBackgroundColor(c);
          table_widget_->setItem(table_widget_->rowCount()-1 , 6, item);

          if (i == cl.current_spectrum)
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
    }
    else  // no peak map
    {
      return; // leave signals blocked
    }

    if (cl.getPeakData()->size() == 1)
    {
      item->setFlags(0);
      return; // leave signals blocked
    }

    table_widget_->blockSignals(false);
    table_widget_->setSortingEnabled(true);
    table_widget_->resizeColumnsToContents ();
    table_widget_->resizeRowsToContents();
  }

  SpectraIdentificationViewWidget::~SpectraIdentificationViewWidget()
  {
  }

}
