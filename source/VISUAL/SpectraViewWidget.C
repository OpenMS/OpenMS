// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/VISUAL/SpectraViewWidget.h>
#include <QtGui/QVBoxLayout>
#include <QtGui/QTreeWidget>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtGui/QHeaderView>
#include <QtGui/QMenu>

namespace OpenMS
{
  SpectraViewWidget::SpectraViewWidget(QWidget* parent)
    : QWidget(parent)
  {
    QVBoxLayout* spectra_widget_layout = new QVBoxLayout(this);
    spectra_treewidget_ = new QTreeWidget(this);
    spectra_treewidget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to open it.");

    //~ no good for huge experiments - omitted:
    //~ spectrum_selection_->setSortingEnabled(true);
    //~ spectrum_selection_->sortByColumn ( 1, Qt::AscendingOrder);

    spectra_treewidget_->setColumnCount(7);/// @improvement make dependend from global "header_labels" to change only once (otherwise changes must be applied in several slots too!)

    spectra_treewidget_->setColumnWidth(0,65);
    spectra_treewidget_->setColumnWidth(1,45);
    spectra_treewidget_->setColumnWidth(2,50);
    spectra_treewidget_->setColumnWidth(3,55);
    spectra_treewidget_->setColumnWidth(4,55);
    spectra_treewidget_->setColumnWidth(5,45);
    spectra_treewidget_->setColumnWidth(6,45);

    ///@improvement write the visibility-status of the columns in toppview.ini and read at start

    QStringList header_labels; /// @improvement make this global to change only once (otherwise changes must be applied in several slots too!)
    header_labels.append(QString("MS level"));
    header_labels.append(QString("index"));
    header_labels.append(QString("RT"));
    header_labels.append(QString("precursor m/z"));
    header_labels.append(QString("dissociation"));
    header_labels.append(QString("scan type"));
    header_labels.append(QString("zoom"));
    spectra_treewidget_->setHeaderLabels(header_labels);

    spectra_treewidget_->setDragEnabled(true);
    spectra_treewidget_->setContextMenuPolicy(Qt::CustomContextMenu);
    spectra_treewidget_->header()->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(spectra_treewidget_, SIGNAL(currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*)), this, SLOT(spectrumSelectionChange_(QTreeWidgetItem*, QTreeWidgetItem*)));
    connect(spectra_treewidget_, SIGNAL(itemDoubleClicked(QTreeWidgetItem*, int)), this, SLOT(spectrumDoubleClicked_(QTreeWidgetItem*, int)));
    connect(spectra_treewidget_, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(spectrumContextMenu_(const QPoint&)));
    connect(spectra_treewidget_->header(), SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(spectrumBrowserHeaderContextMenu_(const QPoint&)));

    spectra_widget_layout->addWidget(spectra_treewidget_);

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
  }

  QTreeWidget* SpectraViewWidget::getTreeWidget()
  {
    return spectra_treewidget_;
  }

  QComboBox* SpectraViewWidget::getComboBox()
  {
    return spectra_combo_box_;
  }

  void SpectraViewWidget::spectrumSelected_(const QString& text)
  {
    QTreeWidget* spectra_view_treewidget = spectra_treewidget_;
    QComboBox* spectra_view_combobox = spectra_combo_box_;
    if(text.size() > 0)
    {
      int col(spectra_view_combobox->currentIndex()+1);
      if(col>5)
      {
        col = 1;
      }
      QList<QTreeWidgetItem *> searched = spectra_view_treewidget->findItems(text, Qt::MatchFixedString /* matchflag exact match */, col );
      QList<QTreeWidgetItem *> selected = spectra_view_treewidget->selectedItems();

      if(searched.size() > 0)
      {
        for(int i = 0; i < selected.size(); ++i)
        {
          selected[i]->setSelected(false);
        }
        spectra_view_treewidget->update();
        int spectrum_index = searched.first()->text(1).toInt();
        searched.first()->setSelected(true);
        spectra_view_treewidget->update();
        spectra_view_treewidget->scrollToItem(searched.first());

        emit spectrumSelected(spectrum_index);
      }
    }
  }

  void SpectraViewWidget::spectrumSelectionChange_(QTreeWidgetItem* current, QTreeWidgetItem* previous)
  {
    /*	test for previous == 0 is important - without it,
        the wrong spectrum will be selected after finishing
        the execution of a TOPP tool on the whole data */
    if (current == 0 || previous == 0)
    {
      return;
    }

    int spectrum_index = current->text(1).toInt();
    emit spectrumSelected(spectrum_index);
  }

  void SpectraViewWidget::spectrumDoubleClicked_(QTreeWidgetItem* current, int)
  {
    if (current == 0)
    {
      return;
    }
    int spectrum_index = current->text(1).toInt();
    emit spectrumDoubleClicked(spectrum_index);
  }

  void SpectraViewWidget::spectrumContextMenu_(const QPoint& pos)
  {
    QTreeWidgetItem* item = spectra_treewidget_->itemAt(pos);
    if (item)
    {
      //create menu
      int spectrum_index = item->text(1).toInt();
      QMenu* context_menu = new QMenu(spectra_treewidget_);
      context_menu->addAction("Show in 1D view");
      context_menu->addAction("Meta data");
      context_menu->addAction("Center here");

      QAction* selected = context_menu->exec(spectra_treewidget_->mapToGlobal(pos));
      if (selected!=0 && selected->text()=="Show in 1D view")
      {
        emit showSpectrumAs1D(spectrum_index);
      }
      else if (selected!=0 && selected->text()=="Meta data")
      {
        emit showSpectrumMetaData(spectrum_index);
      }
      /** TODO
      else if (selected!=0 && selected->text()=="Center here")
      {
        emit centerHere(spectrum_index);
      }
      **/
      delete (context_menu);
    }
  }

  void SpectraViewWidget::spectrumBrowserHeaderContextMenu_(const QPoint& pos)
  {
    //create menu
    QMenu* context_menu = new QMenu(spectra_treewidget_->header());

    QStringList header_labels;
    header_labels.append(QString("MS level"));
    header_labels.append(QString("index"));
    header_labels.append(QString("RT"));
    header_labels.append(QString("precursor m/z"));
    header_labels.append(QString("dissociation"));
    header_labels.append(QString("scan type"));
    header_labels.append(QString("zoom"));
    for(int i = 0; i < header_labels.size(); ++i)
    {
      QAction* tmp = new QAction(header_labels[i],context_menu);
      tmp->setCheckable(true);
      tmp->setChecked(!spectra_treewidget_->isColumnHidden(i));
      context_menu->addAction(tmp);
    }

    //(show and) execute menu
    QAction* selected = context_menu->exec(spectra_treewidget_->mapToGlobal(pos));
    if (selected!=0)
    {
      for(int i = 0; i < header_labels.size(); ++i)
      {
        if(selected->text() == header_labels[i])
        {
          selected->isChecked()? spectra_treewidget_->setColumnHidden(i,false)
            :spectra_treewidget_->setColumnHidden(i,true);
        }
      }
    }
    delete (context_menu);
  }

  void SpectraViewWidget::updateEntries(const LayerData& cl)
  {
    if (!spectra_treewidget_->isVisible())
    {
      return;
    }

    spectra_treewidget_->blockSignals(true);
    spectra_treewidget_->clear();

    QTreeWidgetItem* item = 0;
    QTreeWidgetItem* selected_item = 0;
    QList<QTreeWidgetItem*> toplevel_items;

    if(cl.type == LayerData::DT_PEAK)
    {
      std::vector<QTreeWidgetItem*> parent_stack;
      parent_stack.push_back(0);
      bool fail = false;

      for (Size i = 0; i < cl.getPeakData()->size(); ++i)
      {
        if (i > 0)
        {
          // current MS level = previous MS level + 1 (e.g. current: MS2, previous: MS1)
          if ((*cl.getPeakData())[i].getMSLevel() == (*cl.getPeakData())[i-1].getMSLevel() + 1)
          {
            item = new QTreeWidgetItem(parent_stack.back());
            parent_stack.resize(parent_stack.size()+1);
          }
          // current MS level = previous MS level (e.g. MS2,MS2 or MS1,MS1)
          else if ((*cl.getPeakData())[i].getMSLevel() == (*cl.getPeakData())[i-1].getMSLevel())
          {
            if (parent_stack.size() == 1)
            {
              item = new QTreeWidgetItem((QTreeWidget*)0);
            }
            else
            {
              item = new QTreeWidgetItem(*(parent_stack.end()-2));
            }
          }
          // current MS level < previous MS level (e.g. MS1,MS2)
          else if ((*cl.getPeakData())[i].getMSLevel() < (*cl.getPeakData())[i-1].getMSLevel())
          {
            Int level_diff = (*cl.getPeakData())[i-1].getMSLevel() - (*cl.getPeakData())[i].getMSLevel();
            Size parent_index = 0;
            QTreeWidgetItem* parent = 0;
            if (parent_stack.size() - level_diff >= 2)
            {
              parent_index = parent_stack.size() - level_diff - 1;
              parent = parent_stack[parent_index];

              item = new QTreeWidgetItem(parent, parent_stack[parent_index+1]);
            }
            else
            {
              item = new QTreeWidgetItem((QTreeWidget*)0);
            }
            parent_stack.resize(parent_index+1);
          }
          else
          {
            std::cerr << "Cannot build treelike view for spectrum browser, generating flat list instead." << std::endl;
            fail = true;
            break;
          }
        }
        else
        {
          item = new QTreeWidgetItem((QTreeWidget*)0);
        }

        parent_stack.back() = item;
        if (parent_stack.size() == 1)
        {
          toplevel_items.push_back(item);
        }

        item->setText(0, QString("MS") + QString::number((*cl.getPeakData())[i].getMSLevel()));
        item->setText(1, QString::number(i));
        item->setText(2, QString::number((*cl.getPeakData())[i].getRT()));
        if (!(*cl.getPeakData())[i].getPrecursors().empty())
        {
          item->setText(3,QString::number((*cl.getPeakData())[i].getPrecursors()[0].getMZ()));
          if (!(*cl.getPeakData())[i].getPrecursors().front().getActivationMethods().empty())
          {
            QString t;
            for(std::set<Precursor::ActivationMethod>::const_iterator it = (*cl.getPeakData())[i].getPrecursors().front().getActivationMethods().begin(); it != (*cl.getPeakData())[i].getPrecursors().front().getActivationMethods().end(); ++it)
            {
              if(!t.isEmpty())
              {
                t.append(",");
              }
              t.append(QString::fromStdString((*cl.getPeakData())[i].getPrecursors().front().NamesOfActivationMethod[*((*cl.getPeakData())[i].getPrecursors().front().getActivationMethods().begin())]));
            }
            item->setText(4,t);
          }
          else
          {
            item->setText(4, "-");
          }
        }
        else
        {
          item->setText(3, "-");
          item->setText(4, "-");
        }
        if ((*cl.getPeakData())[i].getInstrumentSettings().getScanMode()>0)
        {
          item->setText(5,QString::fromStdString((*cl.getPeakData())[i].getInstrumentSettings().NamesOfScanMode[(*cl.getPeakData())[i].getInstrumentSettings().getScanMode()]));
        }
        else
        {
          item->setText(5, "-");
        }
        if ((*cl.getPeakData())[i].getInstrumentSettings().getZoomScan())
        {
          item->setText(6,"yes");
        }
        else
        {
          item->setText(6, "no");
        }

        if (i == cl.current_spectrum)
        {
          // just remember it, select later
          selected_item = item;
        }
      }

      if (!fail)
      {
        spectra_treewidget_->addTopLevelItems(toplevel_items);
      }
      else
      {
        // generate flat list instead
        spectra_treewidget_->clear();
        toplevel_items.clear();
        selected_item = 0;
        for (Size i = 0; i < cl.getPeakData()->size(); ++i)
        {
          item = new QTreeWidgetItem((QTreeWidget*)0);
          item->setText(0, QString("MS") + QString::number((*cl.getPeakData())[i].getMSLevel()));
          item->setText(1, QString::number(i));
          item->setText(2, QString::number((*cl.getPeakData())[i].getRT()));
          if (!(*cl.getPeakData())[i].getPrecursors().empty())
          {
            item->setText(3,QString::number((*cl.getPeakData())[i].getPrecursors()[0].getMZ()));
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
              item->setText(4,t);
            }
            else
            {
              item->setText(4, "-");
            }
          }
          else
          {
            item->setText(3, "-");
            item->setText(4, "-");
          }
          if ((*cl.getPeakData())[i].getInstrumentSettings().getScanMode()>0)
          {
            item->setText(5,QString::fromStdString((*cl.getPeakData())[i].getInstrumentSettings().NamesOfScanMode[(*cl.getPeakData())[i].getInstrumentSettings().getScanMode()]));
          }
          else
          {
            item->setText(5, "-");
          }
          if ((*cl.getPeakData())[i].getInstrumentSettings().getZoomScan())
          {
            item->setText(6,"yes");
          }
          else
          {
            item->setText(6, "no");
          }
          toplevel_items.push_back(item);
          if (i == cl.current_spectrum)
          {
            // just remember it, select later
            selected_item = item;
          }
        }
        spectra_treewidget_->addTopLevelItems(toplevel_items);
      }
      if (selected_item)
      {
        // now, select and scroll down to item
        selected_item->setSelected(true);
        spectra_treewidget_->scrollToItem(selected_item);
      }
    }
    else
    {
      item = new QTreeWidgetItem((QTreeWidget*)0);
      item->setText(0, QString("No peak map"));
      item->setText(1, QString("-"));
      item->setText(2, QString("-"));
      item->setText(3, QString::number(0));
      item->setFlags(0);
      spectra_treewidget_->addTopLevelItem(item);
      return; // leave signals blocked
    }

    if (cl.getPeakData()->size() == 1)
    {
      item->setFlags(0);
      return; // leave signals blocked
    }

    spectra_treewidget_->blockSignals(false);
  }

  SpectraViewWidget::~SpectraViewWidget()
  {
  }

}
