// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/VISUAL/SpectraViewWidget.h>
#include <QtGui/QVBoxLayout>
#include <QtGui/QTreeWidget>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtGui/QHeaderView>
#include <QtGui/QMenu>

namespace OpenMS
{
  SpectraViewWidget::SpectraViewWidget(QWidget * parent) :
    QWidget(parent)
  {
    QVBoxLayout * spectra_widget_layout = new QVBoxLayout(this);
    spectra_treewidget_ = new QTreeWidget(this);
    spectra_treewidget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to open it.");

    //~ no good for huge experiments - omitted:
    //~ spectrum_selection_->setSortingEnabled(true);
    //~ spectrum_selection_->sortByColumn ( 1, Qt::AscendingOrder);

    spectra_treewidget_->setColumnCount(7); /// @improvement make dependend from global "header_labels" to change only once (otherwise changes must be applied in several slots too!)

    spectra_treewidget_->setColumnWidth(0, 65);
    spectra_treewidget_->setColumnWidth(1, 45);
    spectra_treewidget_->setColumnWidth(2, 50);
    spectra_treewidget_->setColumnWidth(3, 55);
    spectra_treewidget_->setColumnWidth(4, 55);
    spectra_treewidget_->setColumnWidth(5, 45);
    spectra_treewidget_->setColumnWidth(6, 45);

    ///@improvement write the visibility-status of the columns in toppview.ini and read at start

    QStringList qsl; // names of searchable columns
    qsl << "index" << "RT" << "PC m/z" << "dissociation" << "scan" << "zoom";

    QStringList header_labels; /// @improvement make this global to change only once (otherwise changes must be applied in several slots too!)
    header_labels.append(QString("MS level"));
    header_labels.append(qsl); // all searchable columns
    spectra_treewidget_->setHeaderLabels(header_labels);

    spectra_treewidget_->setDragEnabled(true);
    spectra_treewidget_->setContextMenuPolicy(Qt::CustomContextMenu);
    spectra_treewidget_->header()->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(spectra_treewidget_, SIGNAL(currentItemChanged(QTreeWidgetItem *, QTreeWidgetItem *)), this, SLOT(spectrumSelectionChange_(QTreeWidgetItem *, QTreeWidgetItem *)));
    connect(spectra_treewidget_, SIGNAL(itemDoubleClicked(QTreeWidgetItem *, int)), this, SLOT(spectrumDoubleClicked_(QTreeWidgetItem *)));
    connect(spectra_treewidget_, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(spectrumContextMenu_(const QPoint &)));
    connect(spectra_treewidget_->header(), SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(spectrumBrowserHeaderContextMenu_(const QPoint &)));

    spectra_widget_layout->addWidget(spectra_treewidget_);

    QHBoxLayout * tmp_hbox_layout = new QHBoxLayout();

    spectra_search_box_ = new QLineEdit("<search text>", this);
    spectra_search_box_->setWhatsThis("Search in a certain column. Hits are shown as you type. Press <Enter> to display the first hit.");

    spectra_combo_box_ = new QComboBox(this);
    spectra_combo_box_->addItems(qsl);
    spectra_combo_box_->setWhatsThis("Sets the column in which to search.");

    // search whenever text is typed (and highlight the hits)
    connect(spectra_search_box_, SIGNAL(textEdited(const QString &)), this, SLOT(spectrumSearchText_()));
    // .. show hit upon pressing Enter (internally we search again, since the user could have activated another layer with different selections after last search)
    connect(spectra_search_box_, SIGNAL(returnPressed()), this, SLOT(searchAndShow_()));

    tmp_hbox_layout->addWidget(spectra_search_box_);
    tmp_hbox_layout->addWidget(spectra_combo_box_);
    spectra_widget_layout->addLayout(tmp_hbox_layout);
  }

  QTreeWidget * SpectraViewWidget::getTreeWidget()
  {
    return spectra_treewidget_;
  }

  QComboBox * SpectraViewWidget::getComboBox()
  {
    return spectra_combo_box_;
  }

  void SpectraViewWidget::spectrumSearchText_()
  {
    const QString text = spectra_search_box_->text(); // get text from QLineEdit
    if (text.size() > 0)
    {
      QTreeWidget * spectra_view_treewidget = spectra_treewidget_;
      QComboBox * spectra_view_combobox = spectra_combo_box_;
      int col(spectra_view_combobox->currentIndex() + 1);
      if (col > 5)
      {
        col = 1;
      }

      Qt::MatchFlags matchflags = Qt::MatchFixedString;
      matchflags |=  Qt::MatchRecursive; // match subitems (below top-level)
      if (col != 1)
      {
        // only the index has to be matched exactly
        matchflags = matchflags | Qt::MatchStartsWith;
      }
      QList<QTreeWidgetItem *> searched = spectra_view_treewidget->findItems(text, matchflags, col);

      if (searched.size() > 0)
      {
        //QTreeWidgetItem * olditem = spectra_view_treewidget->currentItem();
        spectra_view_treewidget->clearSelection();
        searched.first()->setSelected(true);
        spectra_view_treewidget->update();
        spectra_view_treewidget->scrollToItem(searched.first());
        //spectrumSelectionChange_(searched.first(), olditem); // updates the plot
      }
    }
  }

  void SpectraViewWidget::spectrumSelectionChange_(QTreeWidgetItem * current, QTreeWidgetItem * previous)
  {
    /*	test for previous == 0 is important - without it,
        the wrong spectrum will be selected after finishing
        the execution of a TOPP tool on the whole data */
    if (current == 0 || previous == 0)
    {
      return;
    }

    int spectrum_index = current->text(1).toInt();

    const QList<QVariant> & res = current->data(0, 0).toList();
    if (res.size() == 0)
    {
      emit spectrumSelected(spectrum_index);
    }
    else
    {
      // open several chromatograms at once
      std::vector<int> chrom_indices;
      for (Int i = 0; i != res.size(); ++i)
      {
        chrom_indices.push_back(res[i].toInt());
      }
      emit spectrumSelected(chrom_indices);
    }

  }

  void SpectraViewWidget::searchAndShow_()
  {
    //QTreeWidgetItem* current = spectra_treewidget_->currentItem();
    spectrumSearchText_(); // update selection first (we might be in a new layer)
    QList<QTreeWidgetItem *> selected = spectra_treewidget_->selectedItems();
    if (selected.size() > 0) spectrumSelectionChange_(selected.first(), selected.first());
  }

  void SpectraViewWidget::spectrumDoubleClicked_(QTreeWidgetItem * current)
  {
    if (current == 0)
    {
      return;
    }
    int spectrum_index = current->text(1).toInt();
    const QList<QVariant> & res = current->data(0, 0).toList();
    if (res.size() == 0)
    {
      emit spectrumDoubleClicked(spectrum_index);
    }
    else
    {
      // open several chromatograms at once
      std::vector<int> chrom_indices;
      for (Int i = 0; i != res.size(); ++i)
      {
        chrom_indices.push_back(res[i].toInt());
      }
      emit spectrumDoubleClicked(chrom_indices);
    }

  }

  void SpectraViewWidget::spectrumContextMenu_(const QPoint & pos)
  {
    QTreeWidgetItem * item = spectra_treewidget_->itemAt(pos);
    if (item)
    {
      //create menu
      int spectrum_index = item->text(1).toInt();
      QMenu * context_menu = new QMenu(spectra_treewidget_);
      context_menu->addAction("Show in 1D view");
      context_menu->addAction("Meta data");
      context_menu->addAction("Center here");

      QAction * selected = context_menu->exec(spectra_treewidget_->mapToGlobal(pos));
      if (selected != 0 && selected->text() == "Show in 1D view")
      {
        std::vector<int> chrom_indices;
        const QList<QVariant> & res = item->data(0, 0).toList();
        if (res.size() == 0)
        {
          emit showSpectrumAs1D(spectrum_index);
        }
        else
        {
          // open several chromatograms at once
          for (Int i = 0; i != res.size(); ++i)
          {
            chrom_indices.push_back(res[i].toInt());
          }
          emit showSpectrumAs1D(chrom_indices);
        }
      }
      else if (selected != 0 && selected->text() == "Meta data")
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

  void SpectraViewWidget::spectrumBrowserHeaderContextMenu_(const QPoint & pos)
  {
    //create menu
    QMenu * context_menu = new QMenu(spectra_treewidget_->header());

    QStringList header_labels;
    header_labels.append(QString("MS level"));
    header_labels.append(QString("index"));
    header_labels.append(QString("RT"));
    header_labels.append(QString("precursor m/z"));
    header_labels.append(QString("dissociation"));
    header_labels.append(QString("scan type"));
    header_labels.append(QString("zoom"));
    for (int i = 0; i < header_labels.size(); ++i)
    {
      QAction * tmp = new QAction(header_labels[i], context_menu);
      tmp->setCheckable(true);
      tmp->setChecked(!spectra_treewidget_->isColumnHidden(i));
      context_menu->addAction(tmp);
    }

    //(show and) execute menu
    QAction * selected = context_menu->exec(spectra_treewidget_->mapToGlobal(pos));
    if (selected != 0)
    {
      for (int i = 0; i < header_labels.size(); ++i)
      {
        if (selected->text() == header_labels[i])
        {
          selected->isChecked() ? spectra_treewidget_->setColumnHidden(i, false)
          : spectra_treewidget_->setColumnHidden(i, true);
        }
      }
    }
    delete (context_menu);
  }

  void SpectraViewWidget::updateEntries(const LayerData & cl)
  {
    if (!spectra_treewidget_->isVisible() || spectra_treewidget_->signalsBlocked())
    {
      return;
    }

    spectra_treewidget_->blockSignals(true);
    spectra_treewidget_->clear();

    QTreeWidgetItem * item = 0;
    QTreeWidgetItem * selected_item = 0;
    QList<QTreeWidgetItem *> toplevel_items;
    bool more_than_one_spectrum = true;

    // Branch if the current layer is a spectrum
    if (cl.type == LayerData::DT_PEAK  && !(cl.chromatogram_flag_set()))
    {
      std::vector<QTreeWidgetItem *> parent_stack;
      parent_stack.push_back(0);
      bool fail = false;

      for (Size i = 0; i < cl.getPeakData()->size(); ++i)
      {
        const MSSpectrum& current_spec = (*cl.getPeakData())[i];

        if (i > 0)
        {
          const MSSpectrum& prev_spec = (*cl.getPeakData())[i-1];
          // current MS level = previous MS level + 1 (e.g. current: MS2, previous: MS1)
          if (current_spec.getMSLevel() == prev_spec.getMSLevel() + 1)
          {
            item = new QTreeWidgetItem(parent_stack.back());
            parent_stack.resize(parent_stack.size() + 1);
          }
          // current MS level = previous MS level (e.g. MS2,MS2 or MS1,MS1)
          else if (current_spec.getMSLevel() == prev_spec.getMSLevel())
          {
            if (parent_stack.size() == 1)
            {
              item = new QTreeWidgetItem((QTreeWidget *)0);
            }
            else
            {
              item = new QTreeWidgetItem(*(parent_stack.end() - 2));
            }
          }
          // current MS level < previous MS level (e.g. MS1,MS2)
          else if (current_spec.getMSLevel() < prev_spec.getMSLevel())
          {
            Int level_diff = prev_spec.getMSLevel() - current_spec.getMSLevel();
            Size parent_index = 0;
            if (parent_stack.size() - level_diff >= 2)
            {
              QTreeWidgetItem * parent = 0;
              parent_index = parent_stack.size() - level_diff - 1;
              parent = parent_stack[parent_index];

              item = new QTreeWidgetItem(parent, parent_stack[parent_index + 1]);
            }
            else
            {
              item = new QTreeWidgetItem((QTreeWidget *)0);
            }
            parent_stack.resize(parent_index + 1);
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
          item = new QTreeWidgetItem((QTreeWidget *)0);
        }

        parent_stack.back() = item;
        if (parent_stack.size() == 1)
        {
          toplevel_items.push_back(item);
        }

        item->setText(0, QString("MS") + QString::number(current_spec.getMSLevel()));
        item->setText(1, QString::number(i));
        item->setText(2, QString::number(current_spec.getRT()));

        const std::vector<Precursor>& current_precursors = current_spec.getPrecursors();

        if (!current_precursors.empty() || current_spec.metaValueExists("analyzer scan offset"))
        {
          double pc_val;
          if (current_spec.metaValueExists("analyzer scan offset"))
          {
            pc_val = current_spec.getMetaValue("analyzer scan offset");
            item->setText(4, "-");
          }
          else 
          {
            const Precursor& current_pc = current_precursors[0];
            pc_val = current_pc.getMZ();
            if (!current_pc.getActivationMethods().empty())
            {
              QString t;
              for (std::set<Precursor::ActivationMethod>::const_iterator it = current_pc.getActivationMethods().begin(); it != current_pc.getActivationMethods().end(); ++it)
              {
                if (!t.isEmpty())
                {
                  t.append(",");
                }
                t.append(QString::fromStdString(current_pc.NamesOfActivationMethod[*(current_pc.getActivationMethods().begin())]));
              }
              item->setText(4, t);
            }
            else
            {
              item->setText(4, "-");
            }
          }
          item->setText(3, QString::number(pc_val));
        }
        else
        {
          item->setText(3, "-");
          item->setText(4, "-");
        }
        if (current_spec.getInstrumentSettings().getScanMode() > 0)
        {
          item->setText(5, QString::fromStdString(current_spec.getInstrumentSettings().NamesOfScanMode[current_spec.getInstrumentSettings().getScanMode()]));
        }
        else
        {
          item->setText(5, "-");
        }
        if (current_spec.getInstrumentSettings().getZoomScan())
        {
          item->setText(6, "yes");
        }
        else
        {
          item->setText(6, "no");
        }
        /*
        std::cout << "adding: ";
	      for (Size k = 0; k != item->columnCount(); ++k)
	      {
          std::cout << item->text(k).toStdString() << " ";
	      }
        std::cout << std::endl;
        */
        if (i == cl.getCurrentSpectrumIndex())
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
          const MSSpectrum& current_spec = (*cl.getPeakData())[i];
          item = new QTreeWidgetItem((QTreeWidget *)0);
          item->setText(0, QString("MS") + QString::number(current_spec.getMSLevel()));
          item->setText(1, QString::number(i));
          item->setText(2, QString::number(current_spec.getRT()));
          if (!current_spec.getPrecursors().empty())
          {
            item->setText(3, QString::number(current_spec.getPrecursors()[0].getMZ()));
            if (!current_spec.getPrecursors().front().getActivationMethods().empty())
            {
              QString t;
              for (std::set<Precursor::ActivationMethod>::const_iterator it = current_spec.getPrecursors().front().getActivationMethods().begin(); it != current_spec.getPrecursors().front().getActivationMethods().end(); ++it)
              {
                if (!t.isEmpty())
                {
                  t.append(",");
                }
                t.append(QString::fromStdString(current_spec.getPrecursors().front().NamesOfActivationMethod[*(current_spec.getPrecursors().front().getActivationMethods().begin())]));
              }
              item->setText(4, t);
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
          if (current_spec.getInstrumentSettings().getScanMode() > 0)
          {
            item->setText(5, QString::fromStdString(current_spec.getInstrumentSettings().NamesOfScanMode[current_spec.getInstrumentSettings().getScanMode()]));
          }
          else
          {
            item->setText(5, "-");
          }
          if (current_spec.getInstrumentSettings().getZoomScan())
          {
            item->setText(6, "yes");
          }
          else
          {
            item->setText(6, "no");
          }
          toplevel_items.push_back(item);
          if (i == cl.getCurrentSpectrumIndex())
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
      if (cl.getPeakData()->size() > 1)
      {
        more_than_one_spectrum = false;
      }
    }
    // Branch if the current layer is a chromatogram (either indicated by its
    // type or by the flag which is set).
    else if (cl.type == LayerData::DT_CHROMATOGRAM || cl.chromatogram_flag_set())
    {

      // We need to redraw the whole Widget because the we have changed all the layers.
      // First we need to figure out which chromatogram was selected and
      // whether multiple ones are selected.

      bool multiple_select = false;
      int this_selected_item = -1;
      if (cl.getPeakData()->size() > 0 && cl.getPeakData()->metaValueExists("multiple_select"))
      {
        multiple_select = cl.getPeakData()->getMetaValue("multiple_select").toBool();
      }
      if (cl.getPeakData()->size() > 0 && cl.getPeakData()->metaValueExists("selected_chromatogram"))
      {
        this_selected_item = (int)cl.getPeakData()->getMetaValue("selected_chromatogram");
      }

      spectra_treewidget_->setColumnCount(5);
      spectra_treewidget_->setColumnWidth(0, 45);
      spectra_treewidget_->setColumnWidth(1, 45);
      spectra_treewidget_->setColumnWidth(2, 80);
      spectra_treewidget_->setColumnWidth(3, 150);
      spectra_treewidget_->setColumnWidth(4, 80);
      spectra_treewidget_->setColumnWidth(5, 80);
      spectra_treewidget_->setColumnWidth(6, 45);
      spectra_treewidget_->setColumnWidth(7, 80);

      // create a different header list
      QStringList header_labels;
      header_labels.append(QString(" type "));
      header_labels.append(QString("index"));
      header_labels.append(QString("m/z"));
      header_labels.append(QString("Description"));
      header_labels.append(QString("rt start"));
      header_labels.append(QString("rt end"));
      header_labels.append(QString("charge"));
      header_labels.append(QString("chromatogram type"));
      spectra_treewidget_->setHeaderLabels(header_labels);

      // create a different combo box
      int curr = spectra_combo_box_->currentIndex();
      QStringList qsl;
      qsl.push_back("index");
      qsl.push_back("m/z");
      qsl.push_back("Description");
      spectra_combo_box_->clear();
      spectra_combo_box_->addItems(qsl);
      spectra_combo_box_->setCurrentIndex(curr);

      LayerData::ExperimentSharedPtrType exp;
      exp = cl.getPeakData();

      if (cl.chromatogram_flag_set())
      {
        exp = cl.getChromatogramData();
      }

      if (exp->getChromatograms().size() > 1)
      {
        more_than_one_spectrum = false;
      }

      // try to retrieve the map from the cache if available
      typedef std::set<Precursor, Precursor::MZLess> PCSetType;
      std::map<Precursor, std::vector<Size>, Precursor::MZLess> map_precursor_to_chrom_idx;
      if (map_precursor_to_chrom_idx_cache_.find((size_t)(exp.get())) != map_precursor_to_chrom_idx_cache_.end())
      {
        map_precursor_to_chrom_idx = map_precursor_to_chrom_idx_cache_[(size_t)(exp.get())];
      }
      else
      {

        // collect all precursor that fall into the mz rt window
        PCSetType precursor_in_rt_mz_window;
        for (std::vector<MSChromatogram >::const_iterator iter = exp->getChromatograms().begin(); iter != exp->getChromatograms().end(); ++iter)
        {
          precursor_in_rt_mz_window.insert(iter->getPrecursor());
        }

        // determine product chromatograms for each precursor
        for (PCSetType::const_iterator pit = precursor_in_rt_mz_window.begin(); pit != precursor_in_rt_mz_window.end(); ++pit)
        {
          for (std::vector<MSChromatogram >::const_iterator iter = exp->getChromatograms().begin(); iter != exp->getChromatograms().end(); ++iter)
          {
            if (iter->getPrecursor() == *pit)
            {
              map_precursor_to_chrom_idx[*pit].push_back(iter - exp->getChromatograms().begin());
            }
          }
        }

        map_precursor_to_chrom_idx_cache_[(size_t)(exp.get())] = map_precursor_to_chrom_idx;
      }

      if (!map_precursor_to_chrom_idx.empty())
      {
        int precursor_idx = 0;
        for (std::map<Precursor, std::vector<Size>, Precursor::MZLess>::iterator mit = map_precursor_to_chrom_idx.begin(); mit != map_precursor_to_chrom_idx.end(); ++mit)
        {
          // Show the peptide sequence if available, otherwise show the m/z and charge only
          QString mz_string = QString::number(mit->first.getMZ());
          QString charge = QString::number(mit->first.getCharge());
          QString description = "";
          if (mit->first.metaValueExists("peptide_sequence"))
          {
            description = String(mit->first.getMetaValue("peptide_sequence")).toQString();
          }
          if (mit->first.metaValueExists("description"))
          {
            description = String(mit->first.getMetaValue("description")).toQString();
          }

          // Show all: iterate over all chromatograms corresponding to the current precursor and add action containing all chromatograms
          QList<QVariant> chroms_idx;
          for (std::vector<Size>::iterator vit = mit->second.begin(); vit != mit->second.end(); ++vit)
          {
            chroms_idx.push_back((unsigned int)*vit);
          }

          bool one_selected = false;

          // Top level precursor entry
          item = new QTreeWidgetItem(0);
          item->setText(0, QString("Peptide"));
          item->setText(1, QString::number(precursor_idx++));
          item->setText(2, mz_string);
          item->setText(3, description);
          //item->setText(4, QString::number(prod_it->second[0].front().getRT()));
          //item->setText(5, QString::number(prod_it->second[0].back().getRT()));
          item->setText(6, QString("-"));
          item->setText(7, charge);
          item->setData(0, 0, chroms_idx);

          toplevel_items.push_back(item);

          // Show single chromatogram: iterate over all chromatograms corresponding to the current precursor and add action for the single chromatogram
          for (std::vector<Size>::iterator vit = mit->second.begin(); vit != mit->second.end(); ++vit)
          {
            const MSChromatogram & current_chromatogram = exp->getChromatograms()[*vit];

            // Children chromatogram entry
            QTreeWidgetItem * sub_item = new QTreeWidgetItem(item);
            if ((int)*vit == this_selected_item)
            {
              one_selected = true;
              selected_item = sub_item;
            }
            QString chrom_description = "ion";
            if (mit->first.metaValueExists("description"))
            {
              chrom_description = String(mit->first.getMetaValue("description")).toQString();
            }

            sub_item->setText(0, QString("Transition"));
            sub_item->setText(1, QString::number((unsigned int)*vit));
            sub_item->setText(2, QString::number(current_chromatogram.getProduct().getMZ()));
            //sub_item->setText(7, QString::number(prod_it->second[0].getProduct().getCharge())); // TODO product charge
            sub_item->setText(3, QString(chrom_description));
            if (! current_chromatogram.empty())
            {
              sub_item->setText(4, QString::number(current_chromatogram.front().getRT()));
              sub_item->setText(5, QString::number(current_chromatogram.back().getRT()));
            }

            switch (current_chromatogram.getChromatogramType())
            {
            case ChromatogramSettings::MASS_CHROMATOGRAM:                         sub_item->setText(6, QString("Mass chromatogram")); break;

            case ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM:            sub_item->setText(6, QString("Total ion chromatogram")); break;

            case ChromatogramSettings::SELECTED_ION_CURRENT_CHROMATOGRAM:         sub_item->setText(6, QString("Selected ion current chromatogram")); break;

            case ChromatogramSettings::BASEPEAK_CHROMATOGRAM:                     sub_item->setText(6, QString("Basepeak chromatogram")); break;

            case ChromatogramSettings::SELECTED_ION_MONITORING_CHROMATOGRAM:      sub_item->setText(6, QString("Selected ion monitoring chromatogram")); break;

            case ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM: sub_item->setText(6, QString("Selected reaction monitoring chromatogram")); break;

            case ChromatogramSettings::ELECTROMAGNETIC_RADIATION_CHROMATOGRAM:    sub_item->setText(6, QString("Electromagnetic radiation chromatogram")); break;

            case ChromatogramSettings::ABSORPTION_CHROMATOGRAM:                   sub_item->setText(6, QString("Absorption chromatogram")); break;

            case ChromatogramSettings::EMISSION_CHROMATOGRAM:                     sub_item->setText(6, QString("Emission chromatogram")); break;

            default:                                                                            sub_item->setText(6, QString("Unknown chromatogram")); break;
            }

          }
          if (one_selected && multiple_select)
          {
            selected_item = item;
          }
        }
      }
      spectra_treewidget_->addTopLevelItems(toplevel_items);

      if (selected_item && this_selected_item != -1)
      {
        // now, select and scroll down to item
        spectra_treewidget_->setCurrentItem(selected_item);
        selected_item->setSelected(true);
        spectra_treewidget_->scrollToItem(selected_item);

        // expand the item if necessary
        if (!multiple_select)
        {
          selected_item->parent()->setExpanded(true);
        }
      }
    }
    // Branch if its neither (just draw an empty item)
    else
    {
      item = new QTreeWidgetItem((QTreeWidget *)0);
      item->setText(0, QString("No peak map"));
      item->setText(1, QString("-"));
      item->setText(2, QString("-"));
      item->setText(3, QString::number(0));
      item->setFlags(0);
      spectra_treewidget_->addTopLevelItem(item);
    }

    if (more_than_one_spectrum && item != 0)
    {
      item->setFlags(0);
    }

    spectra_treewidget_->blockSignals(false);
  }

  SpectraViewWidget::~SpectraViewWidget()
  {
  }

}
