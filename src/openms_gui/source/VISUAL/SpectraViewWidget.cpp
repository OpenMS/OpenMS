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

#include <OpenMS/VISUAL/SpectraViewWidget.h>

#include <OpenMS/CONCEPT/RAIICleanup.h>

#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMenu>


namespace OpenMS
{

  std::vector<int> listToVec(const QList<QVariant>& in)
  {
    std::vector<int> out;
    for (Int i = 0; i != in.size(); ++i)
    {
      out.push_back(in[i].toInt());
    }
    return out;
  }

  SpectraViewWidget::SpectraViewWidget(QWidget * parent) :
    QWidget(parent)
  {
    setObjectName("Scans");
    QVBoxLayout * spectra_widget_layout = new QVBoxLayout(this);
    spectra_treewidget_ = new QTreeWidget(this);
    spectra_treewidget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to show it. "
                                      "Double-clicking might be implemented as well, depending on the data. "
                                      "Context-menus for both the column header and data rows are available by right-clicking.");

    //~ no good for huge experiments - omitted:
    //~ spectrum_selection_->setSortingEnabled(true);
    //~ spectrum_selection_->sortByColumn ( 1, Qt::AscendingOrder);

    spectra_treewidget_->setDragEnabled(true);
    spectra_treewidget_->setContextMenuPolicy(Qt::CustomContextMenu);
    spectra_treewidget_->header()->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(spectra_treewidget_, &QTreeWidget::currentItemChanged, this, &SpectraViewWidget::spectrumSelectionChange_);
    connect(spectra_treewidget_, &QTreeWidget::itemDoubleClicked, this, &SpectraViewWidget::spectrumDoubleClicked_);
    connect(spectra_treewidget_, &QTreeWidget::customContextMenuRequested, this, &SpectraViewWidget::spectrumContextMenu_);
    connect(spectra_treewidget_->header(), &QHeaderView::customContextMenuRequested, this, &SpectraViewWidget::spectrumBrowserHeaderContextMenu_);

    spectra_widget_layout->addWidget(spectra_treewidget_);

    QHBoxLayout* tmp_hbox_layout = new QHBoxLayout();

    spectra_search_box_ = new QLineEdit(this);
    spectra_search_box_->setPlaceholderText("<search text>");
    spectra_search_box_->setWhatsThis("Search in a certain column. Hits are shown as you type. Press <Enter> to display the first hit.");
    spectra_search_box_->setToolTip(spectra_search_box_->whatsThis());

    spectra_combo_box_ = new QComboBox(this);
    spectra_combo_box_->setWhatsThis("Sets the column in which to search.");
    spectra_combo_box_->setToolTip(spectra_combo_box_->whatsThis());


    // search whenever text is typed (and highlight the hits)
    connect(spectra_search_box_, &QLineEdit::textEdited, this, &SpectraViewWidget::spectrumSearchText_);
    // .. show hit upon pressing Enter (internally we search again, since the user could have activated another layer with different selections after last search)
    connect(spectra_search_box_, &QLineEdit::returnPressed, this, &SpectraViewWidget::searchAndShow_);

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

  void SpectraViewWidget::spectrumSearchText_()
  {
    const QString text = spectra_search_box_->text(); // get text from QLineEdit
    if (text.size() > 0)
    {
      Qt::MatchFlags matchflags = Qt::MatchFixedString;
      matchflags |=  Qt::MatchRecursive; // match subitems (below top-level)
      if (spectra_combo_box_->currentText().compare("index", Qt::CaseInsensitive) != 0) // strings not equal
      { // only the index has to be matched exactly
        matchflags = matchflags | Qt::MatchStartsWith;
      }
      QList<QTreeWidgetItem*> searched = spectra_treewidget_->findItems(text, matchflags, spectra_combo_box_->currentIndex());

      if (!searched.isEmpty())
      {
        spectra_treewidget_->clearSelection();
        searched.first()->setSelected(true);
        spectra_treewidget_->update();
        spectra_treewidget_->scrollToItem(searched.first());
      }
    }
  }

  void SpectraViewWidget::spectrumSelectionChange_(QTreeWidgetItem* current, QTreeWidgetItem* previous)
  {
    /*	test for previous == 0 is important - without it,
        the wrong spectrum will be selected after finishing
        the execution of a TOPP tool on the whole data */
    if (current == nullptr || previous == nullptr)
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
    { // open several chromatograms at once
      emit spectrumSelected(listToVec(res));
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
    if (current == nullptr)
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
    { // open several chromatograms at once
      emit spectrumDoubleClicked(listToVec(res));
    }

  }

  void SpectraViewWidget::spectrumContextMenu_(const QPoint& pos)
  {
    QTreeWidgetItem* item = spectra_treewidget_->itemAt(pos);
    if (item)
    {
      //create menu
      int spectrum_index = item->text(1).toInt();
      QMenu context_menu(spectra_treewidget_);
      context_menu.addAction("Show in 1D view", [&]()
      {
        std::vector<int> chrom_indices;
        const QList<QVariant>& res = item->data(0, 0).toList();
        if (res.size() == 0)
        {
          emit showSpectrumAs1D(spectrum_index);
        }
        else
        { // open several chromatograms at once
          emit showSpectrumAs1D(listToVec(res));
        }
      });
      context_menu.addAction("Meta data", [&]() 
      {
        emit showSpectrumMetaData(spectrum_index);
      });
      // todo: context_menu->addAction("Center here", [&]() {emit centerHere(spectrum_index); });

      context_menu.exec(spectra_treewidget_->mapToGlobal(pos));
    }
  }

  void SpectraViewWidget::spectrumBrowserHeaderContextMenu_(const QPoint& pos)
  {
    // allows to hide/show columns
    QMenu context_menu(spectra_treewidget_->header());
    const auto& header = spectra_treewidget_->headerItem();

    for (int i = 0; i < header->columnCount(); ++i)
    {
      auto action = context_menu.addAction(header->text(i), [i, this](){
        spectra_treewidget_->setColumnHidden(i, !spectra_treewidget_->isColumnHidden(i));
      });
      action->setCheckable(true);
      action->setChecked(!spectra_treewidget_->isColumnHidden(i));
    }
    
    // show and execute menu
    context_menu.exec(spectra_treewidget_->mapToGlobal(pos));
  }

  void populateRow_(QTreeWidgetItem* item, const int index, const MSSpectrum& spec)
  {
    item->setText(0, QString("MS") + QString::number(spec.getMSLevel()));
    item->setText(1, QString::number(index));
    item->setText(2, QString::number(spec.getRT()));

    const std::vector<Precursor>& current_precursors = spec.getPrecursors();

    if (!current_precursors.empty() || spec.metaValueExists("analyzer scan offset"))
    {
      double precursor_mz;
      if (spec.metaValueExists("analyzer scan offset"))
      {
        precursor_mz = spec.getMetaValue("analyzer scan offset");
      }
      else
      {
        const Precursor& current_pc = current_precursors[0];
        precursor_mz = current_pc.getMZ();
        item->setText(4, ListUtils::concatenate(current_pc.getActivationMethodsAsString(), ",").toQString());
      }
      item->setText(3, QString::number(precursor_mz));
    }

    item->setText(5, QString::fromStdString(spec.getInstrumentSettings().NamesOfScanMode[spec.getInstrumentSettings().getScanMode()]));
    item->setText(6, (spec.getInstrumentSettings().getZoomScan() ? "yes" : "no"));
  }

  void SpectraViewWidget::updateEntries(const LayerData& cl)
  {
    if (!spectra_treewidget_->isVisible() || spectra_treewidget_->signalsBlocked())
    {
      return;
    }

    spectra_treewidget_->blockSignals(true);
    RAIICleanup clean([&](){ spectra_treewidget_->blockSignals(false); });

    QTreeWidgetItem* item = nullptr;
    QTreeWidgetItem* selected_item = nullptr;
    QList<QTreeWidgetItem*> toplevel_items;
    bool more_than_one_spectrum = true;

    has_data_ = true; // for now ...

    // Branch if the current layer is a spectrum
    if (cl.type == LayerData::DT_PEAK  && !(cl.chromatogram_flag_set()))
    {
      spectra_treewidget_->clear();

      std::vector<QTreeWidgetItem *> parent_stack;
      parent_stack.push_back(nullptr);
      bool fail = false;

      QStringList header_labels;
      header_labels << "MS level" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan" << "zoom";
      spectra_treewidget_->setHeaderLabels(header_labels);
      spectra_treewidget_->setColumnCount(header_labels.size());

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
              item = new QTreeWidgetItem((QTreeWidget *)nullptr);
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
              parent_index = parent_stack.size() - level_diff - 1;
              QTreeWidgetItem * parent = parent_stack[parent_index];
              item = new QTreeWidgetItem(parent, parent_stack[parent_index + 1]);
            }
            else
            {
              item = new QTreeWidgetItem((QTreeWidget *)nullptr);
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
          item = new QTreeWidgetItem((QTreeWidget *)nullptr);
        }

        parent_stack.back() = item;
        if (parent_stack.size() == 1)
        {
          toplevel_items.push_back(item);
        }

        populateRow_(item, i, current_spec);

        if (i == cl.getCurrentSpectrumIndex())
        {
          // just remember it, select later
          selected_item = item;
        }
      }

      if (fail)
      {
        // generate flat list instead
        spectra_treewidget_->clear();
        toplevel_items.clear();
        selected_item = nullptr;
        for (Size i = 0; i < cl.getPeakData()->size(); ++i)
        {
          const MSSpectrum& current_spec = (*cl.getPeakData())[i];
          item = new QTreeWidgetItem((QTreeWidget *)nullptr);
          
          populateRow_(item, i, current_spec);

          toplevel_items.push_back(item);
          if (i == cl.getCurrentSpectrumIndex())
          {
            // just remember it, select later
            selected_item = item;
          }
        }
      }
      spectra_treewidget_->addTopLevelItems(toplevel_items);

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
      LayerData::ConstExperimentSharedPtrType exp = (cl.chromatogram_flag_set()
                                                     ? cl.getChromatogramData()
                                                     : cl.getPeakData());
      
      if (last_peakmap_ == exp.get())
      { // underlying data did not change (which is ALWAYS the chromatograms, never peakdata!)
        // --> Do not update (could be many 10k entries for sqMass data and the lag would be unbearable ...)
        return;
      }
      
      last_peakmap_ = exp.get();
      spectra_treewidget_->clear();
      // New data:
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

      // create a different header list
      QStringList header_labels = QStringList() << " type " << "index" << "m/z" << "Description" << "rt start" << "rt end" << "charge" << "chromatogram type";
      spectra_treewidget_->setHeaderLabels(header_labels);
      spectra_treewidget_->setColumnCount(header_labels.size());
           
      if (exp->getChromatograms().size() > 1)
      {
        more_than_one_spectrum = false;
      }

      // try to retrieve the map from the cache if available
      // TODO: same precursor mass / different precursors are not supported! 
      bool was_cached = map_precursor_to_chrom_idx_cache_.find((size_t)(exp.get())) != map_precursor_to_chrom_idx_cache_.end();
      // create new cache or get the existing one
      std::map<Precursor, std::vector<Size>, Precursor::MZLess>& map_precursor_to_chrom_idx = map_precursor_to_chrom_idx_cache_[(size_t)(exp.get())];
      if (!was_cached)
      { // create cache: collect all precursor that fall into the mz rt window
        for (std::vector<MSChromatogram >::const_iterator iter = exp->getChromatograms().begin(); iter != exp->getChromatograms().end(); ++iter)
        {
          map_precursor_to_chrom_idx[iter->getPrecursor()].push_back(iter - exp->getChromatograms().begin());
        }
      }

      if (!map_precursor_to_chrom_idx.empty())
      {
        int precursor_idx = 0;
        for (auto mit = map_precursor_to_chrom_idx.begin(); mit != map_precursor_to_chrom_idx.end(); ++mit)
        {
          // Show the peptide sequence if available, otherwise show the m/z and charge only
          QString mz_string = QString::number(mit->first.getMZ());
          QString charge = QString::number(mit->first.getCharge());
          QString description = "";
          if (mit->first.metaValueExists("description"))
          {
            description = String(mit->first.getMetaValue("description")).toQString();
          }
          if (mit->first.metaValueExists("peptide_sequence"))
          {
            description = String(mit->first.getMetaValue("peptide_sequence")).toQString();
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
            sub_item->setText(3, QString(chrom_description));
            if (! current_chromatogram.empty())
            {
              sub_item->setText(4, QString::number(current_chromatogram.front().getRT()));
              sub_item->setText(5, QString::number(current_chromatogram.back().getRT()));
            }

            sub_item->setText(6, MSChromatogram::ChromatogramNames[current_chromatogram.getChromatogramType()]);
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
      spectra_treewidget_->setHeaderLabels(QStringList() << "No peak map");
      spectra_treewidget_->setColumnCount(1); // needed, otherwise old column names for column 2, 3, etc are displayed
      has_data_ = false;
    }

    populateSearchBox_();

    if (more_than_one_spectrum && item != nullptr)
    { // not enabled
      item->setFlags(Qt::NoItemFlags);
    }

    // automatically set column width, depending on data
    spectra_treewidget_->header()->setStretchLastSection(false);
    spectra_treewidget_->header()->setSectionResizeMode(QHeaderView::ResizeToContents);
  }

  void SpectraViewWidget::populateSearchBox_()
  {
    const auto& header = spectra_treewidget_->headerItem();
    QStringList header_texts;
    for (int i = 0; i < header->columnCount(); ++i)
    {
      header_texts.push_back(header->text(i));
    }
    int current_index = spectra_combo_box_->currentIndex(); // when repainting we want the index to stay the same
    spectra_combo_box_->clear();
    spectra_combo_box_->addItems(header_texts);
    spectra_combo_box_->setCurrentIndex(current_index);
  }

  void SpectraViewWidget::clear()
  {
    getTreeWidget()->clear();
    getComboBox()->clear();
    has_data_ = false;
  }

}
