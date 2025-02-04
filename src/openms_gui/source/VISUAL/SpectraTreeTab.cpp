// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/SpectraTreeTab.h>

#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/LayerDataChrom.h>
#include <OpenMS/VISUAL/TreeView.h>

#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMenu>

namespace OpenMS
{

  std::vector<int> listToVec(const QList<QVariant>& in)
  {
    std::vector<int> out;
    for (const auto & i : in)
    {
      out.push_back(i.toInt());
    }
    return out;
  }


  QList<QVariant> vecToList(const std::vector<Size>& in)
  {
    QList<QVariant> res;
    for (Size i : in) res.push_back((unsigned int)i);
    return res;
  }


  // Use a namespace to encapsulate names, yet use c-style 'enum' for fast conversion to int.
  // So we can write: 'ClmnPeak::MS_LEVEL', but get implicit conversion to int
  namespace ClmnPeak
  {
    enum HeaderNames
    { // indices into QTableWidget's columns (which start at index 0)
      // note: make sure SPEC_INDEX remains at 1 (and is synced with ClmnChrom::CHROM_INDEX!!!)
      MS_LEVEL, SPEC_INDEX, RT, PRECURSOR_MZ, DISSOCIATION, SCANTYPE, ZOOM, /* last entry --> */ SIZE_OF_HEADERNAMES
    };
    // keep in SYNC with enum HeaderNames
    const QStringList HEADER_NAMES = QStringList()
       << "MS level" << "index" << "RT" << "precursor m/z" << "dissociation" << "scan" << "zoom";
  }
  // Use a namespace to encapsulate names, yet use c-style 'enum' for fast conversion to int.
  // So we can write: 'ClmnPeak::MS_LEVEL', but get implicit conversion to int
  namespace ClmnChrom
  {
    enum HeaderNames
    { // indices into QTableWidget's columns (which start at index 0)
      // note: make sure CHROM_INDEX remains at 1 (and is synced with ClmnPeak::SPEC_INDEX!!!)
      TYPE, CHROM_INDEX, MZ, DESCRIPTION, RT_START, RT_END, CHARGE, CHROM_TYPE, /* last entry --> */ SIZE_OF_HEADERNAMES
    };
    // keep in SYNC with enum HeaderNames
    const QStringList HEADER_NAMES = QStringList()
      << " type" << "index" << "m/z" << "Description" << "rt start" << "rt end" << "charge" << "chromatogram type";
  }

  struct IndexExtrator
  {
    explicit IndexExtrator(const QTreeWidgetItem* item)
      : spectrum_index(item->data(ClmnPeak::SPEC_INDEX, Qt::DisplayRole).toInt()),
        res(item->data(ClmnChrom::TYPE, Qt::UserRole).toList()) // this works, even if the QVariant is invalid (then the list is empty)
    {
    }

    bool hasChromIndices() const
    {
      return !res.empty();
    }

    const int spectrum_index;
    const QList<QVariant> res;
  };


  SpectraTreeTab::SpectraTreeTab(QWidget * parent) :
    QWidget(parent)
  {
    // these must be identical, because there is code which extracts the scan index irrespective of what we show
    assert(ClmnPeak::SPEC_INDEX == ClmnChrom::CHROM_INDEX);

    setObjectName("Scans");
    QVBoxLayout* spectra_widget_layout = new QVBoxLayout(this);
    spectra_treewidget_ = new TreeView(this);
    spectra_treewidget_->setWhatsThis("Spectrum selection bar<BR><BR>Here all spectra of the current experiment are shown. Left-click on a spectrum to show it. "
                                      "Double-clicking might be implemented as well, depending on the data. "
                                      "Context-menus for both the column header and data rows are available by right-clicking.");

    //~ no good for huge experiments - omitted:
    //~ spectrum_selection_->setSortingEnabled(true);
    //~ spectrum_selection_->sortByColumn ( 1, Qt::AscendingOrder);

    spectra_treewidget_->setDragEnabled(true);
    spectra_treewidget_->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(spectra_treewidget_, &QTreeWidget::currentItemChanged, this, &SpectraTreeTab::itemSelectionChange_);
    connect(spectra_treewidget_, &QTreeWidget::itemDoubleClicked, this, &SpectraTreeTab::itemDoubleClicked_);
    connect(spectra_treewidget_, &QTreeWidget::customContextMenuRequested, this, &SpectraTreeTab::spectrumContextMenu_);

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
    connect(spectra_search_box_, &QLineEdit::textEdited, this, &SpectraTreeTab::spectrumSearchText_);
    // .. show hit upon pressing Enter (internally we search again, since the user could have activated another layer with different selections after last search)
    connect(spectra_search_box_, &QLineEdit::returnPressed, this, &SpectraTreeTab::searchAndShow_);

    tmp_hbox_layout->addWidget(spectra_search_box_);
    tmp_hbox_layout->addWidget(spectra_combo_box_);
    spectra_widget_layout->addLayout(tmp_hbox_layout);
  }

  void SpectraTreeTab::spectrumSearchText_()
  {
    const QString& text = spectra_search_box_->text(); // get text from QLineEdit
    if (!text.isEmpty())
    {
      Qt::MatchFlags matchflags = Qt::MatchFixedString;
      matchflags |=  Qt::MatchRecursive; // match subitems (below top-level)
      // 'index' must be named identically for both data types
      assert(ClmnPeak::HEADER_NAMES[ClmnPeak::SPEC_INDEX] == ClmnChrom::HEADER_NAMES[ClmnChrom::CHROM_INDEX]);
      // ... for the following to work:
      if (spectra_combo_box_->currentText() != ClmnPeak::HEADER_NAMES[ClmnPeak::SPEC_INDEX])
      { // only the 'index' has to be matched exactly
        matchflags |= Qt::MatchStartsWith;
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

  bool SpectraTreeTab::getSelectedScan(MSExperiment& exp, LayerDataBase::DataType& current_type) const
  {
    exp.clear(true);
    QTreeWidgetItem* item = spectra_treewidget_->currentItem();
    if (item == nullptr)
    {
      return false;
    }
    // getting the index works for PEAK and CHROM data
    int index = item->data(ClmnPeak::SPEC_INDEX, Qt::DisplayRole).toInt();
    if (spectra_treewidget_->headerItem()->text(ClmnChrom::MZ) == ClmnChrom::HEADER_NAMES[ClmnChrom::MZ])
    { // we currently show chromatogram data
      current_type = LayerDataBase::DT_CHROMATOGRAM;
      exp.addChromatogram(last_peakmap_->getChromatograms()[index]);
    }
    else
    {
      current_type = LayerDataBase::DT_PEAK;
      exp.addSpectrum(last_peakmap_->getSpectra()[index]);
    }
    return true;
  }

  void SpectraTreeTab::itemSelectionChange_(QTreeWidgetItem* current, QTreeWidgetItem* /*previous*/)
  {
    if (current == nullptr)
    {
      return;
    }

    IndexExtrator ie(current);
    if (ie.hasChromIndices())
    { // open several chromatograms at once
      emit chromsSelected(listToVec(ie.res));
    }
    else
    {
      emit spectrumSelected(ie.spectrum_index);
    }
  }

  void SpectraTreeTab::searchAndShow_()
  {
    spectrumSearchText_(); // update selection first (we might be in a new layer)
    QList<QTreeWidgetItem*> selected = spectra_treewidget_->selectedItems();

    // show the first selected item
    if (!selected.empty())
    {
      itemSelectionChange_(selected.first(), selected.first());
    }
  }

  void SpectraTreeTab::itemDoubleClicked_(QTreeWidgetItem* current)
  {
    if (current == nullptr)
    {
      return;
    }
    IndexExtrator ie(current);
    if (!ie.hasChromIndices())
    {
      emit spectrumDoubleClicked(ie.spectrum_index);
    }
    else
    { // open several chromatograms at once
      emit chromsDoubleClicked(listToVec(ie.res));
    }
  }

  void SpectraTreeTab::spectrumContextMenu_(const QPoint& pos)
  {
    QTreeWidgetItem* item = spectra_treewidget_->itemAt(pos);
    if (item)
    {
      // create menu
      IndexExtrator ie(item);
      QMenu context_menu(spectra_treewidget_);
      context_menu.addAction("Show in 1D view", [&]()
      {
        if (!ie.hasChromIndices())
        {
          emit showSpectrumAsNew1D(ie.spectrum_index);
        }
        else
        { // open several chromatograms at once
          emit showChromatogramsAsNew1D(listToVec(ie.res));
        }
      });
      context_menu.addAction("Meta data", [&]() 
      {
        emit showSpectrumMetaData(ie.spectrum_index);
      });

      context_menu.exec(spectra_treewidget_->viewport()->mapToGlobal(pos));
    }
  }


  void populatePeakDataRow_(QTreeWidgetItem* item, const int index, const MSSpectrum& spec)
  {
    item->setText(ClmnPeak::MS_LEVEL, QString("MS") + QString::number(spec.getMSLevel()));
    item->setData(ClmnPeak::SPEC_INDEX, Qt::DisplayRole, index);
    item->setData(ClmnPeak::RT, Qt::DisplayRole, spec.getRT());

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
        item->setText(ClmnPeak::DISSOCIATION, ListUtils::concatenate(current_pc.getActivationMethodsAsString(), ",").toQString());
      }
      item->setData(ClmnPeak::PRECURSOR_MZ, Qt::DisplayRole, precursor_mz);
    }

    item->setText(ClmnPeak::SCANTYPE, QString::fromStdString(spec.getInstrumentSettings().NamesOfScanMode[spec.getInstrumentSettings().getScanMode()]));
    item->setText(ClmnPeak::ZOOM, (spec.getInstrumentSettings().getZoomScan() ? "yes" : "no"));
  }

  bool SpectraTreeTab::hasData(const LayerDataBase* layer)
  {
    if (layer == nullptr)
    {
      return false;
    }
    bool is_peak = layer->type == LayerDataBase::DT_PEAK;
    bool is_chrom = layer->type == LayerDataBase::DT_CHROMATOGRAM;
    return is_peak || is_chrom;
  }

  void SpectraTreeTab::updateEntries(LayerDataBase* layer)
  {
    if (layer == nullptr)
    {
      clear();
      return;
    }
    layer_ = layer;

    if (!spectra_treewidget_->isVisible() || spectra_treewidget_->signalsBlocked())
    {
      return;
    }
    
    spectra_treewidget_->blockSignals(true);
    RAIICleanup clean([&](){
      spectra_treewidget_->blockSignals(false); 
    });

    QTreeWidgetItem* toplevel_item = nullptr;
    QTreeWidgetItem* selected_item = nullptr;
    QList<QTreeWidgetItem*> toplevel_items;
    bool more_than_one_spectrum = true;

    // Branch if the current layer is a spectrum
    if (auto* lp = dynamic_cast<LayerDataPeak*>(layer))
    {
      Size spec_index = std::numeric_limits<Size>::max();
      if (auto* layer_peak1d = dynamic_cast<LayerData1DPeak*>(lp))
      {
        spec_index = layer_peak1d->getCurrentIndex();
      }
      const auto& cl = *lp;
      spectra_treewidget_->clear();

      std::vector<QTreeWidgetItem *> parent_stack;
      parent_stack.push_back(nullptr);
      bool fail = false;
      last_peakmap_ = &(*cl.getPeakData());
      spectra_treewidget_->setHeaders(ClmnPeak::HEADER_NAMES);

      for (Size i = 0; i < cl.getPeakData()->size(); ++i)
      {
        const MSSpectrum& current_spec = (*cl.getPeakData())[i];

        if (i > 0)
        {
          const MSSpectrum& prev_spec = (*cl.getPeakData())[i-1];
          // current MS level = previous MS level + 1 (e.g. current: MS2, previous: MS1)
          if (current_spec.getMSLevel() == prev_spec.getMSLevel() + 1)
          {
            toplevel_item = new QTreeWidgetItem(parent_stack.back());
            parent_stack.resize(parent_stack.size() + 1);
          }
          // current MS level = previous MS level (e.g. MS2,MS2 or MS1,MS1)
          else if (current_spec.getMSLevel() == prev_spec.getMSLevel())
          {
            if (parent_stack.size() == 1)
            {
              toplevel_item = new QTreeWidgetItem();
            }
            else
            {
              toplevel_item = new QTreeWidgetItem(*(parent_stack.end() - 2));
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
              toplevel_item = new QTreeWidgetItem(parent, parent_stack[parent_index + 1]);
            }
            else
            {
              toplevel_item = new QTreeWidgetItem((QTreeWidget *)nullptr);
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
          toplevel_item = new QTreeWidgetItem();
        }

        parent_stack.back() = toplevel_item;
        if (parent_stack.size() == 1)
        {
          toplevel_items.push_back(toplevel_item);
        }

        populatePeakDataRow_(toplevel_item, i, current_spec);

        if (i == spec_index)
        {
          // just remember it, select later
          selected_item = toplevel_item;
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
          toplevel_item = new QTreeWidgetItem();
          
          populatePeakDataRow_(toplevel_item, i, current_spec);

          toplevel_items.push_back(toplevel_item);
          if (i == spec_index)
          {
            // just remember it, select later
            selected_item = toplevel_item;
          }
        }
      }
      spectra_treewidget_->addTopLevelItems(toplevel_items);

      if (selected_item)
      {
        // now, select and scroll down to item
        selected_item->setSelected(true);
        // selected = for mouse clicks and multiselection,
        // current = for arrow navigation and signifies THE ONE active item
        spectra_treewidget_->setCurrentItem(selected_item);
        spectra_treewidget_->scrollToItem(selected_item);
      }
      if (cl.getPeakData()->size() > 1)
      {
        more_than_one_spectrum = false;
      }
    }
    // Branch if the current layer is a chromatogram (either indicated by its
    // type or by the flag which is set).
    else if (auto* lp = dynamic_cast<LayerDataChrom*>(layer))
    {
      const auto cl = *lp;
      LayerDataBase::ConstExperimentSharedPtrType exp = cl.getChromatogramData();
      
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
      if (!cl.getChromatogramData()->empty())
      {
        if (cl.getChromatogramData()->metaValueExists("multiple_select"))
        {
          multiple_select = cl.getChromatogramData()->getMetaValue("multiple_select").toBool();
        }
        if (cl.getChromatogramData()->metaValueExists("selected_chromatogram"))
        {
          this_selected_item = (int)cl.getChromatogramData()->getMetaValue("selected_chromatogram");
        }
      }
      
      // create a header list
      spectra_treewidget_->setHeaders(ClmnChrom::HEADER_NAMES);
           
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
        for (auto it = exp->getChromatograms().cbegin(); it != exp->getChromatograms().cend(); ++it)
        {
          map_precursor_to_chrom_idx[it->getPrecursor()].push_back(it - exp->getChromatograms().begin());
        }
      }

      int precursor_idx = 0;
      for (const auto& [pc, indx] : map_precursor_to_chrom_idx)
      {
        // Top level precursor entry
        toplevel_item = new QTreeWidgetItem();
        toplevel_item->setText(ClmnChrom::TYPE, "Peptide");
        toplevel_item->setData(ClmnChrom::TYPE, Qt::UserRole, vecToList(indx));
        toplevel_item->setData(ClmnChrom::CHROM_INDEX, Qt::DisplayRole, precursor_idx++);
        toplevel_item->setData(ClmnChrom::MZ, Qt::DisplayRole, pc.getMZ());
        // Show the peptide sequence if available, otherwise show the m/z and charge only
        QString description;
        if (pc.metaValueExists("peptide_sequence"))
        {
          description = String(pc.getMetaValue("peptide_sequence")).toQString();
        }
        else if (pc.metaValueExists("description"))
        {
          description = String(pc.getMetaValue("description")).toQString();
        }
        toplevel_item->setText(ClmnChrom::DESCRIPTION, description);
        toplevel_item->setData(ClmnChrom::CHARGE, Qt::DisplayRole, pc.getCharge());

        toplevel_items.push_back(toplevel_item);

        bool one_selected = false;
        // Show single chromatogram: iterate over all chromatograms corresponding to the current precursor and add action for the single chromatogram
        for (const Size chrom_idx : indx)
        {
          const MSChromatogram& current_chromatogram = exp->getChromatograms()[chrom_idx];

          // Children chromatogram entry
          QTreeWidgetItem* sub_item = new QTreeWidgetItem(toplevel_item);
          if ((int)chrom_idx == this_selected_item)
          {
            one_selected = true;
            selected_item = sub_item;
          }
          QString chrom_description = "ion";
          if (pc.metaValueExists("description"))
          {
            chrom_description = String(pc.getMetaValue("description")).toQString();
          }

          sub_item->setText(ClmnChrom::TYPE, "Transition");
          sub_item->setData(ClmnChrom::TYPE, Qt::UserRole, vecToList({chrom_idx}));
          sub_item->setData(ClmnChrom::CHROM_INDEX, Qt::DisplayRole, (unsigned int)chrom_idx);
          sub_item->setData(ClmnChrom::MZ, Qt::DisplayRole, current_chromatogram.getProduct().getMZ());
          sub_item->setText(ClmnChrom::DESCRIPTION, chrom_description);
          if (!current_chromatogram.empty())
          {
            sub_item->setData(ClmnChrom::RT_START, Qt::DisplayRole, current_chromatogram.front().getRT());
            sub_item->setData(ClmnChrom::RT_END, Qt::DisplayRole, current_chromatogram.back().getRT());
          }

          sub_item->setText(ClmnChrom::CHROM_TYPE, MSChromatogram::ChromatogramNames[current_chromatogram.getChromatogramType()]);
        }
        if (one_selected && multiple_select)
        {
          selected_item = toplevel_item;
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
      spectra_treewidget_->setHeaders(QStringList() << "No peak map");
    }

    populateSearchBox_();

    if (more_than_one_spectrum && toplevel_item != nullptr)
    { // not enabled
      toplevel_item->setFlags(Qt::NoItemFlags);
    }

    // automatically set column width, depending on data
    spectra_treewidget_->header()->setStretchLastSection(false);
    spectra_treewidget_->header()->setSectionResizeMode(QHeaderView::ResizeToContents);
  }

  void SpectraTreeTab::populateSearchBox_()
  {
    QStringList headers = spectra_treewidget_->getHeaderNames(WidgetHeader::WITH_INVISIBLE);
    int current_index = spectra_combo_box_->currentIndex(); // when repainting we want the index to stay the same
    spectra_combo_box_->clear();
    spectra_combo_box_->addItems(headers);
    spectra_combo_box_->setCurrentIndex(current_index);
  }

  void SpectraTreeTab::clear()
  {
    spectra_treewidget_->clear();
    spectra_combo_box_->clear();
  }




}
