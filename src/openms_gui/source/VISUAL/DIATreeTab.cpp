// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/DIATreeTab.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/VISUAL/LayerDataChrom.h>
#include <OpenMS/VISUAL/TreeView.h>

namespace OpenMS
{
  // Use a namespace to encapsulate names, yet use c-style 'enum' for fast conversion to int.
  // So we can write: 'Clmn::MS_LEVEL', but get implicit conversion to int
  namespace Clmn
  {
    enum HeaderNames
    { // indices into QTableWidget's columns (which start at index 0)
      ENTITY, INDEX, CHARGE, FULL_NAME, RT_DELTA, QVALUE, /* last entry --> */ SIZE_OF_HEADERNAMES
    };
    // keep in SYNC with enum HeaderNames
    const QStringList HEADER_NAMES = QStringList()
      << "entity" << "index" << "charge" << "full name" << "rt delta" << "q-value";
  }
                                       
  /// given an item, goes up the tree to the root and collects indices in to the OSWData for each level
  OSWIndexTrace getTrace(QTreeWidgetItem* current)
  {
    OSWIndexTrace trace;

    while (current != nullptr)
    {
      OSWHierarchy::Level level = OSWHierarchy::Level(current->data(Clmn::INDEX, Qt::UserRole).toInt());
      int index = current->data(Clmn::INDEX, Qt::DisplayRole).toInt();

      if (trace.lowest == OSWHierarchy::Level::SIZE_OF_VALUES)
      { // set to level of first current
        trace.lowest = level;
      }
      switch (level)
      {
      case OSWHierarchy::Level::PROTEIN:
        trace.idx_prot = index;
        break;
      case OSWHierarchy::Level::PEPTIDE:
        trace.idx_pep = index;
        break;
      case OSWHierarchy::Level::FEATURE:
        trace.idx_feat = index;
        break;
      case OSWHierarchy::Level::TRANSITION:
        trace.idx_trans = index;
        break;
      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
      // up one level
      current = current->parent();
    }
    return trace;
  }

  DIATreeTab::DIATreeTab(QWidget* parent) :
    QWidget(parent)
  {
    setObjectName("DIA OSW View");
    QVBoxLayout* spectra_widget_layout = new QVBoxLayout(this);
    dia_treewidget_ = new TreeView(this);
    dia_treewidget_->setWhatsThis("Protein/Peptide/Transition selection bar<BR><BR>Here all XICs of a DIA experiment are shown. Left-click on a chrom to show it. "
      "Double-clicking might be implemented as well, depending on the data. "
      "Context-menus for both the column header and data rows are available by right-clicking.");

    //~ no good for huge experiments - omitted:
    //~ spectrum_selection_->setSortingEnabled(true);
    //~ spectrum_selection_->sortByColumn ( 1, Qt::AscendingOrder);

    dia_treewidget_->setDragEnabled(true);
    dia_treewidget_->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(dia_treewidget_, &QTreeWidget::currentItemChanged, this, &DIATreeTab::rowSelectionChange_);
    connect(dia_treewidget_, &QTreeWidget::itemClicked, this, &DIATreeTab::rowClicked_);
    connect(dia_treewidget_, &QTreeWidget::itemDoubleClicked, this, &DIATreeTab::rowDoubleClicked_);

    spectra_widget_layout->addWidget(dia_treewidget_);

    QHBoxLayout* tmp_hbox_layout = new QHBoxLayout();

    spectra_search_box_ = new QLineEdit(this);
    spectra_search_box_->setPlaceholderText("<search text>");
    spectra_search_box_->setWhatsThis("Search in a certain column. Hits are shown as you type. Press <Enter> to display the first hit.");
    spectra_search_box_->setToolTip(spectra_search_box_->whatsThis());

    spectra_combo_box_ = new QComboBox(this);
    spectra_combo_box_->setWhatsThis("Sets the column in which to search.");
    spectra_combo_box_->setToolTip(spectra_combo_box_->whatsThis());


    // search whenever text is typed (and highlight the hits)
    connect(spectra_search_box_, &QLineEdit::textEdited, this, &DIATreeTab::spectrumSearchText_);
    // .. show hit upon pressing Enter (internally we search again, since the user could have activated another layer with different selections after last search)
    connect(spectra_search_box_, &QLineEdit::returnPressed, this, &DIATreeTab::searchAndShow_);

    tmp_hbox_layout->addWidget(spectra_search_box_);
    tmp_hbox_layout->addWidget(spectra_combo_box_);
    spectra_widget_layout->addLayout(tmp_hbox_layout);
  }

  /// adds a subtree (with peptides ...) to a given protein 
  void fillProt(const OSWProtein& prot, QTreeWidgetItem* item_prot)
  {
    for (size_t idx_pep = 0; idx_pep < prot.getPeptidePrecursors().size(); ++idx_pep)
    {
      const auto& pep = prot.getPeptidePrecursors()[idx_pep];
      QTreeWidgetItem* item_pep = new QTreeWidgetItem(item_prot);
      item_pep->setData(Clmn::ENTITY, Qt::DisplayRole, OSWHierarchy::LevelName[OSWHierarchy::PEPTIDE]);
      item_pep->setData(Clmn::INDEX, Qt::DisplayRole, (int)idx_pep);
      item_pep->setData(Clmn::INDEX, Qt::UserRole, OSWHierarchy::PEPTIDE); // mark as peptide, so we know how to interpret the display role
      item_pep->setData(Clmn::CHARGE, Qt::DisplayRole, pep.getCharge());
      item_pep->setText(Clmn::FULL_NAME, pep.getSequence().c_str());

      for (size_t idx_feat = 0; idx_feat < pep.getFeatures().size(); ++idx_feat)
      {
        const auto& feat = pep.getFeatures()[idx_feat];
        QTreeWidgetItem* item_feat = new QTreeWidgetItem(item_pep);
        item_feat->setData(Clmn::ENTITY, Qt::DisplayRole, OSWHierarchy::LevelName[OSWHierarchy::FEATURE]);
        item_feat->setData(Clmn::INDEX, Qt::DisplayRole, (int)idx_feat);
        item_feat->setData(Clmn::INDEX, Qt::UserRole, OSWHierarchy::FEATURE); // mark as feature, so we know how to interpret the display role
        item_feat->setData(Clmn::RT_DELTA, Qt::DisplayRole, feat.getRTDelta());
        item_feat->setData(Clmn::QVALUE, Qt::DisplayRole, feat.getQValue());

        for (size_t idx_trans = 0; idx_trans < feat.getTransitionIDs().size(); ++idx_trans)
        {
          QTreeWidgetItem* item_trans = new QTreeWidgetItem(item_feat);
          item_trans->setData(Clmn::ENTITY, Qt::DisplayRole, OSWHierarchy::LevelName[OSWHierarchy::TRANSITION]);
          item_trans->setData(Clmn::INDEX, Qt::DisplayRole, (int)idx_trans);
          item_trans->setData(Clmn::INDEX, Qt::UserRole, OSWHierarchy::TRANSITION); // mark as transition, so we know how to interpret the display role
        }
      }
      //item_prot->addChild(item_pep);
    }
  }

  /// creates a protein subtree (with peptides etc, if available)
  QTreeWidgetItem* createProt(const OSWProtein& prot, int prot_index)
  {
    QTreeWidgetItem* item_prot = new QTreeWidgetItem();
    item_prot->setData(Clmn::ENTITY, Qt::DisplayRole, "protein");
    item_prot->setData(Clmn::INDEX, Qt::DisplayRole, prot_index);
    item_prot->setData(Clmn::INDEX, Qt::UserRole, OSWHierarchy::PROTEIN); // mark as protein, so we know how to interpret the display role
    item_prot->setText(Clmn::FULL_NAME, prot.getAccession().c_str());

    // if possible, fill it already
    fillProt(prot, item_prot);
    
    return item_prot;
  }


  void DIATreeTab::spectrumSearchText_()
  {
    const QString& text = spectra_search_box_->text(); // get text from QLineEdit
    if (!text.isEmpty())
    {
      Qt::MatchFlags matchflags = Qt::MatchFixedString;
      matchflags |= Qt::MatchRecursive; // match subitems (below top-level)
      matchflags |= Qt::MatchStartsWith;

      QList<QTreeWidgetItem*> searched = dia_treewidget_->findItems(text, matchflags, spectra_combo_box_->currentIndex());

      if (!searched.isEmpty())
      {
        dia_treewidget_->clearSelection();
        searched.first()->setSelected(true);
        dia_treewidget_->update();
        dia_treewidget_->scrollToItem(searched.first());
      }
    }
  }


  OSWIndexTrace DIATreeTab::prepareSignal_(QTreeWidgetItem* item)
  {
    OSWIndexTrace tr;
    if (item == nullptr || current_data_ == nullptr)
    {
      return tr;
    }

    tr = getTrace(item);
    switch (tr.lowest)
    {
    case OSWHierarchy::Level::PROTEIN:
      if (item->childCount() == 0)
      { // no peptides... load them
        OSWFile f(current_data_->getSqlSourceFile());
        f.readProtein(*current_data_, tr.idx_prot);
        fillProt(current_data_->getProteins()[tr.idx_prot], item);
      }
      // do nothing else -- showing all transitions for a protein is overwhelming...      
      break;
    case OSWHierarchy::Level::PEPTIDE:
    case OSWHierarchy::Level::FEATURE:
    case OSWHierarchy::Level::TRANSITION:
      break;
    default:
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    return tr;
  }

  void DIATreeTab::rowSelectionChange_(QTreeWidgetItem* current, QTreeWidgetItem* /*previous*/)
  {
    auto tr = prepareSignal_(current);
    if (tr.isSet()) emit entityClicked(tr);
  }

  void DIATreeTab::rowClicked_(QTreeWidgetItem* item, int /*col*/)
  {
    auto tr = prepareSignal_(item);
    if (tr.isSet()) emit entityClicked(tr);
  }

  void DIATreeTab::rowDoubleClicked_(QTreeWidgetItem* item, int /*col*/)
  {
    auto tr = prepareSignal_(item);
    if (tr.isSet())
    {
      entityDoubleClicked(tr);
    }
  }

  void DIATreeTab::searchAndShow_()
  {
    spectrumSearchText_(); // update selection first (we might be in a new layer)
    QList<QTreeWidgetItem*> selected = dia_treewidget_->selectedItems();
    // show the first selected item
    if (selected.size() > 0) rowSelectionChange_(selected.first(), selected.first());
  }


  bool DIATreeTab::hasData(const LayerDataBase* layer)
  {
    if (auto* lp = dynamic_cast<const LayerDataChrom*>(layer))
    {
      OSWData* data = lp->getChromatogramAnnotation().get();
      return (data != nullptr && !data->getProteins().empty());
    }
    return false;
  }

  void DIATreeTab::updateEntries(LayerDataBase* layer)
  {
    if (layer == nullptr)
    {
      clear();
      return;
    }

    if (!dia_treewidget_->isVisible() || dia_treewidget_->signalsBlocked())
    {
      return;
    }
    auto* lp = dynamic_cast<const LayerDataChrom*>(layer);
    if (!lp) return;

    OSWData* data = lp->getChromatogramAnnotation().get();

    if (current_data_ == data)
    {
      // layer data is still the same as last time ..
      // do not repopulate the table for now, since the data should not have changed
      // Note: If we ever need to redraw, the tree's state (which subtress are expanded, which items are selected) will need to be remembered and restored
      return;
    }

    // update last data pointer
    current_data_ = data;

    dia_treewidget_->blockSignals(true);
    RAIICleanup clean([&]() { dia_treewidget_->blockSignals(false); });

    dia_treewidget_->clear();

    dia_treewidget_->setHeaders(Clmn::HEADER_NAMES);


    if (data == nullptr  // DIA tab is active, but the layer has no data to show...
        || data->getProteins().empty())
    {
      dia_treewidget_->setHeaders(QStringList() << "No data");
    }
    else
    {
      for (size_t prot_index = 0; prot_index < data->getProteins().size(); ++prot_index)
      {
        const auto& prot = data->getProteins()[prot_index];
        auto item_prot = createProt(prot, (int)prot_index);
        dia_treewidget_->addTopLevelItem(item_prot);
      }
    }
    
    populateSearchBox_();

    // automatically set column width, depending on data
    dia_treewidget_->header()->setStretchLastSection(false);
    dia_treewidget_->header()->setSectionResizeMode(QHeaderView::ResizeToContents);
  }

  void DIATreeTab::populateSearchBox_()
  {
    QStringList headers = dia_treewidget_->getHeaderNames(WidgetHeader::WITH_INVISIBLE);
    int current_index = spectra_combo_box_->currentIndex(); // when repainting we want the index to stay the same
    spectra_combo_box_->clear();
    spectra_combo_box_->addItems(headers);
    spectra_combo_box_->setCurrentIndex(current_index);
  }

  void DIATreeTab::clear()
  {
    dia_treewidget_->clear();
    spectra_combo_box_->clear();
    current_data_ = nullptr;
  }

}
