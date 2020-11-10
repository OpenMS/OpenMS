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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIATreeTab.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/VISUAL/TreeView.h>

#include <QtWidgets/QMenu>

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

  namespace Entity
  {
    enum Values
    {
      PROTEIN,
      PEPTIDE,
      FEATURE,
      TRANSITION,
      SIZE_OF_VALUES
    };
    const QStringList VALUES = { "protein", "peptide", "feature/peakgroup", "transition" };

  }
 
  /// given an item, goes up the tree to the root and collects indices in to the OSWData for each level
  struct IndexTrace
  {
    int idx_prot = -1;
    int idx_pep = -1;
    int idx_feat = -1;
    int idx_trans = -1;
    Entity::Values lowest = Entity::Values::SIZE_OF_VALUES;

    /// CTor which collects all the information
    IndexTrace(QTreeWidgetItem* current, const OSWData& data)
    {
      while (current != nullptr)
      {
        Entity::Values entity = Entity::Values(current->data(Clmn::INDEX, Qt::UserRole).toInt());
        int index = current->data(Clmn::INDEX, Qt::DisplayRole).toInt();
        
        if (lowest == Entity::Values::SIZE_OF_VALUES)
        { // set to level of first current
          lowest = entity;
        }
        switch (entity)
        {
          case Entity::Values::PROTEIN:
            idx_prot = index;
            break;
          case Entity::Values::PEPTIDE:
            idx_pep = index;
            break;
          case Entity::Values::FEATURE:
            idx_feat = index;  
            break;
          case Entity::Values::TRANSITION:
            idx_trans = index;
            break;
          default:
            throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
        // up one level
        current = current->parent();
      }
    }
  };



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

    connect(dia_treewidget_, &QTreeWidget::itemClicked, this, &DIATreeTab::rowSelectionChange2_);

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
    for (int idx_pep = 0; idx_pep < prot.getPeptidePrecursors().size(); ++idx_pep)
    {
      const auto& pep = prot.getPeptidePrecursors()[idx_pep];
      QTreeWidgetItem* item_pep = new QTreeWidgetItem(item_prot);
      item_pep->setData(Clmn::ENTITY, Qt::DisplayRole, Entity::VALUES[Entity::PEPTIDE]);
      item_pep->setData(Clmn::INDEX, Qt::DisplayRole, idx_pep);
      item_pep->setData(Clmn::INDEX, Qt::UserRole, Entity::PEPTIDE); // mark as peptide, so we know how to interpret the display role
      item_pep->setData(Clmn::CHARGE, Qt::DisplayRole, pep.getCharge());
      item_pep->setText(Clmn::FULL_NAME, pep.getSequence().c_str());

      for (int idx_feat = 0; idx_feat < pep.getFeatures().size(); ++idx_feat)
      {
        const auto& feat = pep.getFeatures()[idx_feat];
        QTreeWidgetItem* item_feat = new QTreeWidgetItem(item_pep);
        item_feat->setData(Clmn::ENTITY, Qt::DisplayRole, Entity::VALUES[Entity::FEATURE]);
        item_feat->setData(Clmn::INDEX, Qt::DisplayRole, idx_feat);
        item_feat->setData(Clmn::INDEX, Qt::UserRole, Entity::FEATURE); // mark as feature, so we know how to interpret the display role
        item_feat->setData(Clmn::RT_DELTA, Qt::DisplayRole, feat.getRTDelta());
        item_feat->setData(Clmn::QVALUE, Qt::DisplayRole, feat.getQValue());

        for (int idx_trans = 0; idx_trans < feat.getTransitionIDs().size(); ++idx_trans)
        {
          const uint trid = feat.getTransitionIDs()[idx_trans];
          QTreeWidgetItem* item_trans = new QTreeWidgetItem(item_feat);
          item_trans->setData(Clmn::ENTITY, Qt::DisplayRole, Entity::VALUES[Entity::TRANSITION]);
          item_trans->setData(Clmn::INDEX, Qt::DisplayRole, idx_trans);
          item_trans->setData(Clmn::INDEX, Qt::UserRole, Entity::TRANSITION); // mark as transition, so we know how to interpret the display role
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
    item_prot->setData(Clmn::INDEX, Qt::UserRole, Entity::PROTEIN); // mark as protein, so we know how to interpret the display role
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

  void DIATreeTab::rowSelectionChange_(QTreeWidgetItem* current, QTreeWidgetItem* previous)
  {
    /*	test for previous == 0 is important - without it,
        the wrong spectrum will be selected after finishing
        the execution of a TOPP tool on the whole data */
    if (current == nullptr || previous == nullptr)
    {
      return;
    }

    OSWData& data = *current_layer_->getChromatogramAnnotation().get();
    std::vector<int> transitions_to_show;

    IndexTrace tr(current, data);
    switch (tr.lowest)
    {
      case Entity::Values::PROTEIN:
        if (current->childCount() == 0)
        { // no peptides... load them
          OSWFile f(data.getSqlSourceFile());
          f.readProtein(data, tr.idx_prot);
        }
        fillProt(data.getProteins()[tr.idx_prot], current);
        // do nothing else -- showing all transitions for a protein is overwhelming...      
        break;
      case Entity::Values::PEPTIDE:
        {
          const auto& prot = data.getProteins()[tr.idx_prot];
          const auto& pep = prot.getPeptidePrecursors()[tr.idx_pep];
          for (const auto& feat : pep.getFeatures())
          {
            const auto& trids = feat.getTransitionIDs();
            transitions_to_show.insert(transitions_to_show.end(), trids.begin(), trids.end());
          }
          break;
        }
      case Entity::Values::FEATURE:
        {
          const auto& prot = data.getProteins()[tr.idx_prot];
          const auto& pep = prot.getPeptidePrecursors()[tr.idx_pep];
          const auto& feat = pep.getFeatures()[tr.idx_feat];
          const auto& trids = feat.getTransitionIDs();
          transitions_to_show.insert(transitions_to_show.end(), trids.begin(), trids.end());
          break;
        }
      case Entity::Values::TRANSITION:
        {
          const auto& prot = data.getProteins()[tr.idx_prot];
          const auto& pep = prot.getPeptidePrecursors()[tr.idx_pep];
          const auto& feat = pep.getFeatures()[tr.idx_feat];
          const auto& trid = feat.getTransitionIDs()[tr.idx_trans];
          transitions_to_show.insert(transitions_to_show.end(), trid);
          break;
        }
      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    std::cerr << "Showing transitions: " << ListUtils::concatenate(transitions_to_show, ", ") << "\n";
    emit transitionSelected(transitions_to_show);
  }

  void DIATreeTab::rowSelectionChange2_(QTreeWidgetItem* item, int col)
  {
    rowSelectionChange_(item, nullptr);
  }

  void DIATreeTab::searchAndShow_()
  {
    spectrumSearchText_(); // update selection first (we might be in a new layer)
    QList<QTreeWidgetItem*> selected = dia_treewidget_->selectedItems();
    // show the first selected item
    if (selected.size() > 0) rowSelectionChange_(selected.first(), selected.first());
  }


  void DIATreeTab::updateEntries(LayerData& cl)
  {
    if (!dia_treewidget_->isVisible() || dia_treewidget_->signalsBlocked())
    {
      return;
    }

    if (current_layer_ == &cl)
    {
      // layer data is still the same as last time ..
      // do not repopulate the table for now, since the data should not have changed
      // Note: If we ever need to redraw, the tree's state (which subtress are expanded, which items are selected) will need to be remembered and restored
      return;
    }

    // remember layer, because we need the OSWData from it, once the user wants to see transition plots...
    current_layer_ = &cl;

    dia_treewidget_->blockSignals(true);
    RAIICleanup clean([&]() { dia_treewidget_->blockSignals(false); });

    dia_treewidget_->clear();

    dia_treewidget_->setHeaders(Clmn::HEADER_NAMES);

    const OSWData* data = cl.getChromatogramAnnotation().get();

    if (data == nullptr  // DIA tab is active, but the layer has no data to show...
        || data->getProteins().empty())
    {
      dia_treewidget_->setHeaders(QStringList() << "No data");
    }
    else
    {
      for (int prot_index = 0; prot_index < data->getProteins().size(); ++prot_index)
      {
        const auto& prot = data->getProteins()[prot_index];
        auto item_prot = createProt(prot, prot_index);
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
  }



}
