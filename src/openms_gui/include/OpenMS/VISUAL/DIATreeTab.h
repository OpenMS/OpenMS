// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <QtWidgets>

#include <OpenMS/VISUAL/DataSelectionTabs.h>
#include <OpenMS/VISUAL/LayerDataBase.h>

class QLineEdit;
class QComboBox;
class QTreeWidget;
class QTreeWidgetItem;

namespace OpenMS
{
  class TreeView;
  struct OSWIndexTrace;

  /**
    @brief Hierarchical visualization and selection of spectra.

    @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI DIATreeTab :
    public QWidget, public DataTabBase
  {
    Q_OBJECT
  public:
    /// Constructor
    DIATreeTab(QWidget* parent = nullptr);
    /// Destructor
    ~DIATreeTab() override = default;

    // docu in base class
    bool hasData(const LayerDataBase* layer) override;

    /// refresh the table using data from @p cl
    /// @param cl Layer with OSW data; cannot be const, since we might read missing protein data from source on demand
    void updateEntries(LayerDataBase* cl) override;

    /// remove all visible data
    void clear() override;

  signals:
    /// emitted when a protein, peptide, feature or transition was selected
    void entityClicked(const OSWIndexTrace& trace);
    /// emitted when a protein, peptide, feature or transition was double-clicked
    void entityDoubleClicked(const OSWIndexTrace& trace);

  private:
    QLineEdit* spectra_search_box_ = nullptr;
    QComboBox* spectra_combo_box_ = nullptr;
    TreeView* dia_treewidget_ = nullptr;

    /// points to the data which is currently shown
    /// Useful to avoid useless repaintings, which would loose the open/close state of internal tree nodes and selected items
    OSWData* current_data_ = nullptr;

    /** 
      @brief convert a tree item to a pointer into an OSWData structure

      @param item The tree item (protein, peptide,...) that was clicked
      @return The index into the current OSWData @p current_data_
    **/
    OSWIndexTrace prepareSignal_(QTreeWidgetItem* item);

  private slots:
    /// fill the search-combo-box with current column header names
    void populateSearchBox_();
    /// searches for rows containing a search text (from spectra_search_box_); called when text search box is used
    void spectrumSearchText_();
    /// emits entityClicked() for all subitems
    void rowSelectionChange_(QTreeWidgetItem*, QTreeWidgetItem*);
    /// emits entityClicked() for all subitems
    void rowClicked_(QTreeWidgetItem*, int col);
    /// emits entityDoubleClicked() for all subitems
    void rowDoubleClicked_(QTreeWidgetItem*, int col);
    /// searches using text box and plots the spectrum
    void searchAndShow_();
  };
}

