// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
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
  /**
    @brief Hierarchical visualization and selection of spectra.

    @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI SpectraTreeTab :
    public QWidget,
    public DataTabBase
  {
    Q_OBJECT
public:
    /// Constructor
    SpectraTreeTab(QWidget * parent = nullptr);

    /// Destructor
    ~SpectraTreeTab() override = default;

    /// docu in base class
    bool hasData(const LayerDataBase* layer) override;

    /// refresh the table using data from @p cl
    void updateEntries(LayerDataBase* cl) override;
    
    /// remove all visible data
    void clear() override;

    /// Return a copy of the currently selected spectrum/chrom (for drag'n'drop to new window)
    /// and store it either as Spectrum or Chromatogram in @p exp (all other data is cleared)
    /// If no spectrum/chrom is selected, false is returned and @p exp is empty
    /// 
    /// @param[out] exp The currently active spec/chrom
    /// @param current_type Either DT_PEAK or DT_CHROMATOGRAM, depending on what is currently shown
    /// @return true if a spec/chrom is currently active
    bool getSelectedScan(MSExperiment& exp, LayerDataBase::DataType& current_type) const;

signals:
    void spectrumSelected(int);
    void chromsSelected(std::vector<int> indices);
    void spectrumDoubleClicked(int);
    void chromsDoubleClicked(std::vector<int> indices);
    void showSpectrumAsNew1D(int);
    void showChromatogramsAsNew1D(std::vector<int> indices);
    void showSpectrumMetaData(int);
private:
    QLineEdit* spectra_search_box_ = nullptr;
    QComboBox* spectra_combo_box_ = nullptr;
    TreeView* spectra_treewidget_ = nullptr;
    LayerDataBase* layer_ = nullptr;
    /// cache to store mapping of chromatogram precursors to chromatogram indices
    std::map<size_t, std::map<Precursor, std::vector<Size>, Precursor::MZLess> > map_precursor_to_chrom_idx_cache_;
    /// remember the last PeakMap that we used to fill the spectra list (and avoid rebuilding it)
    const PeakMap* last_peakmap_ = nullptr;

private slots:

    /// fill the search-combo-box with current column header names
    void populateSearchBox_();
    /// searches for rows containing a search text (from spectra_search_box_); called when text search box is used
    void spectrumSearchText_();
    /// emits spectrumSelected() for PEAK or chromsSelected() for CHROM data
    void itemSelectionChange_(QTreeWidgetItem *, QTreeWidgetItem *);
    /// searches using text box and plots the spectrum
    void searchAndShow_(); 
    /// called upon double click on an item; emits spectrumDoubleClicked() or chromsDoubleClicked() after some checking
    void itemDoubleClicked_(QTreeWidgetItem *); 
    /// Display context menu; allows to open metadata window
    void spectrumContextMenu_(const QPoint &);
  };
}

