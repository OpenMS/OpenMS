// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/DataSelectionTabs.h>
#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/TableView.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <QtWidgets>
#include <QCheckBox>

#include <unordered_map>
#include <vector>

namespace OpenMS
{
  /**
    @brief Tabular visualization / selection of identified spectra.

    @htmlinclude OpenMS_SpectraIDViewTab.parameters
  */
  class OPENMS_GUI_DLLAPI SpectraIDViewTab :
    public QWidget,
    public DefaultParamHandler,
    public DataTabBase
  {
    Q_OBJECT
  public:
    /// Constructor
    SpectraIDViewTab(const Param& preferences, QWidget* parent = nullptr);
    /// Destructor
    ~SpectraIDViewTab() override = default;

    // docu in base class
    bool hasData(const LayerDataBase* layer) override;

    /// set layer data and create table anew; if given a nullptr or the layer is not LayerDataPeak, behaves as clear()
    void updateEntries(LayerDataBase* model) override;
    /// get layer data
    LayerDataBase* getLayer();

    /// clears all visible data from table widget and voids the layer
    void clear() override;

    /// Helper member to block outgoing signals
    bool ignore_update = false; 
  
  protected slots:
    /// Rebuild table entries
    void updateEntries_();
    /// Rebuild protein table entries
    void updateProteinEntries_(int spec_cell_row_idx);
    /// Switch horizontal or vertical layout of the PSM and Proteintable 
    void switchOrientation_();
  signals:
    /// request to show a specific spectrum, and (if available) a specific pepId + pepHit in there (otherwise -1, -1)
    void spectrumSelected(int spectrum_index, int pep_id_index, int pep_hit_index);
    /// request to unshow a spectrum
    void spectrumDeselected(int spectrum_index);
    /// request to zoom into a 1D spec
    void requestVisibleArea1D(double lower_mz, double upper_mz);

  private:
    /// partially fill the bottom-most row  
    void fillRow_(const MSSpectrum& spectrum, const int spec_index, const QColor& background_color);
    /// extract the required part of the accession 
    static QString extractNumFromAccession_(const QString& listItem);
    /// open browser to navigate to uniport site with accession
    void openUniProtSiteWithAccession_(const QString& accession);

    class SelfResizingTableView_ : TableView
    {
      void resizeEvent(QResizeEvent * event) override;
    };

    LayerDataPeak* layer_ = nullptr;
    QCheckBox* hide_no_identification_ = nullptr;
    QCheckBox* create_rows_for_commmon_metavalue_ = nullptr;
    TableView* table_widget_ = nullptr;
    TableView* protein_table_widget_ = nullptr;
    QTableWidget* fragment_window_ = nullptr;
    QSplitter* tables_splitter_ = nullptr;
    bool is_first_time_loading_ = true;
    std::unordered_map<String, std::vector<const PeptideIdentification*>> protein_to_peptide_id_map;

  private slots:
    /// Saves the (potentially filtered) IDs as an idXML or mzIdentML file
    void saveIDs_();
    /// update PeptideIdentification / PeptideHits, when data in the table changes (status of checkboxes)
    void updatedSingleCell_(QTableWidgetItem* item);
    /// Cell clicked in table_widget; emits which spectrum (row) was clicked, and may show additional data
    void currentCellChanged_(int row, int column, int old_row, int old_column);

    /// Create 'protein accession to peptide identification' map using C++ STL unordered_map
    void createProteinToPeptideIDMap_();

    /// Cell selected or deselected: this is only used to check for deselection, rest happens in currentCellChanged_
    void currentSpectraSelectionChanged_();

    /// update ProteinHits, when data in the table changes (status of checkboxes)
    void updatedSingleProteinCell_(QTableWidgetItem* /*item*/);
    /// Protein Cell clicked in protein_table_widget; emits which protein (row) was clicked, and may show additional data
    void proteinCellClicked_(int row, int column);
  };
}
