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

#pragma once

#include <OpenMS/VISUAL/DataSelectionTabs.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/VISUAL/TableView.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <QtWidgets>
#include <QLineEdit>
#include <QComboBox>
#include <QTableWidget>
#include <QCheckBox>

namespace OpenMS
{
  /**
    @brief Tabular visualization / selection of identified spectra.

    @htmlinclude OpenMS_DigestSimulation.parameters
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
    bool hasData(const LayerData* layer) override;

    /// set layer data and create table anew; if given a nullptr, behaves as clear()
    void updateEntries(LayerData* model) override;
    /// get layer data
    LayerData* getLayer();

    /// clears all visible data from table widget and voids the layer
    void clear() override;

    /// Helper member to block outgoing signals
    bool ignore_update = false; 
  
  protected slots:
    /// Rebuild table entries
    void updateEntries_();
  signals:
    /// request to show a specific spectrum, and (if available) a specific pepId + pepHit in there (otherwise -1, -1)
    void spectrumSelected(int spectrum_index, int pep_id_index, int pep_hit_index);
    /// request to unshow a spectrum
    void spectrumDeselected(int spectrum_index);
    /// request to zoom into a 1D spec
    void requestVisibleArea1D(double lower_mz, double upper_mz);

  private:
   /// partially fill the bottom-most row  
   void fillRow_(const MSSpectrum& spectrum, const int spec_index, const QColor background_color);

    LayerData* layer_ = nullptr;
    QCheckBox* hide_no_identification_ = nullptr;
    QCheckBox* create_rows_for_commmon_metavalue_ = nullptr;
    TableView* table_widget_ = nullptr;
    QTableWidget* fragment_window_ = nullptr;
    bool is_ms1_shown_ = false;
  
  private slots:
    /// Saves the (potentially filtered) IDs as an idXML or mzIdentML file
    void saveIDs_();
    /// update PeptideIdentification / PeptideHits, when data in the table changes (status of checkboxes)
    void updatedSingleCell_(QTableWidgetItem* item);
    /// Cell clicked in table_widget; emits which spectrum (row) was clicked, and may show additional data
    void currentCellChanged_(int row, int column, int old_row, int old_column);
  };
}
