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
  class SpectraIdentificationViewWidget :
    public QWidget,
    public DefaultParamHandler
  {
    Q_OBJECT
public:
    /// Constructor
    SpectraIdentificationViewWidget(const Param& preferences, QWidget* parent = nullptr);
    /// Destructor
    ~SpectraIdentificationViewWidget() override = default;

    /// set layer data and create table anew
    void setLayer(LayerData* model);
    /// get layer data
    LayerData* getLayer();

    /// clears all visible data from table widget and void the layer
    void clear();

    /// Helper member to block outgoing signals
    bool ignore_update = false; 
  
  protected slots:
    /// Rebuild table entries
    void updateEntries();
signals:
    void spectrumSelected(int spectrum_index, int pep_id_index, int pep_hit_index);
    void spectrumDeselected(int spectrum_index);
    void spectrumDoubleClicked(int spectrum_index);
    void showSpectrumAs1D(int spectrum_index);
    void showSpectrumMetaData(int spectrum_index);
    void requestVisibleArea1D(double lower_mz, double upper_mz);
private:

   /// partially fill the bottom-most row  
   void fillRow_(const MSSpectrum& spectrum, int index, const QColor background_color);

    LayerData* layer_ = nullptr;
    QCheckBox* hide_no_identification_ = nullptr;
    QCheckBox* create_rows_for_commmon_metavalue_ = nullptr;
    TableView* table_widget_ = nullptr;
    QTableWidget* fragment_window_ = nullptr;
    bool is_ms1_shown_ = false;
private slots:
    /// Emits spectrumSelected with the current spectrum index
    void spectrumSelectionChange_(QTableWidgetItem*, QTableWidgetItem*);
    /// Saves the (potentially filtered) IDs as an idXML or mzIdentML file
    void saveIDs_();
    /// update PeptideIdentification / PeptideHits, when data in the table changes (status of checkboxes)
    void updatedSingleCell_(QTableWidgetItem* item);
    /// Cell clicked in table_widget
    void cellClicked_(int row, int column);
  };
}
