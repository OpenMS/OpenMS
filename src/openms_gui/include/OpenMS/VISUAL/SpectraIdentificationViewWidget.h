// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

#ifndef OPENMS_VISUAL_SPECTRAIDENTIFICATIONVIEWWIDGET_H
#define OPENMS_VISUAL_SPECTRAIDENTIFICATIONVIEWWIDGET_H

#include <QWidget>
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>
#include <QtGui/QTableWidget>
#include <QtGui/QCheckBox>

#include <OpenMS/VISUAL/LayerData.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

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
    SpectraIdentificationViewWidget(const Param & preferences, QWidget * parent = 0);
    /// Destructor
    virtual ~SpectraIdentificationViewWidget();
    /// Attach model
    void attachLayer(LayerData * model);
    /// Helper function to block outgoing signals
    bool ignore_update;

    // access the table widget
    QTableWidget* getTableWidget();
public slots:
    /// Rebuild table entries
    void updateEntries();
signals:
    void spectrumSelected(int);
    void spectrumDeselected(int);
    void spectrumDoubleClicked(int);
    void showSpectrumAs1D(int);
    void showSpectrumMetaData(int);
    void requestVisibleArea1D(double, double);
private:
    LayerData * layer_;
    QCheckBox * hide_no_identification_;
    QCheckBox * create_rows_for_commmon_metavalue_;
    QTableWidget * table_widget_;
    bool is_ms1_shown_;
private slots:
    /// Emits spectrumSelected with the current spectrum index
    void spectrumSelectionChange_(QTableWidgetItem *, QTableWidgetItem *);
    /// Export table entries as csv
    void exportEntries_();
    /// Saves the (potentially filtered) idXML
    void saveIdXML_();
    /// Display header context menu
    void headerContextMenu_(const QPoint &);
    /// Cell clicked in table_widget
    void cellClicked_(int row, int column);
  };
}

#endif // OPENMS_VISUAL_SPECTRAIDENTIFICATIONVIEWWIDGET_H
