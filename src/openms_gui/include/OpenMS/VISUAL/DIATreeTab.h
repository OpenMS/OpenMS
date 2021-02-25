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

#pragma once

#include <QtWidgets>

#include <OpenMS/VISUAL/DataSelectionTabs.h>
#include <OpenMS/VISUAL/LayerData.h>

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
    ~DIATreeTab() = default;

    // docu in base class
    bool hasData(const LayerData* layer) override;

    /// refresh the table using data from @p cl
    /// @param cl Layer with OSW data; cannot be const, since we might read missing protein data from source on demand
    void updateEntries(LayerData* cl) override;
    /// remove all visible data
    void clear();

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

