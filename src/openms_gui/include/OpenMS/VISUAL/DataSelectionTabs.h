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

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/KERNEL/StandardTypes.h>

#include <QTabWidget>

namespace OpenMS
{
  class DIATreeTab;
  class LayerData;
  class SpectraTreeTab;
  class SpectraIDViewTab;
  class TVDIATreeTabController;
  class TVIdentificationViewController;
  class TVSpectraViewController;
  class TOPPViewBase;

  /// all tabs need to implement this interface
  class OPENMS_GUI_DLLAPI DataTabBase
  {
  public:
    /// given a layer, determine if the tab could use it to show data (useful to decide if the tab should be enabled/disabled)
    /// If a nullptr is given, it HAS to return false!
    virtual bool hasData(const LayerData* layer) = 0;

    /// populate the tab using date from @p layer
    /// Should handle nullptr well (by calling clear())
    virtual void updateEntries(LayerData* layer) = 0;

    /// explicitly show no data at all
    virtual void clear() = 0;
  };

  /**
    @brief A tabbed view, to browse lists of spectra or identifications
    
  */
  class OPENMS_GUI_DLLAPI DataSelectionTabs
    : public QTabWidget
  {
    Q_OBJECT

  public:
    enum TAB_INDEX
    {
      SPECTRA_IDX = 0,  ///< first tab
      IDENT_IDX = 1,    ///< second tab
      DIAOSW_IDX = 2,   ///< third tab
      SIZE_OF_TAB_INDEX    
    };

    /// Default constructor
    DataSelectionTabs(QWidget* parent, TOPPViewBase* tv);

    /// Update items in the tabs according to the currently selected layer.
    /// Tabs which have data to show are automatically enabled. Others are disabled.
    /// If the currently visible tab would have to data to show, we pick the highest (rightmost) tab
    /// which has data and show that instead
    void update();

    /// invoked when user changes the active tab to @p tab_index
    void currentTabChanged(int tab_index);
    
    /// forwards to the TOPPView*Behaviour classes, to show a certain spectrum in 1D
    void showSpectrumAsNew1D(int index);
    
    /// forwards to the TOPPView*Behaviour classes, to show a certain set of chromatograms in 1D
    void showChromatogramsAsNew1D(const std::vector<int>& indices);

    /// double-click on disabled identification view
    /// --> enables it and creates an empty identification structure
    void tabBarDoubleClicked(int tab_index);

    SpectraIDViewTab* getSpectraIDViewTab();
  signals:

  private:
    ///@name Spectrum selection widgets
    //@{
    SpectraTreeTab* spectra_view_widget_;
    SpectraIDViewTab* id_view_widget_;
    DIATreeTab* dia_widget_;
    //@}
    std::vector< DataTabBase* > tab_ptrs_; ///< holds pointers to all of the above tabs, for iteration purposes

    /// TOPPView behavior for the spectra view
    TVSpectraViewController* spectraview_controller_;
    /// TOPPView behavior for the identification view
    TVIdentificationViewController* idview_controller_;
    /// TOPPView behavior for the DIA view
    TVDIATreeTabController* diatab_controller_;

    /// pointer to base class to access some members (going signal/slot would be cleaner)
    TOPPViewBase* tv_;
  };

} //namespace
