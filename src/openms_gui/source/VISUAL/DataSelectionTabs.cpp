// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

// OpenMS_GUI config
#include <OpenMS/VISUAL/DataSelectionTabs.h>

#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/VISUAL/DIATreeTab.h>
#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/SpectraTreeTab.h>
#include <OpenMS/VISUAL/SpectraIDViewTab.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>
#include <OpenMS/VISUAL/Plot2DCanvas.h>
#include <OpenMS/VISUAL/TVDIATreeTabController.h>
#include <OpenMS/VISUAL/TVSpectraViewController.h>
#include <OpenMS/VISUAL/TVIdentificationViewController.h>

namespace OpenMS
{
  /// enable and show the @p which tab


  /// double-click on disabled identification view
  /// --> enables it and creates an empty identification structure


  /// Default constructor

  DataSelectionTabs::DataSelectionTabs(QWidget* parent, TOPPViewBase* tv)
    : QTabWidget(parent),
    spectra_view_widget_(new SpectraTreeTab(this)),
    id_view_widget_(new SpectraIDViewTab(Param(), this)),
    dia_widget_(new DIATreeTab(this)),
    tab_ptrs_{ spectra_view_widget_, id_view_widget_, dia_widget_ },   // make sure to add new tabs here!
    spectraview_controller_(new TVSpectraViewController(tv)),
    idview_controller_(new TVIdentificationViewController(tv, id_view_widget_)),
    diatab_controller_(new TVDIATreeTabController(tv)),
    tv_(tv)
  {
    // Hook-up controller and views for spectra
    connect(spectra_view_widget_, &SpectraTreeTab::showSpectrumMetaData, tv, &TOPPViewBase::showSpectrumMetaData);
    connect(spectra_view_widget_, &SpectraTreeTab::showSpectrumAsNew1D, spectraview_controller_, &TVSpectraViewController::showSpectrumAsNew1D);
    connect(spectra_view_widget_, &SpectraTreeTab::showChromatogramsAsNew1D, spectraview_controller_, &TVSpectraViewController::showChromatogramsAsNew1D);
    connect(spectra_view_widget_, &SpectraTreeTab::spectrumSelected, spectraview_controller_, CONNECTCAST(TVSpectraViewController, activate1DSpectrum, (int)));
    connect(spectra_view_widget_, &SpectraTreeTab::chromsSelected, spectraview_controller_, CONNECTCAST(TVSpectraViewController, activate1DSpectrum, (const std::vector<int>&)));
    connect(spectra_view_widget_, &SpectraTreeTab::spectrumDoubleClicked, spectraview_controller_, &TVSpectraViewController::showSpectrumAsNew1D);
    connect(spectra_view_widget_, &SpectraTreeTab::chromsDoubleClicked, spectraview_controller_, &TVSpectraViewController::showChromatogramsAsNew1D);

    // Hook-up controller and views for identification
    connect(id_view_widget_, &SpectraIDViewTab::spectrumDeselected, idview_controller_, &TVIdentificationViewController::deactivate1DSpectrum);
    connect(id_view_widget_, &SpectraIDViewTab::spectrumSelected, idview_controller_, CONNECTCAST(TVIdentificationViewController, activate1DSpectrum, (int, int, int)));
    connect(id_view_widget_, &SpectraIDViewTab::requestVisibleArea1D, idview_controller_, &TVIdentificationViewController::setVisibleArea1D);

    // Hook-up controller and views for DIA
    connect(dia_widget_, &DIATreeTab::entityClicked, diatab_controller_, &TVDIATreeTabController::showChromatograms);
    connect(dia_widget_, &DIATreeTab::entityDoubleClicked, diatab_controller_, &TVDIATreeTabController::showChromatogramsAsNew1D);

    int index;
    index = addTab(spectra_view_widget_, spectra_view_widget_->objectName());
    if (index != SPECTRA_IDX)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Tab index is expected to be 0");
    }
    index = addTab(id_view_widget_, id_view_widget_->objectName());
    if (index != IDENT_IDX)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Tab index is expected to be 1");
    }
    index = addTab(dia_widget_, dia_widget_->objectName());
    if (index != DIAOSW_IDX)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Tab index is expected to be 2");
    }
    // make sure initialization was correct
    assert(tabBar()->count() == (int)tab_ptrs_.size());

    // switch between different view tabs
    connect(this, &QTabWidget::currentChanged, this, &DataSelectionTabs::currentTabChanged);
    connect(this, &QTabWidget::tabBarDoubleClicked, this, &DataSelectionTabs::tabBarDoubleClicked);
  }

  DataSelectionTabs::~DataSelectionTabs()
  {
    delete spectraview_controller_;
    delete idview_controller_;
    delete diatab_controller_;
  }

  LayerDataBase* getCurrentLayerData(TOPPViewBase* tv)
  {
    PlotCanvas* cc = tv->getActiveCanvas();
    if (cc == nullptr)
    {
      return nullptr;
    }
    if (cc->getCurrentLayerIndex() == Size(-1))
    {
      return nullptr;
    }
    return &(cc->getCurrentLayer());
  }

  // called externally
  // and internally by signals
  void DataSelectionTabs::callUpdateEntries()
  {
    // prevent infinite loop when calling 'setTabEnabled' -> currentTabChanged() -> update()
    this->blockSignals(true);
    RAIICleanup cleanup([&]()
    {
      this->blockSignals(false);
    });

    auto layer_ptr = getCurrentLayerData(tv_); // can be nullptr

    // becomes true if the currently visible tab has no data
    bool auto_select = false; 
    // the order is important here. On auto-select, we will pick the highest one which has data to show!
    Size highest_data_index = 0; // will pick spectra_view_widget_ if layer_ptr==nullptr
    for (Size i = 0; i < tab_ptrs_.size(); ++i)
    {
      auto widget = dynamic_cast<QWidget*>(tab_ptrs_[i]);
      bool has_data = tab_ptrs_[i]->hasData(layer_ptr);
      setTabEnabled(i, has_data); // enable/disable depending on data
      if (has_data)
      {
        highest_data_index = i;
      }
      if (!has_data && // the currently visible tab has no data --> select a new tab
          widget->isVisible())
      {
        auto_select = true;
      }
    }
    // pick the highest tab which has data
    if (auto_select)
    { 
      setCurrentIndex(highest_data_index);
    }
    Size current_index = currentIndex();
    // update the currently visible tab (might be disabled if no data is shown)
    tab_ptrs_[current_index]->updateEntries(layer_ptr);
  }

  void DataSelectionTabs::currentTabChanged(int tab_index)
  {
    // set new behavior
    switch (tab_index)
    {
    case SPECTRA_IDX:
      idview_controller_->deactivateBehavior(); // finalize old behavior
      diatab_controller_->deactivateBehavior();
      spectraview_controller_->activateBehavior(); // initialize new behavior
      break;
    case IDENT_IDX:
      spectraview_controller_->deactivateBehavior();
      diatab_controller_->deactivateBehavior();
      if (tv_->getActive2DWidget()) // currently, 2D window is open
      {
        idview_controller_->showSpectrumAsNew1D(0);
      }
      idview_controller_->activateBehavior();
      break;
    case DIAOSW_IDX:
      idview_controller_->deactivateBehavior(); // finalize old behavior
      spectraview_controller_->deactivateBehavior();
      diatab_controller_->activateBehavior(); // initialize new behavior
      break;
    default:
      std::cerr << "Error: tab_index " << tab_index << " is invalid\n";
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    callUpdateEntries(); //TODO actually this is overkill. Why would you load the entire table again
    // when you only switched tabs? The TabView should get notified when the layer data changes, so it only
    // updates when necessary...
    // The only thing that maybe needs to happen when switching tabs is to sync the index across the tables in the different tabs.
    // which is the only reason why we need to actually use callUpdateEntries here.
    // At least we reduced it to only updateEntries during tab switch, not EVERY update() [e.g. when resizing, refocussing...]
  }

  void DataSelectionTabs::showSpectrumAsNew1D(int index)
  {
    Plot1DWidget* widget_1d = tv_->getActive1DWidget();
    Plot2DWidget* widget_2d = tv_->getActive2DWidget();

    if (widget_1d || widget_2d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_controller_->showSpectrumAsNew1D(index);
      }

      if (id_view_widget_->isVisible())
      {
        idview_controller_->showSpectrumAsNew1D(index);
      }
    }
  }

  void DataSelectionTabs::showChromatogramsAsNew1D(const std::vector<int>& indices)
  {
    Plot1DWidget* widget_1d = tv_->getActive1DWidget();
    Plot2DWidget* widget_2d = tv_->getActive2DWidget();

    if (widget_1d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_controller_->showChromatogramsAsNew1D(indices);
      }
    }
    else if (widget_2d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_controller_->showChromatogramsAsNew1D(indices);
      }
    }
  }

  void DataSelectionTabs::tabBarDoubleClicked(int tab_index)
  {
    if (!tv_->getActivePlotWidget())
    {
      return;
    }
    switch (tab_index)
    {
    case IDENT_IDX:
      if (!isTabEnabled(IDENT_IDX))
      {
        setTabEnabled(IDENT_IDX, true); // enable identification view

        spectraview_controller_->deactivateBehavior();
        if (tv_->getActive2DWidget()) // currently 2D window is open
        {
          idview_controller_->showSpectrumAsNew1D(0);
        }
        idview_controller_->activateBehavior();

        // TODO: check this triggers update!
        setCurrentIndex(IDENT_IDX); // switch to identification view --> triggers currentTabChanged() slot
      }
    case SPECTRA_IDX:
    default:
      break;
    }

    // update here?
  }

  SpectraIDViewTab* DataSelectionTabs::getSpectraIDViewTab()
  {
    return id_view_widget_;
  }
} //namespace

