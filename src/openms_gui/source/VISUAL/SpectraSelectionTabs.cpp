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
#include <OpenMS/VISUAL/SpectraSelectionTabs.h>

#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/SpectraViewWidget.h>
#include <OpenMS/VISUAL/SpectraIdentificationViewWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/TOPPViewSpectraViewBehavior.h>
#include <OpenMS/VISUAL/TOPPViewIdentificationViewBehavior.h>

namespace OpenMS
{
  /// enable and show the @p which tab


  /// double-click on disabled identification view
  /// --> enables it and creates an empty identification structure


  /// Default constructor

  inline SpectraSelectionTabs::SpectraSelectionTabs(QWidget* parent, TOPPViewBase* tv)
    : QTabWidget(parent),
    spectra_view_widget_(new SpectraViewWidget(this)),
    id_view_widget_(new SpectraIdentificationViewWidget(Param(), this)),
    spectraview_behavior_(new TOPPViewSpectraViewBehavior(tv)),
    idview_behaviour_(new TOPPViewIdentificationViewBehavior(tv, id_view_widget_)),
    tv_(tv)
  {
    // Hook-up controller and views for spectra inspection
    connect(spectra_view_widget_, &SpectraViewWidget::showSpectrumMetaData, tv, &TOPPViewBase::showSpectrumMetaData);
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, showSpectrumAs1D, (int)), this, CONNECTCAST(SpectraSelectionTabs, showSpectrumAs1D, (int)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, showSpectrumAs1D, (std::vector<int>)), this, CONNECTCAST(SpectraSelectionTabs, showSpectrumAs1D, (std::vector<int>)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, spectrumSelected, (int)), spectraview_behavior_, CONNECTCAST(TOPPViewSpectraViewBehavior, activate1DSpectrum, (int)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, spectrumSelected, (std::vector<int>)), spectraview_behavior_, CONNECTCAST(TOPPViewSpectraViewBehavior, activate1DSpectrum, (const std::vector<int>&)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, spectrumDoubleClicked, (int)), this, CONNECTCAST(SpectraSelectionTabs, showSpectrumAs1D, (int)));
    connect(spectra_view_widget_, CONNECTCAST(SpectraViewWidget, spectrumDoubleClicked, (std::vector<int>)), this, CONNECTCAST(SpectraSelectionTabs, showSpectrumAs1D, (std::vector<int>)));

    // Hook-up controller and views for identification inspection
    connect(id_view_widget_, &SpectraIdentificationViewWidget::spectrumDeselected, idview_behaviour_, &TOPPViewIdentificationViewBehavior::deactivate1DSpectrum);
    connect(id_view_widget_, &SpectraIdentificationViewWidget::showSpectrumAs1D, this, CONNECTCAST(SpectraSelectionTabs, showSpectrumAs1D, (int)));
    connect(id_view_widget_, &SpectraIdentificationViewWidget::spectrumSelected, idview_behaviour_, CONNECTCAST(TOPPViewIdentificationViewBehavior, activate1DSpectrum, (int, int, int)));
    connect(id_view_widget_, &SpectraIdentificationViewWidget::requestVisibleArea1D, idview_behaviour_, &TOPPViewIdentificationViewBehavior::setVisibleArea1D);

    int index;
    index = addTab(spectra_view_widget_, spectra_view_widget_->objectName());
    OPENMS_PRECONDITION(index == SPECTRA_IDX, "Tab index is expected to be 0");
    index = addTab(id_view_widget_, id_view_widget_->objectName());
    OPENMS_PRECONDITION(index == IDENT_IDX, "Tab index is expected to be 1");
    setTabEnabled(SPECTRA_IDX, false);
    setTabEnabled(IDENT_IDX, false);

    // switch between different view tabs
    connect(this, &QTabWidget::currentChanged, this, &SpectraSelectionTabs::currentTabChanged);
    connect(this, &QTabWidget::tabBarDoubleClicked, this, &SpectraSelectionTabs::tabBarDoubleClicked);
  }

  inline void SpectraSelectionTabs::update()
  {
    // prevent infinite loop when calling 'setTabEnabled' -> currentTabChanged() -> update()
    this->blockSignals(true);
    RAIICleanup cleanup([&]()
    {
      this->blockSignals(false);
    });

    SpectrumCanvas* cc = tv_->getActiveCanvas();
    int layer_row = (cc == nullptr ? -1 : (int)cc->activeLayerIndex());

    if (cc == nullptr || layer_row == -1)
    {
      spectra_view_widget_->clear();
      id_view_widget_->clear();
      setTabEnabled(SPECTRA_IDX, true);
      setTabEnabled(IDENT_IDX, false);
      return;
    }


    if (spectra_view_widget_->isVisible())
    {
      spectra_view_widget_->updateEntries(cc->getCurrentLayer());
    }

    if (id_view_widget_->isVisible())
    {
      if (&cc->getCurrentLayer() != id_view_widget_->getLayer())
      {
        id_view_widget_->setLayer(&cc->getCurrentLayer());
      }
    }


  }

  inline void SpectraSelectionTabs::currentTabChanged(int tab_index)
  {
    // set new behavior
    switch (tab_index)
    {
    case SPECTRA_IDX:
      idview_behaviour_->deactivateBehavior(); // finalize old behavior
      spectraview_behavior_->activateBehavior(); // initialize new behavior
      break;
    case IDENT_IDX:
      spectraview_behavior_->deactivateBehavior();
      if (tv_->getActive2DWidget()) // currently 2D window is open
      {
        showSpectrumAs1D(0);
      }
      idview_behaviour_->activateBehavior();
      break;
    default:
      std::cerr << "Error: tab_index " << tab_index << endl;
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    update();
  }

  inline void SpectraSelectionTabs::showSpectrumAs1D(int index)
  {
    Spectrum1DWidget* widget_1d = tv_->getActive1DWidget();
    Spectrum2DWidget* widget_2d = tv_->getActive2DWidget();

    if (widget_1d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_->showSpectrumAs1D(index);
      }

      if (id_view_widget_->isVisible())
      {
        idview_behaviour_->showSpectrumAs1D(index);
      }
    }
    else if (widget_2d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_->showSpectrumAs1D(index);
      }

      if (id_view_widget_->isVisible())
      {
        idview_behaviour_->showSpectrumAs1D(index);
      }
    }
  }

  inline void SpectraSelectionTabs::showSpectrumAs1D(std::vector<int> indices)
  {
    Spectrum1DWidget* widget_1d = tv_->getActive1DWidget();
    Spectrum2DWidget* widget_2d = tv_->getActive2DWidget();

    if (widget_1d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_->showSpectrumAs1D(indices);
      }
    }
    else if (widget_2d)
    {
      if (spectra_view_widget_->isVisible())
      {
        spectraview_behavior_->showSpectrumAs1D(indices);
      }
    }
  }

  inline void SpectraSelectionTabs::setTabEnabled(int index, bool b)
  {
    if (index == 0 && b == false)
    {
      std::cerr << "HA!";
    }
    QTabWidget::setTabEnabled(index, b);
  }

  inline void SpectraSelectionTabs::tabBarDoubleClicked(int tab_index)
  {
    if (!tv_->getActiveSpectrumWidget()) return;

    switch (tab_index)
    {
    case IDENT_IDX:
      if (!isTabEnabled(IDENT_IDX))
      {
        setTabEnabled(IDENT_IDX, true); // enable identification view

        spectraview_behavior_->deactivateBehavior();
        if (tv_->getActive2DWidget()) // currently 2D window is open
        {
          showSpectrumAs1D(0);
        }
        idview_behaviour_->activateBehavior();

        // TODO: check this triggers update!
        setCurrentIndex(IDENT_IDX); // switch to identification view --> triggers currentTabChanged() slot
      }
    case SPECTRA_IDX:
    default:
      break;
    }

    // update here?
  }

  inline void SpectraSelectionTabs::show(TAB_INDEX which)
  {
    setTabEnabled(which, true);
    setCurrentIndex(which);
  }
  inline SpectraIdentificationViewWidget* SpectraSelectionTabs::getSpectraIdentificationViewWidget()
  {
    return id_view_widget_;
  }
} //namespace

