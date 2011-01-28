// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/TOPPViewSpectraViewBehavior.h>

#include <QtGui/QMessageBox>
#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  TOPPViewSpectraViewBehavior::TOPPViewSpectraViewBehavior(TOPPViewBase* parent):
      tv_(parent)
  {
  }

  void TOPPViewSpectraViewBehavior::showSpectrumAs1D(int index)
  {    
    // basic behavior 1
    const LayerData& layer = tv_->getActiveCanvas()->getCurrentLayer();
    ExperimentSharedPtrType exp_sptr = layer.getPeakData();

   //open new 1D widget
   Spectrum1DWidget* w = new Spectrum1DWidget(tv_->getSpectrumParameters(1), (QWidget*)tv_->getWorkspace());

   //add data
   if (!w->canvas()->addLayer(exp_sptr, layer.filename) || (Size)index >= w->canvas()->getCurrentLayer().getPeakData()->size())
   {
     return;
   }

   w->canvas()->activateSpectrum(index);

   // set relative (%) view of visible area
   w->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);

   //for MS1 spectra set visible area to visible area in 2D view.
   UInt ms_level = w->canvas()->getCurrentLayer().getCurrentSpectrum().getMSLevel();
   if (ms_level == 1)
   {
     // set visible aree to visible area in 2D view
     w->canvas()->setVisibleArea(tv_->getActiveCanvas()->getVisibleArea());
   }

   // basic behavior 2
   String caption = layer.name;
   w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);

   tv_->showSpectrumWidgetInWindow(w,caption);
   tv_->updateLayerBar();
   tv_->updateViewBar();
   tv_->updateFilterBar();
   tv_->updateMenu();
  }

  void TOPPViewSpectraViewBehavior::activate1DSpectrum(int index)
  {
    Spectrum1DWidget* widget_1d = tv_->getActive1DWidget();
    widget_1d->canvas()->activateSpectrum(index);
  }

  void TOPPViewSpectraViewBehavior::deactivate1DSpectrum(int /* spectrum_index */)
  {
    // no special handling of spectrum deactivation needed
  }

  void TOPPViewSpectraViewBehavior::activateBehavior()
  {
    // no special handling of activation
  }

  void TOPPViewSpectraViewBehavior::deactivateBehavior()
  {
    // no special handling of deactivation
  }

} // OpenMS
