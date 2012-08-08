// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/TOPPViewSpectraViewBehavior.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>

#include <QtGui/QMessageBox>
#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  TOPPViewSpectraViewBehavior::TOPPViewSpectraViewBehavior(TOPPViewBase * parent) :
    tv_(parent)
  {
  }

  String caption;

  void TOPPViewSpectraViewBehavior::showSpectrumAs1D(int index)
  {

    // basic behavior 1
    const LayerData & layer = tv_->getActiveCanvas()->getCurrentLayer();
    ExperimentSharedPtrType exp_sptr = layer.getPeakData();

    // open new 1D widget
    Spectrum1DWidget * w = new Spectrum1DWidget(tv_->getSpectrumParameters(1), (QWidget *)tv_->getWorkspace());

    if (layer.type == LayerData::DT_CHROMATOGRAM)
    {
      // create managed pointer to experiment data
      ExperimentType * chrom_exp = new ExperimentType();
      ExperimentSharedPtrType chrom_exp_sptr(chrom_exp);
      chrom_exp->setMetaValue("is_chromatogram", "true"); //this is a hack to store that we have chromatogram data

      const ExperimentType exp = *exp_sptr;
      SpectrumType spectrum;
      const MSChromatogram<ChromatogramPeak> & current_chrom = exp.getChromatograms()[index];
      for (Size i = 0; i != current_chrom.size(); ++i)
      {
        const ChromatogramPeak & cpeak = current_chrom[i];
        Peak1D peak1d;
        peak1d.setMZ(cpeak.getRT());
        peak1d.setIntensity(cpeak.getIntensity());
        spectrum.push_back(peak1d);
      }
      chrom_exp_sptr->push_back(spectrum);

      caption = layer.name + "[" + index + "]";
      //add chromatogram data as peak spectrum
      if (!w->canvas()->addLayer(chrom_exp_sptr, layer.filename))
      {
        return;
      }
      w->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);

    }
    else if (layer.type == LayerData::DT_PEAK)
    {
      caption = layer.name;

      //add data
      if (!w->canvas()->addLayer(exp_sptr, layer.filename) || (Size)index >= w->canvas()->getCurrentLayer().getPeakData()->size())
      {
        return;
      }
    }

    w->canvas()->activateSpectrum(index);

    // set relative (%) view of visible area
    w->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);

    if (layer.type == LayerData::DT_PEAK)
    {
      //for MS1 spectra set visible area to visible area in 2D view.
      UInt ms_level = w->canvas()->getCurrentLayer().getCurrentSpectrum().getMSLevel();
      if (ms_level == 1)
      {
        // set visible aree to visible area in 2D view
        w->canvas()->setVisibleArea(tv_->getActiveCanvas()->getVisibleArea());
      }
    }
    else if (layer.type == LayerData::DT_CHROMATOGRAM)
    {
      // set visible area to visible area in 2D view
      // switch X/Y because now we want to have RT on the x-axis and not m/z
      DRange<2> visible_area = tv_->getActiveCanvas()->getVisibleArea();
      int tmp_x1 = visible_area.minX();
      int tmp_x2 = visible_area.maxX();
      visible_area.setMinX(visible_area.minY());
      visible_area.setMaxX(visible_area.maxY());
      visible_area.setMinY(tmp_x1);
      visible_area.setMaxY(tmp_x2);
      w->canvas()->setVisibleArea(visible_area);

      w->canvas()->getCurrentLayer().getChromatogramData() = exp_sptr; // save the original chromatogram data so that we can access it later
    }

    // basic behavior 2
    String caption = layer.name;
    w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);

    tv_->showSpectrumWidgetInWindow(w, caption);
    tv_->updateLayerBar();
    tv_->updateViewBar();
    tv_->updateFilterBar();
    tv_->updateMenu();
  }

  void TOPPViewSpectraViewBehavior::showSpectrumAs1D(std::vector<int, std::allocator<int> > indices)
  {

    // basic behavior 1
    const LayerData & layer = tv_->getActiveCanvas()->getCurrentLayer();
    ExperimentSharedPtrType exp_sptr = layer.getPeakData();

    // string for naming the different chromatogram layers with their index
    String chromatogram_caption;
    // string for naming the tab title with the indices of the chromatograms
    caption = layer.name;

    //open new 1D widget
    Spectrum1DWidget * w = new Spectrum1DWidget(tv_->getSpectrumParameters(1), (QWidget *)tv_->getWorkspace());

    for (Size index = 0; index != indices.size(); ++index)
    {
      // create managed pointer to experiment data
      ExperimentType * chrom_exp = new ExperimentType();
      ExperimentSharedPtrType chrom_exp_sptr(chrom_exp);
      chrom_exp->setMetaValue("is_chromatogram", "true"); //this is a hack to store that we have chromatogram data

      const ExperimentType exp = *exp_sptr;
      SpectrumType spectrum;
      if (layer.type == LayerData::DT_CHROMATOGRAM)
      {

        const MSChromatogram<ChromatogramPeak> & current_chrom = exp.getChromatograms()[indices[index]];
        for (Size i = 0; i != current_chrom.size(); ++i)
        {
          const ChromatogramPeak & cpeak = current_chrom[i];
          Peak1D peak1d;
          peak1d.setMZ(cpeak.getRT());
          peak1d.setIntensity(cpeak.getIntensity());
          spectrum.push_back(peak1d);
        }
        chrom_exp_sptr->push_back(spectrum);

        caption = caption + " [" + indices[index] + "];";
        chromatogram_caption = layer.name + "[" + indices[index] + "]";
        //add chromatogram data as peak spectrum
        if (!w->canvas()->addLayer(chrom_exp_sptr, layer.filename))
        {
          return;
        }
        w->canvas()->activateSpectrum(*indices.begin());
        w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), chromatogram_caption);
        w->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);

        w->canvas()->getCurrentLayer().getChromatogramData() = exp_sptr; // save the original chromatogram data so that we can access it later

        // set visible area to visible area in 2D view
        // switch X/Y because now we want to have RT on the x-axis and not m/z
        DRange<2> visible_area = tv_->getActiveCanvas()->getVisibleArea();
        int tmp_x1 = visible_area.minX();
        int tmp_x2 = visible_area.maxX();
        visible_area.setMinX(visible_area.minY());
        visible_area.setMaxX(visible_area.maxY());
        visible_area.setMinY(tmp_x1);
        visible_area.setMaxY(tmp_x2);
        w->canvas()->setVisibleArea(visible_area);
      }
    }


    // set relative (%) view of visible area
    w->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);

    // basic behavior 2

    tv_->showSpectrumWidgetInWindow(w, caption);
    tv_->updateLayerBar();
    tv_->updateViewBar();
    tv_->updateFilterBar();
    tv_->updateMenu();
  }

  void TOPPViewSpectraViewBehavior::activate1DSpectrum(int index)
  {
    Spectrum1DWidget * widget_1d = tv_->getActive1DWidget();
    widget_1d->canvas()->activateSpectrum(index);
    const LayerData & layer = tv_->getActiveCanvas()->getCurrentLayer();

    // If we have a chromatogram, we cannot just simply activate this spectrum.
    // we have to do much more work, e.g. creating a new experiment with the
    // new spectrum.
    if (layer.chromatogram_flag_set())
    {
      // first get chromatogram data
      ExperimentSharedPtrType exp_sptr = widget_1d->canvas()->getCurrentLayer().getChromatogramData();
      const LayerData & layer = widget_1d->canvas()->getCurrentLayer();
      String fname = layer.filename;
      String lname = layer.name;

      Size layercount = widget_1d->canvas()->getLayerCount();
      for (Size i = 0; i != layercount; ++i)
      {
        widget_1d->canvas()->removeLayer(0); // remove layer 0 until there are no more layers
      }

      ExperimentType * chrom_exp = new ExperimentType();
      ExperimentSharedPtrType chrom_exp_sptr(chrom_exp);

      const ExperimentType exp = *exp_sptr;
      SpectrumType spectrum;
      const MSChromatogram<ChromatogramPeak> & current_chrom = exp.getChromatograms()[index];
      for (Size i = 0; i != current_chrom.size(); ++i)
      {
        const ChromatogramPeak & cpeak = current_chrom[i];
        Peak1D peak1d;
        peak1d.setMZ(cpeak.getRT());
        peak1d.setIntensity(cpeak.getIntensity());
        spectrum.push_back(peak1d);
      }
      chrom_exp_sptr->push_back(spectrum);
      caption = lname + "[" + index + "]";
      //add chromatogram data as peak spectrum
      if (!widget_1d->canvas()->addLayer(chrom_exp_sptr, fname))
      {
        return;
      }
      widget_1d->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      widget_1d->canvas()->activateSpectrum(index);

      widget_1d->canvas()->getCurrentLayer().name = caption;
      widget_1d->canvas()->getCurrentLayer().filename = fname;
      widget_1d->canvas()->getCurrentLayer().getChromatogramData() = exp_sptr; // save the original chromatogram data so that we can access it later
      //this is a hack to store that we have chromatogram data, that we selected multiple ones and which one we selected
      widget_1d->canvas()->getCurrentLayer().getPeakData()->setMetaValue("is_chromatogram", "true");
      widget_1d->canvas()->getCurrentLayer().getPeakData()->setMetaValue("multiple_select", "false");
      widget_1d->canvas()->getCurrentLayer().getPeakData()->setMetaValue("selected_chromatogram", index);

      tv_->updateLayerBar();
      tv_->updateViewBar();
      tv_->updateFilterBar();
      tv_->updateMenu();
    }
  }

  void TOPPViewSpectraViewBehavior::activate1DSpectrum(std::vector<int, std::allocator<int> > indices)
  {
    Spectrum1DWidget * widget_1d = tv_->getActive1DWidget();

    // If we have a chromatogram, we cannot just simply activate this spectrum.
    // we have to do much more work, e.g. creating a new experiment with the
    // new spectrum.
    const LayerData & layer = widget_1d->canvas()->getCurrentLayer();
    String fname = layer.filename;
    if (layer.chromatogram_flag_set())
    {

      // first get chromatogram data
      ExperimentSharedPtrType exp_sptr = widget_1d->canvas()->getCurrentLayer().getChromatogramData();

      Size layercount = widget_1d->canvas()->getLayerCount();
      for (Size i = 0; i != layercount; ++i)
      {
        widget_1d->canvas()->removeLayer(0); // remove layer 0 until there are no more layers
      }

      for (Size index = 0; index != indices.size(); ++index)
      {

        ExperimentType * chrom_exp = new ExperimentType();
        ExperimentSharedPtrType chrom_exp_sptr(chrom_exp);
        SpectrumType spectrum;
        const MSChromatogram<ChromatogramPeak> & current_chrom = exp_sptr->getChromatograms()[indices[index]];

        for (Size i = 0; i != current_chrom.size(); ++i)
        {
          const ChromatogramPeak & cpeak = current_chrom[i];
          Peak1D peak1d;
          peak1d.setMZ(cpeak.getRT());
          peak1d.setIntensity(cpeak.getIntensity());
          spectrum.push_back(peak1d);
        }
        chrom_exp_sptr->push_back(spectrum);
        caption = fname + "[" + indices[index] + "]";
        if (current_chrom.getPrecursor().metaValueExists("peptide_sequence"))
        {
          caption = String(current_chrom.getPrecursor().getMetaValue("peptide_sequence")) + "[" + indices[index] + "]";
        }
        //add chromatogram data as peak spectrum
        if (!widget_1d->canvas()->addLayer(chrom_exp_sptr, fname))
        {
          return;
        }
        widget_1d->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
        widget_1d->canvas()->setIntensityMode(Spectrum1DCanvas::IM_NONE);
        widget_1d->canvas()->activateSpectrum(indices[index]);
        //widget_1d->canvas()->activateSpectrum(0);

        widget_1d->canvas()->getCurrentLayer().name = caption;
        widget_1d->canvas()->getCurrentLayer().filename = fname;
        widget_1d->canvas()->getCurrentLayer().getChromatogramData() = exp_sptr; // save the original chromatogram data so that we can access it later
        //this is a hack to store that we have chromatogram data, that we selected multiple ones and which one we selected
        widget_1d->canvas()->getCurrentLayer().getPeakData()->setMetaValue("is_chromatogram", "true");
        widget_1d->canvas()->getCurrentLayer().getPeakData()->setMetaValue("multiple_select", "true");
        widget_1d->canvas()->getCurrentLayer().getPeakData()->setMetaValue("selected_chromatogram", indices[index]);
      }

      tv_->updateLayerBar();
      tv_->updateViewBar();
      tv_->updateFilterBar();
      tv_->updateMenu();
    }
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
