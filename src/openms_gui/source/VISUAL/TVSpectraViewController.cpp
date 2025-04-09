// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TVSpectraViewController.h>

#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>

#include <QtWidgets/QMessageBox>
#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  TVSpectraViewController::TVSpectraViewController(TOPPViewBase* parent):
    TVControllerBase(parent)
  {  
  }

  void TVSpectraViewController::showSpectrumAsNew1D(int index)
  {
    // basic behavior 1
    LayerDataBase& layer = tv_->getActiveCanvas()->getCurrentLayer();

    // create new 1D widget; if we return due to error, the widget will be cleaned up automatically
    unique_ptr<Plot1DWidget> wp(new Plot1DWidget(tv_->getCanvasParameters(1), DIM::Y, (QWidget*)tv_->getWorkspace()));
    Plot1DWidget* w = wp.get();
    
    // copy data from current layer (keeps the TYPE and underlying data identical)
    if (!w->canvas()->addLayer(layer.to1DLayer()))
    {
      // Behavior if its neither (user may have clicked on an empty tree or a
      // dummy entry as drawn by SpectraTreeTab::updateEntries)
      QMessageBox::critical(w, "Error", "Cannot open data that is neither chromatogram nor spectrum data. Aborting!");
      return;
    }

    w->canvas()->activateSpectrum(index);

    // set visible area to visible area in 2D view
    w->canvas()->setVisibleArea(tv_->getActiveCanvas()->getVisibleArea());

    // set relative (%) view of visible area
    w->canvas()->setIntensityMode(PlotCanvas::IM_SNAP);

    tv_->showPlotWidgetInWindow(wp.release());
    tv_->updateLayerBar();
    tv_->updateViewBar();
    tv_->updateFilterBar();
    tv_->updateMenu();
  }

  bool add1DChromLayers(const std::vector<int>& indices,
                        Plot1DWidget* target, 
                        const LayerDataDefs::ExperimentSharedPtrType& chrom_exp_sptr,
                        const LayerDataDefs::ODExperimentSharedPtrType& ondisc_sptr,
                        const OSWDataSharedPtrType& chrom_annotation,
                        const String& layer_basename,
                        const String& filename)
  {
    //
    for (const auto& index : indices)
    {
      // get caption (either chromatogram idx or peptide sequence, if available)
      String basename_suffix;
      if (chrom_exp_sptr->metaValueExists("peptide_sequence"))
      {
        basename_suffix = String(chrom_exp_sptr->getMetaValue("peptide_sequence"));
      }
      ((basename_suffix += "[") += index) += "]";

      // add chromatogram data
      if (!target->canvas()->addChromLayer(chrom_exp_sptr, ondisc_sptr, chrom_annotation, index, filename, layer_basename, basename_suffix))
      {
        return false;
      }
    }
    return true;
  }

  void TVSpectraViewController::showChromatogramsAsNew1D(const std::vector<int>& indices)
  {
    // show multiple spectra together is only used for chromatograms directly
    // where multiple (SRM) traces are shown together
    auto layer_chrom = dynamic_cast<LayerDataChrom*>(&tv_->getActiveCanvas()->getCurrentLayer());
    if (!layer_chrom) return;

    auto exp_sptr = layer_chrom->getChromatogramData();
    auto ondisc_sptr = layer_chrom->getOnDiscPeakData();

    // open new 1D widget
    auto* w = new Plot1DWidget(tv_->getCanvasParameters(1), DIM::Y, (QWidget *)tv_->getWorkspace());
    // use RT + intensity mapping
    w->setMapper({{DIM_UNIT::RT, DIM_UNIT::INT}});

    if (!add1DChromLayers(indices, w, layer_chrom->getChromatogramData(), layer_chrom->getOnDiscPeakData(),
                     layer_chrom->getChromatogramAnnotation(), layer_chrom->getName(), layer_chrom->filename))
    {
      return;
    }
    // set relative (%) view of visible area (recalcs snap factor)
    w->canvas()->setIntensityMode(PlotCanvas::IM_SNAP);

    tv_->showPlotWidgetInWindow(w);
    tv_->updateBarsAndMenus();
  }

  // called by SpectraTreeTab::spectrumSelected()
  void TVSpectraViewController::activate1DSpectrum(int index)
  {
    Plot1DWidget* widget_1d = tv_->getActive1DWidget();

    // return if no active 1D widget is present or no layers are present (e.g. the addPeakLayer call failed)
    if (widget_1d == nullptr) return;
    if (widget_1d->canvas()->getLayerCount() == 0) return;

    widget_1d->canvas()->activateSpectrum(index);
  }

  // called by SpectraTreeTab::chromsSelected()
  void TVSpectraViewController::activate1DSpectrum(const std::vector<int>& indices)
  {
    Plot1DWidget * widget_1d = tv_->getActive1DWidget();

    // return if no active 1D widget is present or no layers are present (e.g. the addPeakLayer call failed)
    if (widget_1d == nullptr) return;
    if (widget_1d->canvas()->getLayerCount() == 0) return;

    const auto* layer = dynamic_cast<LayerDataChrom*>(&widget_1d->canvas()->getCurrentLayer());
    if (!layer) return;

    auto chrom_sptr = layer->getChromatogramData();
    auto ondisc_sptr = layer->getOnDiscPeakData();
    auto annotation = layer->getChromatogramAnnotation();
    const String basename = layer->getName();
    const String filename = layer->filename;
    widget_1d->canvas()->removeLayers(); // this actually deletes layer
    layer = nullptr;                     // ... make sure its not used any more

    widget_1d->canvas()->blockSignals(true);
    RAIICleanup clean([&]() {widget_1d->canvas()->blockSignals(false); });
    
    if (!add1DChromLayers(indices, widget_1d, chrom_sptr, ondisc_sptr, annotation, basename,
                          filename))
    {
      return;
    }
    // set relative (%) view of visible area (recalcs snap factor)
    widget_1d->canvas()->setIntensityMode(PlotCanvas::IM_SNAP);

    tv_->updateBarsAndMenus(); // needed since we blocked update above (to avoid repeated layer updates for many layers!)
  }

  void TVSpectraViewController::deactivate1DSpectrum(int /* spectrum_index */)
  {
    // no special handling of spectrum deactivation needed
  }

} // OpenMS

