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
    LayerData & layer = const_cast<LayerData&>(tv_->getActiveCanvas()->getCurrentLayer());
    ExperimentSharedPtrType exp_sptr = layer.getPeakDataMuteable();
    LayerData::ODExperimentSharedPtrType od_exp_sptr = layer.getOnDiscPeakData();
    auto ondisc_sptr = layer.getOnDiscPeakData();

    // create new 1D widget; if we return due to error, the widget will be cleaned up automatically
    unique_ptr<Plot1DWidget> wp(new Plot1DWidget(tv_->getSpectrumParameters(1), (QWidget*)tv_->getWorkspace()));
    Plot1DWidget* w = wp.get();

    if (layer.type == LayerData::DT_CHROMATOGRAM)
    {
      // fix legend and set layer name
      String caption = layer.getName() + "[" + index + "]";
      w->xAxis()->setLegend(PlotWidget::RT_AXIS_TITLE);

      // add chromatogram data as peak spectrum
      if (!w->canvas()->addChromLayer(exp_sptr, ondisc_sptr, layer.getChromatogramAnnotation(), index, layer.filename, caption, false))
      {
        return;
      }
      w->canvas()->activateSpectrum(index);
    }
    else if (layer.type == LayerData::DT_PEAK)
    {
      String caption = layer.getName();

      // add data
      if (!w->canvas()->addLayer(exp_sptr, od_exp_sptr, layer.filename) || (Size)index >= w->canvas()->getCurrentLayer().getPeakData()->size())
      {
        return;
      }
      w->canvas()->activateSpectrum(index);
    }
    else
    {
      // Behavior if its neither (user may have clicked on an empty tree or a
      // dummy entry as drawn by SpectraTreeTab::updateEntries)
      QMessageBox::critical(w, "Error", "Cannot open data that is neither chromatogram nor spectrum data. Aborting!");
      return;
    }

    // set relative (%) view of visible area
    w->canvas()->setIntensityMode(PlotCanvas::IM_SNAP);

    if (layer.type == LayerData::DT_PEAK)
    {
      //for MS1 spectra set visible area to visible area in 2D view.
      UInt ms_level = w->canvas()->getCurrentLayer().getCurrentSpectrum().getMSLevel();
      if (ms_level == 1)
      {
        // set visible area to visible area in 2D view
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
    }

    // basic behavior 2
    String caption = layer.getName();
    w->canvas()->setLayerName(w->canvas()->getCurrentLayerIndex(), caption);

    tv_->showPlotWidgetInWindow(w, caption);
    wp.release();

    tv_->updateLayerBar();
    tv_->updateViewBar();
    tv_->updateFilterBar();
    tv_->updateMenu();
  }

  void TVSpectraViewController::showChromatogramsAsNew1D(const std::vector<int>& indices)
  {

    // basic behavior 1

    // show multiple spectra together is only used for chromatograms directly
    // where multiple (SRM) traces are shown together
    LayerData & layer = const_cast<LayerData&>(tv_->getActiveCanvas()->getCurrentLayer());
    ExperimentSharedPtrType exp_sptr = layer.getPeakDataMuteable();
    auto ondisc_sptr = layer.getOnDiscPeakData();

    // string for naming the different chromatogram layers with their index
    String chromatogram_caption;
    // string for naming the tab title with the indices of the chromatograms
    String caption = layer.getName();

    //open new 1D widget
    Plot1DWidget * w = new Plot1DWidget(tv_->getSpectrumParameters(1), (QWidget *)tv_->getWorkspace());
    // fix legend if its a chromatogram
    w->xAxis()->setLegend(PlotWidget::RT_AXIS_TITLE);

    for (const auto& index : indices)
    {
      if (layer.type == LayerData::DT_CHROMATOGRAM)
      {
        // fix legend and set layer name
        caption += String(" [") + index + "];";
        chromatogram_caption = layer.getName() + "[" + index + "]";

        // add chromatogram data as peak spectrum
        if (!w->canvas()->addChromLayer(exp_sptr, ondisc_sptr, layer.getChromatogramAnnotation(), index, layer.filename, chromatogram_caption, true))
        {
          return;
        }

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
    w->canvas()->setIntensityMode(PlotCanvas::IM_SNAP);

    // basic behavior 2

    tv_->showPlotWidgetInWindow(w, caption);
    tv_->updateBarsAndMenus();
  }

  // called by SpectraTreeTab::spectrumSelected()
  void TVSpectraViewController::activate1DSpectrum(int index)
  {
    Plot1DWidget* widget_1d = tv_->getActive1DWidget();

    // return if no active 1D widget is present or no layers are present (e.g. the addLayer call failed)
    if (widget_1d == nullptr) return;
    if (widget_1d->canvas()->getLayerCount() == 0) return;

    LayerData& layer = widget_1d->canvas()->getCurrentLayer();

    // If we have a chromatogram, we cannot just simply activate this spectrum.
    // we have to do much more work, e.g. creating a new experiment with the
    // new spectrum.
    if (!layer.chromatogram_flag_set())
    {
      widget_1d->canvas()->activateSpectrum(index);
    }
    else 
    {
      widget_1d->canvas()->blockSignals(true);
      RAIICleanup clean([&]() {widget_1d->canvas()->blockSignals(false); });

      // first get raw data (the full experiment with all chromatograms), we
      // only need to grab the one with the desired index
      ExperimentSharedPtrType exp_sptr = layer.getChromatogramData();
      auto ondisc_sptr = layer.getOnDiscPeakData();

      widget_1d->canvas()->removeLayers();

      // fix legend and set layer name
      String fname = layer.filename;
      String caption = fname + "[" + index + "]";

      // add chromatogram data as peak spectrum and update other controls
      widget_1d->canvas()->addChromLayer(exp_sptr, ondisc_sptr, layer.getChromatogramAnnotation(), index, fname, caption, false);

      tv_->updateBarsAndMenus(); // needed since we blocked update above (to avoid repeated layer updates for many layers!)
    }
  }

  // called by SpectraTreeTab::chromsSelected()
  void TVSpectraViewController::activate1DSpectrum(const std::vector<int>& indices)
  {
    Plot1DWidget * widget_1d = tv_->getActive1DWidget();

    // return if no active 1D widget is present or no layers are present (e.g. the addLayer call failed)
    if (widget_1d == nullptr) return;
    if (widget_1d->canvas()->getLayerCount() == 0) return;

    const LayerData& layer = widget_1d->canvas()->getCurrentLayer();
    // If we have a chromatogram, we cannot just simply activate this spectrum.
    // we have to do much more work, e.g. creating a new experiment with the
    // new spectrum.
    if (layer.chromatogram_flag_set())
    {
      // first get raw data (the full experiment with all chromatograms), we
      // only need to grab the ones with the desired indices
      ExperimentSharedPtrType exp_sptr = layer.getChromatogramData();
      auto ondisc_sptr = layer.getOnDiscPeakData();

      widget_1d->canvas()->removeLayers();

      widget_1d->canvas()->blockSignals(true);
      RAIICleanup clean([&]() {widget_1d->canvas()->blockSignals(false); });
      String fname = layer.filename;
      for (const auto& index : indices)
      {
        // get caption (either chromatogram idx or peptide sequence, if available)
        String caption = fname;
        if (exp_sptr->metaValueExists("peptide_sequence"))
        {
          caption = String(exp_sptr->getMetaValue("peptide_sequence"));
        }
        ((caption += "[") += index) += "]";
        // add chromatogram data as peak spectrum
        widget_1d->canvas()->addChromLayer(exp_sptr, ondisc_sptr, layer.getChromatogramAnnotation(), index, fname, caption, true);
      }

      tv_->updateBarsAndMenus(); // needed since we blocked update above (to avoid repeated layer updates for many layers!)
    }
  }

  void TVSpectraViewController::deactivate1DSpectrum(int /* spectrum_index */)
  {
    // no special handling of spectrum deactivation needed
  }

} // OpenMS

