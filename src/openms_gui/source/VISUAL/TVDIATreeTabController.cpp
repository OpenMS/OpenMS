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

#include <OpenMS/VISUAL/TVDIATreeTabController.h>

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
  TVDIATreeTabController::TVDIATreeTabController(TOPPViewBase* parent):
    TVControllerBase(parent)
  {  
  }

  void TVDIATreeTabController::showTransitionsAs1D(const std::vector<int>& /*indices*/)
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


    // set relative (%) view of visible area
    w->canvas()->setIntensityMode(PlotCanvas::IM_SNAP);

    // basic behavior 2

    tv_->showPlotWidgetInWindow(w, caption);
    tv_->updateBarsAndMenus();
  }

  // called by SpectraTreeTab::chromsSelected()
  void TVDIATreeTabController::activate1DTransitions(const std::vector<int>& /*indices*/)
  {
    Plot1DWidget * widget_1d = tv_->getActive1DWidget();

    // return if no active 1D widget is present or no layers are present (e.g. the addLayer call failed)
    if (widget_1d == nullptr) return;
    if (widget_1d->canvas()->getLayerCount() == 0) return;

    const LayerData& layer = widget_1d->canvas()->getCurrentLayer();
    if (layer.getChromatogramAnnotation().get())
    {
      widget_1d->canvas()->removeLayers();
      widget_1d->canvas()->blockSignals(true);
      RAIICleanup clean([&]() {widget_1d->canvas()->blockSignals(false); });
      /*
      String fname = layer.filename;
      for (const auto& index : indices)
      {
        ExperimentSharedPtrType chrom_exp_sptr = prepareChromatogram(index, exp_sptr, ondisc_sptr);

        // get caption (either chromatogram idx or peptide sequence, if available)
        String caption = fname;
        if (chrom_exp_sptr->metaValueExists("peptide_sequence"))
        {
          caption = String(chrom_exp_sptr->getMetaValue("peptide_sequence"));
        }
        ((caption += "[") += index) += "]";
        // add chromatogram data as peak spectrum
        widget_1d->canvas()->addChromLayer(chrom_exp_sptr, ondisc_sptr, fname, caption, exp_sptr, index, true);
      }
           */
      tv_->updateBarsAndMenus(); // needed since we blocked update above (to avoid repeated layer updates for many layers!)
    }
  }

  void TVDIATreeTabController::deactivate1DTransitions(const std::vector<int>& /*indices*/)
  {
    // no special handling of spectrum deactivation needed
  }

} // OpenMS

