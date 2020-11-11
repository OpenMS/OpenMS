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
#include <OpenMS/DATASTRUCTURES/OSWData.h>
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

  void TVDIATreeTabController::showData(const OSWIndexTrace& /*trace*/)
  {
  /*
    // basic behavior 1

    // show multiple spectra together is only used for chromatograms directly
    // where multiple (SRM) traces are shown together
    LayerData& layer = const_cast<LayerData&>(tv_->getActiveCanvas()->getCurrentLayer());
    ExperimentSharedPtrType exp_sptr = layer.getPeakDataMuteable();

    OSWData* data = layer.getChromatogramAnnotation().get();
    if (data == nullptr)
    { // no OSWData available ... strange...
      return;
    }

    //open new 1D widget
    Plot1DWidget* w = new Plot1DWidget(tv_->getSpectrumParameters(1), (QWidget*)tv_->getWorkspace());

    // string for naming the different chromatogram layers with their index
    String chromatogram_caption;
    // string for naming the tab title with the indices of the chromatograms
    String caption = layer.getName();
    switch (trace.lowest)
    {
    case OSWHierarchy::Level::PROTEIN:
      // do nothing else -- showing all transitions for a protein is overwhelming...      
      break;
    case OSWHierarchy::Level::PEPTIDE:
    {
      const auto& prot = data.getProteins()[tr.idx_prot];
      const auto& pep = prot.getPeptidePrecursors()[tr.idx_pep];
      for (const auto& feat : pep.getFeatures())
      {
        const auto& trids = feat.getTransitionIDs();
        transitions_to_show.insert(transitions_to_show.end(), trids.begin(), trids.end());
      }
      break;
    }
    case OSWHierarchy::Level::FEATURE:
    {
      const auto& prot = data.getProteins()[tr.idx_prot];
      const auto& pep = prot.getPeptidePrecursors()[tr.idx_pep];
      const auto& feat = pep.getFeatures()[tr.idx_feat];
      const auto& trids = feat.getTransitionIDs();
      transitions_to_show.insert(transitions_to_show.end(), trids.begin(), trids.end());
      break;
    }
    case OSWHierarchy::Level::TRANSITION:
    {
      const auto& prot = data.getProteins()[tr.idx_prot];
      const auto& pep = prot.getPeptidePrecursors()[tr.idx_pep];
      const auto& feat = pep.getFeatures()[tr.idx_feat];
      const auto& trid = feat.getTransitionIDs()[tr.idx_trans];
      transitions_to_show.insert(transitions_to_show.end(), trid);
      break;
    }
    default:
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    // add data and return if something went wrong
    if (!w->canvas()->addLayer(exp_sptr, od_exp_sptr, layer.filename)
      || (Size)spectrum_index >= w->canvas()->getCurrentLayer().getPeakData()->size())
    {
      return;
    }

    
    // fix legend if its a chromatogram
    w->xAxis()->setLegend(PlotWidget::RT_AXIS_TITLE);


    // set relative (%) view of visible area
    w->canvas()->setIntensityMode(PlotCanvas::IM_SNAP);

    tv_->showPlotWidgetInWindow(w, caption);
    tv_->updateBarsAndMenus();

    */
  }


} // OpenMS

