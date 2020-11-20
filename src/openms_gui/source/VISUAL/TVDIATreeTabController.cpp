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
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DVerticalLineItem.h>


#include <QtWidgets/QMessageBox>
#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{


  typedef LayerData::ExperimentSharedPtrType ExperimentSharedPtrType;
  typedef LayerData::ConstExperimentSharedPtrType ConstExperimentSharedPtrType;
  typedef LayerData::ODExperimentSharedPtrType ODExperimentSharedPtrType;
  typedef LayerData::OSWDataSharedPtrType OSWDataSharedPtrType;

  /// represents all the information we need from a layer
  /// We cannot use a full layer, because the original layer might get destroyed in the process...
  struct MiniLayer
  {
    ExperimentSharedPtrType full_chrom_exp_sptr;
    ODExperimentSharedPtrType ondisc_sptr;
    OSWDataSharedPtrType annot_sptr;
    String filename;
    String layername;

    MiniLayer(LayerData& layer)
    : full_chrom_exp_sptr(getFullChromData(layer)),
      ondisc_sptr(layer.getOnDiscPeakData()),
      annot_sptr(layer.getChromatogramAnnotation()),
      filename(layer.filename),
      layername(layer.getName())
    {
    }

    /// get the full chromExperiment
    /// Could be backed up in layer.getChromatogramData() (if layer.getPeakDataMuteable() shows converted chroms already; as we will do below)
    /// ... or layer.getChromatogramData() is empty and thus layer.getPeakDataMuteable() is the original chrom data
    static ExperimentSharedPtrType getFullChromData(LayerData& layer)
    {
      ExperimentSharedPtrType exp_sptr(layer.getChromatogramData().get() == nullptr ||
                                       layer.getChromatogramData().get()->getNrChromatograms() == 0
                                       ? layer.getPeakDataMuteable() : layer.getChromatogramData());
      return exp_sptr;
    }
  };

  

  TVDIATreeTabController::TVDIATreeTabController(TOPPViewBase* parent):
    TVControllerBase(parent)
  {  
  }

  bool addTransitionAsLayer(Plot1DWidget* w, 
                            MiniLayer ml,
                            const int chrom_index,
                            const OSWPeakGroup* feature = nullptr)
  {
    String chrom_caption = FileHandler::stripExtension(File::basename(ml.filename)) + "[" + chrom_index + "]";

    // add data and return if something went wrong
    if (!w->canvas()->addChromLayer(ml.full_chrom_exp_sptr, ml.ondisc_sptr, ml.annot_sptr, chrom_index, ml.filename, chrom_caption, false))
    {
      return false;
    }
    // add boundaries
    if (feature)
    {
      double width = feature->getRTRightWidth() - feature->getRTLeftWidth();
      double center = feature->getRTLeftWidth() + width / 2;
      String ann = String("RT: ") + String(feature->getRTExperimental(), false) + "\ndRT: " + String(feature->getRTDelta(), false) + "\nQ-Value: " + String(feature->getQValue(), false);
      Annotation1DItem* item = new Annotation1DVerticalLineItem(center, width, 30, QColor("invalid"), ann.toQString());
      item->setSelected(false);
      w->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
    }


    w->canvas()->activateSpectrum(0, false);
    return true;
  }


  void TVDIATreeTabController::showChromatogramsAsNew1D(const OSWIndexTrace& trace)
  {
    LayerData& layer = const_cast<LayerData&>(tv_->getActiveCanvas()->getCurrentLayer());
    MiniLayer ml(layer);
    // create new 1D widget; if we return due to error, the widget will be cleaned up
    unique_ptr<Plot1DWidget> w(new Plot1DWidget(tv_->getSpectrumParameters(1), (QWidget*)tv_->getWorkspace()));

    if (showChromatogramsInCanvas_(trace, ml, w.get()))
    { // success!
      tv_->showPlotWidgetInWindow(w.get(), ml.layername);
      w.release(); // do NOT delete the widget; tv_ owns it now ...
      tv_->updateBarsAndMenus();
    }
  }

  void TVDIATreeTabController::showChromatograms(const OSWIndexTrace& trace)
  {
    Plot1DWidget* w = tv_->getActive1DWidget();
    if (w == nullptr)
    { // currently not a 1d widget... ignore the signal
      return;
    }
    MiniLayer ml(w->canvas()->getCurrentLayer());
    // clear all layers
    w->canvas()->removeLayers();
    // add new layers
    if (showChromatogramsInCanvas_(trace, ml, w))
    {
      tv_->updateBarsAndMenus();
    }
  }

  bool TVDIATreeTabController::showChromatogramsInCanvas_(const OSWIndexTrace& trace, MiniLayer& ml, Plot1DWidget* w)
  {

    OSWData* data = ml.annot_sptr.get();
    if (data == nullptr)
    { // no OSWData available ... strange...
      return false;
    }

    switch (trace.lowest)
    {
    case OSWHierarchy::Level::PROTEIN:
    {
      const auto& prot = data->getProteins()[trace.idx_prot];
      // show only the first peptide for now...
      const auto& pep = prot.getPeptidePrecursors()[0];
      for (const auto& feat : pep.getFeatures())
      {
        const auto& trids = feat.getTransitionIDs();
        for (UInt trid : trids)
        {
          if (!addTransitionAsLayer(w, ml, (Size)trid, &feat))
          { // something went wrong. abort
            return false;
          }
        }
      }
      break;
    }
    case OSWHierarchy::Level::PEPTIDE:
    {
      const auto& prot = data->getProteins()[trace.idx_prot];
      const auto& pep = prot.getPeptidePrecursors()[trace.idx_pep];
      for (const auto& feat : pep.getFeatures())
      {
        const auto& trids = feat.getTransitionIDs();
        for (UInt trid : trids)
        {
          if (!addTransitionAsLayer(w, ml, (Size)trid, &feat))
          { // something went wrong. abort
            return false;
          }
        }
      }
      break;
    }
    case OSWHierarchy::Level::FEATURE:
    {
      const auto& prot = data->getProteins()[trace.idx_prot];
      const auto& pep = prot.getPeptidePrecursors()[trace.idx_pep];
      const auto& feat = pep.getFeatures()[trace.idx_feat];
      const auto& trids = feat.getTransitionIDs();
      for (UInt trid : trids)
      {
        if (!addTransitionAsLayer(w, ml, (Size)trid, &feat))
        { // something went wrong. abort
          return false;
        }
      }
      break;
    }
    case OSWHierarchy::Level::TRANSITION:
    {
      const auto& prot = data->getProteins()[trace.idx_prot];
      const auto& pep = prot.getPeptidePrecursors()[trace.idx_pep];
      const auto& feat = pep.getFeatures()[trace.idx_feat];
      const auto& trid = feat.getTransitionIDs()[trace.idx_trans];
      if (!addTransitionAsLayer(w, ml, (Size)trid, &feat))
      { // something went wrong. abort
        return false;
      }
      break;
    }
    default:
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    return true;
  }


} // OpenMS

