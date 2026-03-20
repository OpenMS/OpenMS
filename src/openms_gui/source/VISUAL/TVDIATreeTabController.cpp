// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/LayerDataChrom.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DVerticalLineItem.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{


  typedef LayerDataBase::ExperimentSharedPtrType ExperimentSharedPtrType;
  typedef LayerDataBase::ConstExperimentSharedPtrType ConstExperimentSharedPtrType;
  typedef LayerDataBase::ODExperimentSharedPtrType ODExperimentSharedPtrType;
  typedef LayerDataBase::OSWDataSharedPtrType OSWDataSharedPtrType;

  /// represents all the information we need from a chromatogram layer
  /// We cannot use a full layer, because the original layer might get destroyed in the process...
  struct MiniLayer
  {
    ExperimentSharedPtrType full_chrom_exp_sptr;
    ODExperimentSharedPtrType ondisc_sptr;
    OSWDataSharedPtrType annot_sptr;
    String filename;
    String layername;

    explicit MiniLayer(LayerDataChrom& layer)
    : full_chrom_exp_sptr(layer.getChromatogramData()),
      ondisc_sptr(layer.getOnDiscPeakData()),
      annot_sptr(layer.getChromatogramAnnotation()),
      filename(layer.filename),
      layername(layer.getName())
    {
    }
  };


  bool addTransitionAsLayer(Plot1DWidget* w, 
                            const MiniLayer& ml,
                            const int transition_id,
                            std::set<UInt32>& transitions_seen)
  {
    if (transitions_seen.find(transition_id) != transitions_seen.end())
    { // duplicate .. do not show
      return true;
    }
    transitions_seen.insert(transition_id);

    // convert from native id to chrom_index
    int chrom_index = ml.annot_sptr->fromNativeID(transition_id);

    // add data and return if something went wrong
    if (!w->canvas()->addChromLayer(ml.full_chrom_exp_sptr, ml.ondisc_sptr, ml.annot_sptr,
                                    chrom_index, ml.filename, 
                                    FileHandler::stripExtension(File::basename(ml.filename)),
                                    String("[") + transition_id + "]"))
    {
      return false;
    }
    w->canvas()->activateSpectrum(chrom_index, false);
    return true;
  }

  void addFeatures(Plot1DWidget* w, std::vector<OSWPeakGroup>& features)
  {
    // nothing to do...
    if (features.empty())
    {
      return;
    }
    // sort features by left RT
    std::sort(features.begin(), features.end(), [](const OSWPeakGroup& a, const OSWPeakGroup& b)
    {
      return a.getRTLeftWidth() < b.getRTLeftWidth();
    });
    const OSWPeakGroup* best_feature = &features[0];
    auto findBestFeature = [&best_feature](const OSWPeakGroup& f)
    {
      if (best_feature->getQValue() > f.getQValue())
      {
        best_feature = &f;
      }
    };
    std::for_each(features.begin(), features.end(), findBestFeature);
    if (best_feature->getQValue() == -1)
    { // no q-values are annotated. make them all grey.
      best_feature = nullptr;
    }
    GUIHelpers::OverlapDetector od(3); // three y-levels for showing annotation


    // show feature boundaries
    for (const auto& feature : features)
    {
      auto width = feature.getRTRightWidth() - feature.getRTLeftWidth();
      auto center = feature.getRTLeftWidth() + width / 2;
      String ann = String("RT:\n ") + String(feature.getRTExperimental(), false) + "\ndRT:\n " + String(feature.getRTDelta(), false) + "\nQ:\n " + String(feature.getQValue(), false);
      QColor col = GUIHelpers::ColorBrewer::Distinct().values[(best_feature == &feature) 
                          ? GUIHelpers::ColorBrewer::Distinct::LightGreen
                          : GUIHelpers::ColorBrewer::Distinct::LightGrey];
      Annotation1DVerticalLineItem* item = new Annotation1DVerticalLineItem(center, width, 150, false, col, ann.toQString());
      item->setSelected(false);
      auto text_size = item->getTextRect(); // this is in px units (Qt widget coordinates)
      // translate to axis units (our native 'data'):
      auto p_text = w->canvas()->widgetToDataDistance(text_size.width(), 0);
      auto chunk = od.placeItem(feature.getRTLeftWidth(), feature.getRTLeftWidth() + p_text.getX());
      item->setTextOffset(chunk * text_size.height());

      w->canvas()->getCurrentLayer().getCurrentAnnotations().push_back(item);
    }
    // paint the expected RT once
    auto expected_RT = features[0].getRTExperimental() - features[0].getRTDelta();
    Annotation1DItem* item = new Annotation1DVerticalLineItem(expected_RT, 3, 200, true, Qt::darkGreen, "");
    item->setSelected(false);
    w->canvas()->getCurrentLayer().getCurrentAnnotations().push_back(item);
  }



  TVDIATreeTabController::TVDIATreeTabController(TOPPViewBase* parent) :
    TVControllerBase(parent)
  {
  }

  void TVDIATreeTabController::showChromatogramsAsNew1D(const OSWIndexTrace& trace)
  {
    auto* layer_ptr = dynamic_cast<LayerDataChrom*>(&tv_->getActiveCanvas()->getCurrentLayer());
    if (!layer_ptr)
    { // not a chrom layer?
      std::cerr << __FILE__ << ": " << __LINE__ << " showChromatograms() invoked on Non-Chrom layer... weird..\n";
      return;
    }
    MiniLayer ml(*layer_ptr);
    // create new 1D widget; if we return due to error, the widget will be cleaned up
    unique_ptr<Plot1DWidget> w(new Plot1DWidget(tv_->getCanvasParameters(1), DIM::Y, (QWidget*)tv_->getWorkspace()));

    if (showChromatogramsInCanvas_(trace, ml, w.get()))
    { // success!
      tv_->showPlotWidgetInWindow(w.get());
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
    auto* layer_ptr = dynamic_cast<LayerDataChrom*>(&w->canvas()->getCurrentLayer());
    if (!layer_ptr)
    { // not a chrom layer?
      std::cerr << __FILE__ << ": " << __LINE__ << " showChromatograms() invoked on Non-Chrom layer... weird..\n";
      return;
    }
    MiniLayer ml(*layer_ptr);
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

    std::set<UInt32> transitions_seen;
    std::vector<OSWPeakGroup> features;

    switch (trace.lowest)
    {
    case OSWHierarchy::Level::PROTEIN:
    {
      const auto& prot = data->getProteins()[trace.idx_prot];
      // show only the first peptide for now...
      const auto& pep = prot.getPeptidePrecursors()[0];
      features = pep.getFeatures();
      for (const auto& feat : pep.getFeatures())
      {
        const auto& trids = feat.getTransitionIDs();
        for (UInt trid : trids)
        {
          if (!addTransitionAsLayer(w, ml, (Size)trid, transitions_seen))
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
      features = pep.getFeatures();
      for (const auto& feat : pep.getFeatures())
      {
        const auto& trids = feat.getTransitionIDs();
        for (UInt trid : trids)
        {
          if (!addTransitionAsLayer(w, ml, (Size)trid, transitions_seen))
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
      features = { feat };
      const auto& trids = feat.getTransitionIDs();
      for (UInt trid : trids)
      {
        if (!addTransitionAsLayer(w, ml, (Size)trid, transitions_seen))
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
      if (!addTransitionAsLayer(w, ml, (Size)trid, transitions_seen))
      { // something went wrong. abort
        return false;
      }
      break;
    }
    default:
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    // add bars for all identified features
    addFeatures(w, features);

    return true;
  }


} // OpenMS

