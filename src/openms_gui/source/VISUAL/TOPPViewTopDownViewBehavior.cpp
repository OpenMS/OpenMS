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

#include <OpenMS/VISUAL/TOPPViewTopDownViewBehavior.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DCaret.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/TopDownViewWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DVerticalLineItem.h>
#include <boost/range/adaptor/reversed.hpp>
#include <QtWidgets/QMessageBox>
#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

//#define DEBUG_IDENTIFICATION_VIEW

namespace OpenMS
{
  TOPPViewTopDownViewBehavior::TOPPViewTopDownViewBehavior(TOPPViewBase* parent, TopDownViewWidget* spec_id_view) :
    tv_(parent),
    spec_id_view_(spec_id_view)
  {
  }

  void TOPPViewTopDownViewBehavior::showSpectrumAs1D(int spectrum_index)
  {
    LayerData & layer = const_cast<LayerData&>(tv_->getActiveCanvas()->getCurrentLayer());
    LayerData::ExperimentSharedPtrType exp_sptr = layer.getPeakDataMuteable();
    LayerData::ODExperimentSharedPtrType od_exp_sptr = layer.getOnDiscPeakData();

    if (layer.type == LayerData::DT_PEAK)
    {
      // open new 1D widget with the current default parameters
      Spectrum1DWidget* w = new Spectrum1DWidget(tv_->getSpectrumParameters(1), (QWidget*)tv_->getWorkspace());

      // add data and return if something went wrong
      if (!w->canvas()->addLayer(exp_sptr, od_exp_sptr, layer.filename)
        || (Size)spectrum_index >= w->canvas()->getCurrentLayer().getPeakData()->size())
      {
        return;
      }

      w->canvas()->activateSpectrum(spectrum_index);

      // set relative (%) view of visible area
      w->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);

      // for MS1 spectra set visible area to visible area in 2D view.
      UInt ms_level = w->canvas()->getCurrentLayer().getCurrentSpectrum().getMSLevel();
      if (ms_level == 1)
      {
        // set visible area to visible area in 2D view
        SpectrumCanvas::AreaType a = tv_->getActiveCanvas()->getVisibleArea();
        w->canvas()->setVisibleArea(a);
      }

      String caption = layer.getName();
      w->canvas()->setLayerName(w->canvas()->getCurrentLayerIndex(), caption);

      tv_->showSpectrumWidgetInWindow(w, caption);

      // mass annotation
      addPeakAnnotations_(vector<double>()); // TODO: values

      tv_->updateLayerBar(); // todo replace
      tv_->updateViewBar();
      tv_->updateFilterBar();
      tv_->updateMenu();
    }
    else
    {
       throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }    
  }

  void TOPPViewTopDownViewBehavior::addPeakAnnotations_(const vector<double>& masses)
  {
    // called anew for every click on a spectrum
    LayerData& current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();

    if (current_layer.getCurrentSpectrum().empty())
    {
      OPENMS_LOG_WARN << "Spectrum is empty! Nothing to annotate!" << endl;
      return;
    }

    if (!current_layer.getCurrentSpectrum().isSorted())
    {
      QMessageBox::warning(tv_, "Error", "The spectrum is not sorted! Aborting!");
      return;
    }


    for (double mass : masses)
    {
      /*
       // mass precision to match a peak's m/z to a feature m/z
       // m/z values of features are usually an average over multiple scans...
       const double ppm = 0.5;
       Size peak_idx = current_layer.getCurrentSpectrum().findNearest(mass); // TODO that
       // m/z fits ?
       if (Math::getPPMAbs(mass, current_layer.getCurrentSpectrum()[peak_idx].getMZ()) > ppm) continue; 
       // TODO: add stuff like the new VerticalAnnotationItem or Annotation1DCaret* first_dit(nullptr);
     
        Annotation1DCaret* ditem = new Annotation1DCaret(points,
                                                         QString(),
                                                         cols[i],
                                                         current_layer.param.getValue("peak_color").toQString());
        ditem->setSelected(false);
        temporary_annotations_.push_back(ditem); // for removal (no ownership)
        current_layer.getCurrentAnnotations().push_front(ditem); // for visualization (ownership)
      */

      auto ai = new Annotation1DVerticalLineItem(mass, Qt::green, QString::number(mass, 'f', 3));
      temporary_annotations_.push_back(ai); // for removal (no ownership)
      current_layer.getCurrentAnnotations().push_front(ai); // for visualization (ownership)

    }
  }

  void TOPPViewTopDownViewBehavior::activate1DSpectrum(int spectrum_index)
  {
    Spectrum1DWidget* widget_1D = tv_->getActive1DWidget();

    // return if no active 1D widget is present
    if (widget_1D == nullptr) { return; }

    widget_1D->canvas()->activateSpectrum(spectrum_index);
    LayerData& current_layer = widget_1D->canvas()->getCurrentLayer();

    if (current_layer.type == LayerData::DT_PEAK)
    {
      addPeakAnnotations_({100,200,300,400,500}); // TODO: for now just hardcode some masses
    } // end DT_PEAK
    else
    {
       throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
  }
  
  // TODO: same code for IdentificationView -> refactor
  void TOPPViewTopDownViewBehavior::addPrecursorLabels1D_(const vector<Precursor>& pcs)
  {
    LayerData& current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();

    if (current_layer.type == LayerData::DT_PEAK)
    {
      const auto& spectrum = current_layer.getCurrentSpectrum();

      for (vector<Precursor>::const_iterator it = pcs.begin(); it != pcs.end(); ++it)
      {
        // determine start and stop of isolation window
        double center_mz = it->metaValueExists("isolation window target m/z") ?
          double(it->getMetaValue("isolation window target m/z")) : it->getMZ();
        double isolation_window_lower_mz = center_mz - it->getIsolationWindowLowerOffset();
        double isolation_window_upper_mz = center_mz + it->getIsolationWindowUpperOffset();

        // determine maximum peak intensity in isolation window
        auto vbegin = spectrum.MZBegin(isolation_window_lower_mz);
        auto vend = spectrum.MZEnd(isolation_window_upper_mz);

        double max_intensity = (numeric_limits<double>::min)();
        for (; vbegin != vend; ++vbegin)
        {
          if (vbegin->getIntensity() > max_intensity)
          {
            max_intensity = vbegin->getIntensity();
          }
        }

        // DPosition<2> precursor_position = DPosition<2>(it->getMZ(), max_intensity);
        DPosition<2> lower_position = DPosition<2>(isolation_window_lower_mz, max_intensity);
        DPosition<2> upper_position = DPosition<2>(isolation_window_upper_mz, max_intensity);

        Annotation1DDistanceItem* item = new Annotation1DDistanceItem(QString::number(it->getCharge()), lower_position, upper_position);
        // add additional tick at precursor target position (e.g. to show if isolation window is asymmetric)
        vector<double> ticks;
        ticks.push_back(it->getMZ());
        item->setTicks(ticks);
        item->setSelected(false);

        temporary_annotations_.push_back(item); // for removal (no ownership)
        current_layer.getCurrentAnnotations().push_front(item); // for visualization (ownership)
      }
    }
    else
    {
       throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }
  }

  // TODO: same code for IdentificationView -> refactor and move to LayerData
  void TOPPViewTopDownViewBehavior::removeTemporaryAnnotations_(Size spectrum_index)
  {
#ifdef DEBUG_IDENTIFICATION_VIEW
    cout << "removePrecursorLabels1D_ " << spectrum_index << endl;
#endif
    // Delete annotations added by IdentificationView (but not user added annotations)
    LayerData& current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();
    const vector<Annotation1DItem*>& cas = temporary_annotations_;
    Annotations1DContainer& las = current_layer.getAnnotations(spectrum_index);
    for (vector<Annotation1DItem*>::const_iterator it = cas.begin(); it != cas.end(); ++it)
    {
      Annotations1DContainer::iterator i = find(las.begin(), las.end(), *it);
      if (i != las.end())
      {
        delete(*i);
        las.erase(i);
      }
    }
    temporary_annotations_.clear();
  }

  // TODO: same code for IdentificationView? -> refactor and move to LayerData
  void TOPPViewTopDownViewBehavior::removeGraphicalPeakAnnotations_(int spectrum_index)
  {
    Spectrum1DWidget* widget_1D = tv_->getActive1DWidget();
    LayerData& current_layer = widget_1D->canvas()->getCurrentLayer();

    #ifdef DEBUG_IDENTIFICATION_VIEW
          cout << "Removing peak annotations." << endl;
    #endif
    // remove all graphical peak annotations as these will be recreated from the stored peak annotations
    Annotations1DContainer& las = current_layer.getAnnotations(spectrum_index);
    auto new_end = remove_if(las.begin(), las.end(),
                            [](const Annotation1DItem* a)
                            {
                              #ifdef DEBUG_IDENTIFICATION_VIEW
                              cout << a->getText().toStdString() << endl;
                              #endif
                              return dynamic_cast<const Annotation1DPeakItem*>(a) != nullptr;
                            });
    las.erase(new_end, las.end());

    return;
  }

  void TOPPViewTopDownViewBehavior::deactivate1DSpectrum(int spectrum_index)
  {
    // Retrieve active 1D widget
    Spectrum1DWidget* widget_1D = tv_->getActive1DWidget();

    // Return if none present
    if (widget_1D == nullptr) return;

    LayerData& current_layer = widget_1D->canvas()->getCurrentLayer();

    // Return if no valid peak layer attached
    if (current_layer.getPeakData()->size() == 0 || current_layer.type != LayerData::DT_PEAK) { return; }

    MSSpectrum& spectrum = (*current_layer.getPeakDataMuteable())[spectrum_index];
    int ms_level = spectrum.getMSLevel();
    if (ms_level == 2)
    {
      removeGraphicalPeakAnnotations_(spectrum_index);
    }

    removeTemporaryAnnotations_(spectrum_index);

    widget_1D->canvas()->setTextBox(QString());
  }

  void TOPPViewTopDownViewBehavior::activateBehavior()
  {
    // set relative (%) view of visible area
    Spectrum1DWidget* w = tv_->getActive1DWidget();
    if (w == nullptr) return;
    w->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);
  }

  void TOPPViewTopDownViewBehavior::deactivateBehavior()
  {
    Spectrum1DWidget* widget_1D = tv_->getActive1DWidget();    
    if (widget_1D == nullptr) return; // return if no active 1D widget is present

    // remove precusor labels, theoretical spectra and trigger repaint
    LayerData& cl = tv_->getActive1DWidget()->canvas()->getCurrentLayer();
    removeTemporaryAnnotations_(cl.getCurrentSpectrumIndex());

    tv_->getActive1DWidget()->canvas()->repaint();
  }

  void TOPPViewTopDownViewBehavior::setVisibleArea1D(double l, double h)
  {
    Spectrum1DWidget* widget_1D = tv_->getActive1DWidget();

    // return if no active 1D widget is present
    if (widget_1D == nullptr) return;

    DRange<2> range = tv_->getActive1DWidget()->canvas()->getVisibleArea();
    range.setMinX(l);
    range.setMaxX(h);
    tv_->getActive1DWidget()->canvas()->setVisibleArea(range);
    tv_->getActive1DWidget()->canvas()->repaint();
  }

}
