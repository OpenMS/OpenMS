// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/TOPPViewSpectraViewBehavior.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>


// preprocessing and filtering
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>

#include <QtWidgets/QMessageBox>
#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

  LayerData::ExperimentSharedPtrType prepareChromatogram(Size index, LayerData::ExperimentSharedPtrType exp_sptr, LayerData::ODExperimentSharedPtrType ondisc_sptr)
  {
    // create a managed pointer fill it with a spectrum containing the chromatographic data
    LayerData::ExperimentSharedPtrType chrom_exp_sptr(new LayerData::ExperimentType());
    chrom_exp_sptr->setMetaValue("is_chromatogram", "true"); //this is a hack to store that we have chromatogram data
    LayerData::ExperimentType::SpectrumType spectrum;

    // retrieve chromatogram (either from in-memory or on-disc representation)
    MSChromatogram current_chrom;
    current_chrom = exp_sptr->getChromatograms()[index];
    if (current_chrom.empty() )
    {
      current_chrom = ondisc_sptr->getChromatogram(index);
    }

    // fill "dummy" spectrum with chromatogram data
    for (Size i = 0; i != current_chrom.size(); ++i)
    {
      const ChromatogramPeak & cpeak = current_chrom[i];
      Peak1D peak1d;
      peak1d.setMZ(cpeak.getRT());
      peak1d.setIntensity(cpeak.getIntensity());
      spectrum.push_back(peak1d);
    }

    // Add at least one data point to the chromatogram, otherwise
    // "addLayer" will fail and a segfault occurs later
    if (current_chrom.empty()) 
    {
      Peak1D peak1d(-1, 0);
      spectrum.push_back(peak1d);
    }

    // store peptide_sequence if available
    if (current_chrom.getPrecursor().metaValueExists("peptide_sequence"))
    {
      chrom_exp_sptr->setMetaValue("peptide_sequence", current_chrom.getPrecursor().getMetaValue("peptide_sequence"));
    }

    chrom_exp_sptr->addSpectrum(spectrum);
    return chrom_exp_sptr;
  }

  TOPPViewSpectraViewBehavior::TOPPViewSpectraViewBehavior(TOPPViewBase * parent) :
    tv_(parent)
  {
  }

  String caption;

  void TOPPViewSpectraViewBehavior::addPeakMZToImportantPeaks_()
  {
    SpectrumCanvas* current_canvas = tv_->getActive1DWidget()->canvas();
    LayerData& current_layer = current_canvas->getCurrentLayer();
    const SpectrumType& current_spectrum = current_layer.getCurrentSpectrum();
    Size current_spectrum_index = current_layer.getCurrentSpectrumIndex();

    // find interesting peaks
    MSSpectrum spec{current_spectrum};

    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakSpectrum(spec);

    // deisotope so we don't consider higher isotopic peaks
    Deisotoper::deisotopeAndSingleCharge(spec, 
      50,     // tolerance
      true,   // ppm 
      1, 6,   // min / max charge 
      false,  // keep only deisotoped
      3, 10,  // min / max isopeaks 
      false,  // don't convert fragment m/z to mono-charge
      true);  // annotate charge in integer data array

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 50.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 2, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "slide", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);
    window_mower_filter.filterPeakSpectrum(spec);

    NLargest nlargest_filter(10);  // maximum number of annotated m/z's
    nlargest_filter.filterPeakSpectrum(spec);
    spec.sortByPosition(); // nlargest changes order

#ifdef DEBUG_IDENTIFICATION_VIEW
    cout << "Interesting peaks: " << current_peak_index << endl;
#endif    
    for (Size i = 0; i != spec.size(); ++i)
    {
      Size current_peak_index = current_spectrum.findNearest(spec[i].getMZ());

      QString s = QString::number(current_spectrum[current_peak_index].getMZ());
      if (!spec.getIntegerDataArrays().empty() 
        && spec.getIntegerDataArrays()[0].size() == spec.size())
      {
        int charge = spec.getIntegerDataArrays()[0][i]; 
        // TODO: handle negative mode

        // here we explicitly also annotate singly charged ions to distinguish them from unknown charge (0)
        if (charge != 0) 
        {
          s += charge == 1 ? "<sup>+</sup>" : "<sup>" + QString::number(charge) + "+</sup>";
        }
      }
      
      PeakIndex pi(current_spectrum_index, current_peak_index);
      Annotation1DItem* item = tv_->getActive1DWidget()->canvas()->addPeakAnnotation(pi, s, Qt::darkGray);
      temporary_annotations_.push_back(item); // save label for later removal
    }
  }

  void TOPPViewSpectraViewBehavior::removeTemporaryAnnotations_(Size spectrum_index)
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

  void TOPPViewSpectraViewBehavior::showSpectrumAs1D(int index)
  {
    // basic behavior 1
    LayerData & layer = const_cast<LayerData&>(tv_->getActiveCanvas()->getCurrentLayer());
    ExperimentSharedPtrType exp_sptr = layer.getPeakDataMuteable();
    LayerData::ODExperimentSharedPtrType od_exp_sptr = layer.getOnDiscPeakData();
    auto ondisc_sptr = layer.getOnDiscPeakData();

    // open new 1D widget
    Spectrum1DWidget * w = new Spectrum1DWidget(tv_->getSpectrumParameters(1), (QWidget *)tv_->getWorkspace());

    if (layer.type == LayerData::DT_CHROMATOGRAM)
    {
      ExperimentSharedPtrType chrom_exp_sptr = prepareChromatogram(index, exp_sptr, ondisc_sptr);

      // fix legend and set layer name
      caption = layer.name + "[" + index + "]";
      w->xAxis()->setLegend("Time [s]");

      // add chromatogram data as peak spectrum
      if (!w->canvas()->addLayer(chrom_exp_sptr, ondisc_sptr, layer.filename))
      {
        return;
      }
      w->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      w->canvas()->activateSpectrum(0);
    }
    else if (layer.type == LayerData::DT_PEAK)
    {
      caption = layer.name;

      // add data
      if (!w->canvas()->addLayer(exp_sptr, od_exp_sptr, layer.filename) 
        || (Size)index >= w->canvas()->getCurrentLayer().getPeakData()->size())
      {
        return;
      }

      w->canvas()->activateSpectrum(index);
    }
    else
    {
      // Behavior if its neither (user may have clicked on an empty tree or a
      // dummy entry as drawn by SpectraViewWidget::updateEntries)
      QMessageBox::critical(w, "Error", "Cannot open data that is neither chromatogram nor spectrum data. Aborting!");
      return;
    }

    // set relative (%) view of visible area
    w->canvas()->setIntensityMode(SpectrumCanvas::IM_SNAP);

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

    // show multiple spectra together is only used for chromatograms directly
    // where multiple (SRM) traces are shown together
    LayerData & layer = const_cast<LayerData&>(tv_->getActiveCanvas()->getCurrentLayer());
    ExperimentSharedPtrType exp_sptr = layer.getPeakDataMuteable();
    auto ondisc_sptr = layer.getOnDiscPeakData();

    // string for naming the different chromatogram layers with their index
    String chromatogram_caption;
    // string for naming the tab title with the indices of the chromatograms
    caption = layer.name;

    //open new 1D widget
    Spectrum1DWidget * w = new Spectrum1DWidget(tv_->getSpectrumParameters(1), (QWidget *)tv_->getWorkspace());
    // fix legend if its a chromatogram
    w->xAxis()->setLegend("Time [s]");

    for (const auto& index : indices)
    {
      if (layer.type == LayerData::DT_CHROMATOGRAM)
      {
        ExperimentSharedPtrType chrom_exp_sptr = prepareChromatogram(index, exp_sptr, ondisc_sptr);

        // fix legend and set layer name
        caption = caption + " [" + index + "];";
        chromatogram_caption = layer.name + "[" + index + "]";

        // add chromatogram data as peak spectrum
        if (!w->canvas()->addLayer(chrom_exp_sptr, ondisc_sptr, layer.filename))
        {
          return;
        }
        w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), chromatogram_caption);
        w->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);

        w->canvas()->getCurrentLayer().getChromatogramData() = exp_sptr; // save the original chromatogram data so that we can access it later

        //this is a hack to store that we have chromatogram data, that we selected multiple ones and which one we selected
        w->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("is_chromatogram", "true");
        w->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("multiple_select", "true");
        w->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("selected_chromatogram", index);

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

    // return if no active 1D widget is present or no layers are present (e.g. the addLayer call failed)
    if (widget_1d == nullptr 
      || widget_1d->canvas() == nullptr 
      || widget_1d->canvas()->getLayerCount() == 0)
    {
      return;
    }

    // remove old annotations (if present)
    SpectrumCanvas* current_canvas = widget_1d->canvas();
    LayerData& current_layer = current_canvas->getCurrentLayer();
    Size current_spectrum_index = current_layer.getCurrentSpectrumIndex();

    removeTemporaryAnnotations_(current_spectrum_index);

    // activate new spectrum
    widget_1d->canvas()->activateSpectrum(index);
    LayerData & layer = const_cast<LayerData&>(tv_->getActiveCanvas()->getCurrentLayer());

    addPeakMZToImportantPeaks_();

    // If we have a chromatogram, we cannot just simply activate this spectrum.
    // we have to do much more work, e.g. creating a new experiment with the
    // new spectrum.
    if (layer.chromatogram_flag_set())
    {
      // first get raw data (the full experiment with all chromatograms), we
      // only need to grab the ones with the desired indices
      ExperimentSharedPtrType exp_sptr = widget_1d->canvas()->getCurrentLayer().getChromatogramData();
      auto ondisc_sptr = layer.getOnDiscPeakData();

      const LayerData & layer = widget_1d->canvas()->getCurrentLayer();
      String fname = layer.filename;
      String lname = layer.name;

      Size layercount = widget_1d->canvas()->getLayerCount();
      for (Size i = 0; i != layercount; ++i)
      {
        widget_1d->canvas()->removeLayer(0); // remove layer 0 until there are no more layers
      }

      ExperimentSharedPtrType chrom_exp_sptr = prepareChromatogram(index, exp_sptr, ondisc_sptr);

      // fix legend and set layer name
      caption = lname + "[" + index + "]";

      // add chromatogram data as peak spectrum
      if (!widget_1d->canvas()->addLayer(chrom_exp_sptr, ondisc_sptr, fname))
      {
        return;
      }

      widget_1d->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
      widget_1d->canvas()->setIntensityMode(Spectrum1DCanvas::IM_NONE);

      widget_1d->canvas()->getCurrentLayer().name = caption;
      widget_1d->canvas()->getCurrentLayer().filename = fname;
      widget_1d->canvas()->getCurrentLayer().getChromatogramData() = exp_sptr; // save the original chromatogram data so that we can access it later
      //this is a hack to store that we have chromatogram data, that we selected multiple ones and which one we selected
      widget_1d->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("is_chromatogram", "true");
      widget_1d->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("multiple_select", "false");
      widget_1d->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("selected_chromatogram", index);

      tv_->updateLayerBar();
      tv_->updateViewBar();
      tv_->updateFilterBar();
      tv_->updateMenu();
    }
  }

  void TOPPViewSpectraViewBehavior::activate1DSpectrum(std::vector<int, std::allocator<int> > indices)
  {
    Spectrum1DWidget * widget_1d = tv_->getActive1DWidget();

    // return if no active 1D widget is present or no layers are present (e.g. the addLayer call failed)
    if (widget_1d == nullptr) return;
    if (widget_1d->canvas()->getLayerCount() == 0) return;

    // If we have a chromatogram, we cannot just simply activate this spectrum.
    // we have to do much more work, e.g. creating a new experiment with the
    // new spectrum.
    const LayerData & layer = widget_1d->canvas()->getCurrentLayer();
    String fname = layer.filename;
    if (layer.chromatogram_flag_set())
    {
      // first get raw data (the full experiment with all chromatograms), we
      // only need to grab the ones with the desired indices
      ExperimentSharedPtrType exp_sptr = widget_1d->canvas()->getCurrentLayer().getChromatogramData();
      auto ondisc_sptr = layer.getOnDiscPeakData();

      Size layercount = widget_1d->canvas()->getLayerCount();
      for (Size i = 0; i != layercount; ++i)
      {
        widget_1d->canvas()->removeLayer(0); // remove layer 0 until there are no more layers
      }

      for (const auto& index : indices)
      {
        ExperimentSharedPtrType chrom_exp_sptr = prepareChromatogram(index, exp_sptr, ondisc_sptr);

        // get caption (either chromatogram idx or peptide sequence, if available)
        caption = fname + "[" + index + "]";
        if (chrom_exp_sptr->metaValueExists("peptide_sequence"))
        {
          caption = String(chrom_exp_sptr->getMetaValue("peptide_sequence")) + "[" + index + "]";
        }

        // add chromatogram data as peak spectrum
        if (!widget_1d->canvas()->addLayer(chrom_exp_sptr, ondisc_sptr, fname))
        {
          return;
        }

        widget_1d->canvas()->setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
        widget_1d->canvas()->setIntensityMode(Spectrum1DCanvas::IM_NONE);

        widget_1d->canvas()->getCurrentLayer().name = caption;
        widget_1d->canvas()->getCurrentLayer().filename = fname;
        widget_1d->canvas()->getCurrentLayer().getChromatogramData() = exp_sptr; // save the original chromatogram data so that we can access it later
        //this is a hack to store that we have chromatogram data, that we selected multiple ones and which one we selected
        widget_1d->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("is_chromatogram", "true");
        widget_1d->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("multiple_select", "true");
        widget_1d->canvas()->getCurrentLayer().getPeakDataMuteable()->setMetaValue("selected_chromatogram", index);
      }

      tv_->updateLayerBar();
      tv_->updateViewBar();
      tv_->updateFilterBar();
      tv_->updateMenu();
    }
  }

  void TOPPViewSpectraViewBehavior::deactivate1DSpectrum(int spectrum_index)
  {
    // no special handling of activation
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

