// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/TOPPViewIdentificationViewBehavior.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DCaret.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <QtGui/QMessageBox>
#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  TOPPViewIdentificationViewBehavior::TOPPViewIdentificationViewBehavior(TOPPViewBase * parent) :
    tv_(parent)
  {
  }

  void TOPPViewIdentificationViewBehavior::showSpectrumAs1D(int index)
  {
    // basic behavior 1
    const LayerData & layer = tv_->getActiveCanvas()->getCurrentLayer();
    ExperimentSharedPtrType exp_sptr = layer.getPeakData();


    if (layer.type == LayerData::DT_PEAK)
    {
      // open new 1D widget with the current default parameters
      Spectrum1DWidget* w = new Spectrum1DWidget(tv_->getSpectrumParameters(1), (QWidget *)tv_->getWorkspace());
      // add data
      if (!w->canvas()->addLayer(exp_sptr, layer.filename) || (Size)index >= w->canvas()->getCurrentLayer().getPeakData()->size())
      {
        return;
      }

      w->canvas()->activateSpectrum(index);

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

      String caption = layer.name;
      w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);

      tv_->showSpectrumWidgetInWindow(w, caption);

      // special behavior
      const vector<PeptideIdentification>& pi = w->canvas()->getCurrentLayer().getCurrentSpectrum().getPeptideIdentifications();
      if (!pi.empty())
      {
        // mass fingerprint annotation of name etc 
        if (ms_level == 1) addPeakAnnotations_(pi);
        PeptideHit hit;
        if (IDFilter().getBestHit(pi, false, hit))
        {
          addTheoreticalSpectrumLayer_(hit);
        }
        else
        {
          LOG_ERROR << "Spectrum has no hits" << std::endl;
        }
      }

      tv_->updateLayerBar();
      tv_->updateViewBar();
      tv_->updateFilterBar();
      tv_->updateMenu();
    }
    else if (layer.type == LayerData::DT_CHROMATOGRAM)
    {

    }

  }

  void TOPPViewIdentificationViewBehavior::addPeakAnnotations_(const std::vector<PeptideIdentification>& ph)
  {
    // called anew for every click on a spectrum
    LayerData& current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();

    if (current_layer.getCurrentSpectrum().empty())
    {
      LOG_WARN << "Spectrum is empty! Nothing to annotate!" << std::endl;
    }

    // mass precision to match a peak's m/z to a feature m/z
    // m/z values of features are usually an average over multiple scans...
    double ppm = 0.5;

    std::vector<QColor> cols;
    cols.push_back(Qt::blue);
    cols.push_back(Qt::green);
    cols.push_back(Qt::red);
    cols.push_back(Qt::gray);
    cols.push_back(Qt::darkYellow);


    if (!current_layer.getCurrentSpectrum().isSorted())
    {
      QMessageBox::warning(tv_, "Error", "The spectrum is not sorted! Aborting!");
      return;
    }
    for (std::vector<PeptideIdentification>::const_iterator it = ph.begin();
                                                            it!= ph.end();
                                                            ++it)
    {
      if (!it->hasMZ()) continue;
      double mz = it->getMZ();
      Size peak_idx = current_layer.getCurrentSpectrum().findNearest(mz);
        
      // m/z fits ?
      if ( abs(mz - current_layer.getCurrentSpectrum()[peak_idx].getMZ()) / mz * 1e6 > ppm) continue;

      double peak_int = current_layer.getCurrentSpectrum()[peak_idx].getIntensity();

      Annotation1DCaret* first_dit(0);
      // we could have many many hits for different compounds which have the exact same sum formula... so first group by sum formula
      std::map<String, StringList> formula_to_names;
      for (std::vector< PeptideHit >::const_iterator ith = it->getHits().begin();
                                                      ith!= it->getHits().end();
                                                      ++ith)
      {
        if (ith->metaValueExists("identifier") && ith->metaValueExists("chemical_formula"))
        {
          String name = ith->getMetaValue("identifier");
          if (name.length() > 20) 
          {
            name = name.substr(0, 17) + "...";
          }
          formula_to_names[ith->getMetaValue("chemical_formula")].push_back(name);
        }
        else
        {
          StringList msg;
          if (!ith->metaValueExists("identifier")) msg.push_back("identifier");
          if (!ith->metaValueExists("chemical_formula")) msg.push_back("chemical_formula");
          LOG_WARN << "Missing meta-value(s): " << ListUtils::concatenate(msg, ", ") << ". Cannot annotate!\n";
        }
      }

      // assemble annotation (each formula gets a paragraph)
      String text = "<html><body>";
      Size i(0);
      for (std::map<String, StringList>::iterator ith = formula_to_names.begin();
                                                  ith!= formula_to_names.end();
                                                  ++ith)
      {
        if (++i >= 4)
        { // at this point, this is the 4th entry.. which we don't show any more...
          text += String("<b><span style=\"color:") + cols[i].name() + "\">..." + Size(std::distance(formula_to_names.begin(), formula_to_names.end()) - 4 + 1) + " more</span></b><br>";
          break;
        }
        text += String("<b><span style=\"color:") + cols[i].name() + "\">" + ith->first + "</span></b><br>\n";
        // carets for isotope profile
        EmpiricalFormula ef(ith->first);
        IsotopeDistribution id = ef.getIsotopeDistribution(3); // three isotopes at most
        double int_factor = peak_int / id.begin()->second;
        Annotation1DCaret::PositionsType points;
        Size itic(0);
        for (IsotopeDistribution::ConstIterator iti = id.begin(); iti != id.end(); ++iti)
        {
          points.push_back(Annotation1DCaret::PointType(mz + itic*Constants::C13C12_MASSDIFF_U, iti->second * int_factor));
          ++itic;
        }
        Annotation1DCaret* ditem = new Annotation1DCaret(points,
                                                          QString(),
                                                          cols[i]);
        ditem->setSelected(false);
        temporary_annotations_.push_back(ditem); // for removal (no ownership)
        current_layer.getCurrentAnnotations().push_front(ditem); // for visualization (ownership)
        if (first_dit==0) first_dit = ditem; // remember first item (we append the text, when ready)

        // list of compound names  (shorten if required)
        if (ith->second.size() > 3)
        {
          Size s = ith->second.size();
          ith->second[3] = String("...") + (s-3) + " more";
          ith->second.resize(4);
        }
        text += " - " + ListUtils::concatenate(ith->second, "<br> - ") + "<br>\n";
      }
      text += "</body></html>";
      if (first_dit!=0)
      {
        first_dit->setRichText(text.toQString());
      }
    }
  }

  void TOPPViewIdentificationViewBehavior::activate1DSpectrum(int index)
  {
    Spectrum1DWidget * widget_1D = tv_->getActive1DWidget();
    widget_1D->canvas()->activateSpectrum(index);
    const LayerData & current_layer = widget_1D->canvas()->getCurrentLayer();

    if (current_layer.type == LayerData::DT_PEAK)
    {
      UInt ms_level = current_layer.getCurrentSpectrum().getMSLevel();

      if (ms_level == 2) // show theoretical spectrum with automatic alignment
      {
        vector<PeptideIdentification> pi = current_layer.getCurrentSpectrum().getPeptideIdentifications();
        if (!pi.empty())
        {
          PeptideHit hit;
          if (IDFilter().getBestHit(pi, false, hit)) addTheoreticalSpectrumLayer_(hit);
          else LOG_ERROR << "Spectrum has no hits\n";
        }
      }
      else if (ms_level == 1)   // show precursor locations
      {
        const vector<PeptideIdentification>& pi = current_layer.getCurrentSpectrum().getPeptideIdentifications();
        addPeakAnnotations_(pi);

        vector<Precursor> precursors;
        // collect all MS2 spectra precursor till next MS1 spectrum is encountered
        for (Size i = index + 1; i < current_layer.getPeakData()->size(); ++i)
        {
          if ((*current_layer.getPeakData())[i].getMSLevel() == 1)
          {
            break;
          }
          // skip MS2 without precursor
          if ((*current_layer.getPeakData())[i].getPrecursors().empty())
          {
            continue;
          }
          // there should be only one precursor per MS2 spectrum.
          vector<Precursor> pcs = (*current_layer.getPeakData())[i].getPrecursors();
          copy(pcs.begin(), pcs.end(), back_inserter(precursors));
        }
        addPrecursorLabels1D_(precursors);
      }
    } // end DT_PEAK
    else if (current_layer.type == LayerData::DT_CHROMATOGRAM)
    {

    }
  }

  void TOPPViewIdentificationViewBehavior::addPrecursorLabels1D_(const vector<Precursor> & pcs)
  {
    LayerData & current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();

    if (current_layer.type == LayerData::DT_PEAK)
    {
      const SpectrumType& spectrum = current_layer.getCurrentSpectrum();

      for (vector<Precursor>::const_iterator it = pcs.begin(); it != pcs.end(); ++it)
      {
        // determine start and stop of isolation window
        double isolation_window_lower_mz = it->getMZ() - it->getIsolationWindowLowerOffset();
        double isolation_window_upper_mz = it->getMZ() + it->getIsolationWindowUpperOffset();

        // determine maximum peak intensity in isolation window
        SpectrumType::const_iterator vbegin = spectrum.MZBegin(isolation_window_lower_mz);
        SpectrumType::const_iterator vend = spectrum.MZEnd(isolation_window_upper_mz);

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

        Annotation1DDistanceItem * item = new Annotation1DDistanceItem(QString::number(it->getCharge()), lower_position, upper_position);
        // add additional tick at precursor target position (e.g. to show if isolation window is assymetric)
        vector<double> ticks;
        ticks.push_back(it->getMZ());
        item->setTicks(ticks);
        item->setSelected(false);

        temporary_annotations_.push_back(item); // for removal (no ownership)
        current_layer.getCurrentAnnotations().push_front(item); // for visualisation (ownership)
      }
    }
    else if (current_layer.type == LayerData::DT_CHROMATOGRAM)
    {

    }
  }

  /// Behavior for activate1DSpectrum
  void TOPPViewIdentificationViewBehavior::activate1DSpectrum(std::vector<int, std::allocator<int> >)
  {
  }

  void TOPPViewIdentificationViewBehavior::removeTemporaryAnnotations_(Size spectrum_index)
  {
#ifdef DEBUG_IDENTIFICATION_VIEW
    cout << "removePrecursorLabels1D_ " << spectrum_index << endl;
#endif
    // Delete annotations added by IdentificationView (but not user added annotations)
    LayerData & current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();
    const vector<Annotation1DItem *> & cas = temporary_annotations_;
    Annotations1DContainer & las = current_layer.getAnnotations(spectrum_index);
    for (vector<Annotation1DItem *>::const_iterator it = cas.begin(); it != cas.end(); ++it)
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

  void TOPPViewIdentificationViewBehavior::addTheoreticalSpectrumLayer_(const PeptideHit & ph)
  {
    SpectrumCanvas * current_canvas = tv_->getActive1DWidget()->canvas();
    LayerData & current_layer = current_canvas->getCurrentLayer();
    SpectrumType & current_spectrum = current_layer.getCurrentSpectrum();

    AASequence aa_sequence = ph.getSequence();

    // get measured spectrum indices and spectrum
    Size current_spectrum_layer_index = current_canvas->activeLayerIndex();
    Size current_spectrum_index = current_layer.getCurrentSpectrumIndex();

    const Param & tv_params = tv_->getParameters();

    RichPeakSpectrum rich_spec;
    TheoreticalSpectrumGenerator generator;
    Param p;
    p.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");

    p.setValue("max_isotope", tv_params.getValue("preferences:idview:max_isotope"), "Number of isotopic peaks");
    p.setValue("add_losses", tv_params.getValue("preferences:idview:add_losses"), "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    p.setValue("add_isotopes", tv_params.getValue("preferences:idview:add_isotopes"), "If set to 1 isotope peaks of the product ion peaks are added");
    p.setValue("add_abundant_immonium_ions", tv_params.getValue("preferences:idview:add_abundant_immonium_ions"), "Add most abundant immonium ions");

    p.setValue("a_intensity", current_spectrum.getMaxInt() * (double)tv_params.getValue("preferences:idview:a_intensity"), "Intensity of the a-ions");
    p.setValue("b_intensity", current_spectrum.getMaxInt() * (double)tv_params.getValue("preferences:idview:b_intensity"), "Intensity of the b-ions");
    p.setValue("c_intensity", current_spectrum.getMaxInt() * (double)tv_params.getValue("preferences:idview:c_intensity"), "Intensity of the c-ions");
    p.setValue("x_intensity", current_spectrum.getMaxInt() * (double)tv_params.getValue("preferences:idview:x_intensity"), "Intensity of the x-ions");
    p.setValue("y_intensity", current_spectrum.getMaxInt() * (double)tv_params.getValue("preferences:idview:y_intensity"), "Intensity of the y-ions");
    p.setValue("z_intensity", current_spectrum.getMaxInt() * (double)tv_params.getValue("preferences:idview:z_intensity"), "Intensity of the z-ions");
    p.setValue("relative_loss_intensity", tv_params.getValue("preferences:idview:relative_loss_intensity"), "Intensity of loss ions, in relation to the intact ion intensity");
    generator.setParameters(p);

    try
    {
      Int max_charge = max(1, ph.getCharge()); // at least generate charge 1 if no charge (0) is annotated

      // generate mass ladder for each charge state
      for (Int charge = 1; charge <= max_charge; ++charge)
      {
        if (tv_params.getValue("preferences:idview:show_a_ions").toBool()) // "A-ions"
        {
          generator.addPeaks(rich_spec, aa_sequence, Residue::AIon, charge);
        }
        if (tv_params.getValue("preferences:idview:show_b_ions").toBool()) // "B-ions"
        {
          generator.addPeaks(rich_spec, aa_sequence, Residue::BIon, charge);
        }
        if (tv_params.getValue("preferences:idview:show_c_ions").toBool()) // "C-ions"
        {
          generator.addPeaks(rich_spec, aa_sequence, Residue::CIon, charge);
        }
        if (tv_params.getValue("preferences:idview:show_x_ions").toBool()) // "X-ions"
        {
          generator.addPeaks(rich_spec, aa_sequence, Residue::XIon, charge);
        }
        if (tv_params.getValue("preferences:idview:show_y_ions").toBool()) // "Y-ions"
        {
          generator.addPeaks(rich_spec, aa_sequence, Residue::YIon, charge);
        }
        if (tv_params.getValue("preferences:idview:show_z_ions").toBool()) // "Z-ions"
        {
          generator.addPeaks(rich_spec, aa_sequence, Residue::ZIon, charge);
        }
        if (tv_params.getValue("preferences:idview:show_precursor").toBool()) // "Precursor"
        {
          generator.addPrecursorPeaks(rich_spec, aa_sequence, charge);
        }
      }
      if (tv_params.getValue("preferences:idview:add_abundant_immonium_ions").toBool()) // "abundant Immonium-ions"
      {
        generator.addAbundantImmoniumIons(rich_spec);
      }
    }
    catch (Exception::BaseException & e)
    {
      QMessageBox::warning(tv_, "Error", QString("Spectrum generation failed! (") + e.what() + "). Please report this to the developers (specify what input you used)!");
      return;
    }

    // convert rich spectrum to simple spectrum
    PeakSpectrum new_spec;
    for (RichPeakSpectrum::Iterator it = rich_spec.begin(); it != rich_spec.end(); ++it)
    {
      new_spec.push_back(static_cast<Peak1D>(*it));
    }

    PeakMap new_exp;
    new_exp.addSpectrum(new_spec);
    ExperimentSharedPtrType new_exp_sptr(new PeakMap(new_exp));
    FeatureMapSharedPtrType f_dummy(new FeatureMapType());
    ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
    vector<PeptideIdentification> p_dummy;

    // Block update events for identification widget
    tv_->getSpectraIdentificationViewWidget()->ignore_update = true;

    String layer_caption = aa_sequence.toString().toQString() + QString(" (identification view)");
    tv_->addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, LayerData::DT_CHROMATOGRAM, false, false, false, "", layer_caption.toQString());

    // get layer index of new layer
    Size theoretical_spectrum_layer_index = tv_->getActive1DWidget()->canvas()->activeLayerIndex();

    // kind of a hack to check whether adding the layer was successful
    if (current_spectrum_layer_index != theoretical_spectrum_layer_index)
    {
      // Ensure theoretical spectrum is drawn as dashed sticks
      tv_->setDrawMode1D(Spectrum1DCanvas::DM_PEAKS);
      tv_->getActive1DWidget()->canvas()->setCurrentLayerPeakPenStyle(Qt::DashLine);

      // Add ion names as annotations to the theoretical spectrum
      for (RichPeakSpectrum::Iterator it = rich_spec.begin(); it != rich_spec.end(); ++it)
      {
        if (it->getMetaValue("IonName") != DataValue::EMPTY)
        {
          DPosition<2> position = DPosition<2>(it->getMZ(), it->getIntensity());
          QString s(((string)it->getMetaValue("IonName")).c_str());

          if (s.at(0) == 'y')
          {
            Annotation1DItem * item = new Annotation1DPeakItem(position, s, Qt::darkRed);
            item->setSelected(false);
            tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
          }
          else if (s.at(0) == 'b')
          {
            Annotation1DItem * item = new Annotation1DPeakItem(position, s, Qt::darkGreen);
            item->setSelected(false);
            tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
          }
        }
      }

      // remove theoretical and activate real data layer and spectrum
      tv_->getActive1DWidget()->canvas()->changeVisibility(theoretical_spectrum_layer_index, false);
      tv_->getActive1DWidget()->canvas()->activateLayer(current_spectrum_layer_index);
      tv_->getActive1DWidget()->canvas()->getCurrentLayer().setCurrentSpectrumIndex(current_spectrum_index);

      // zoom to maximum visible area in real data (as theoretical might be much larger and therefor squeezes the interesting part)
      DRange<2> visible_area = tv_->getActive1DWidget()->canvas()->getVisibleArea();
      double min_mz = tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentSpectrum().getMin()[0];
      double max_mz = tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentSpectrum().getMax()[0];
      double delta_mz = max_mz - min_mz;
      visible_area.setMin(min_mz - 0.1 * delta_mz);
      visible_area.setMax(max_mz + 0.1 * delta_mz);

      tv_->getActive1DWidget()->canvas()->setVisibleArea(visible_area);

      // spectra alignment
      Param param;

      double tolerance = tv_params.getValue("preferences:idview:tolerance");

      param.setValue("tolerance", tolerance, "Defines the absolute (in Da) or relative (in ppm) tolerance in the alignment");
      tv_->getActive1DWidget()->performAlignment(current_spectrum_layer_index, theoretical_spectrum_layer_index, param);

      std::vector<std::pair<Size, Size> > aligned_peak_indices = tv_->getActive1DWidget()->canvas()->getAlignedPeaksIndices();

      // annotate original spectrum with ions and sequence
      for (Size i = 0; i != aligned_peak_indices.size(); ++i)
      {
        PeakIndex pi(current_spectrum_index, aligned_peak_indices[i].first);
        QString s(((string)rich_spec[aligned_peak_indices[i].second].getMetaValue("IonName")).c_str());
        QString ion_nr_string = s;

        if (s.at(0) == 'y')
        {
          ion_nr_string.replace("y", "");
          ion_nr_string.replace("+", "");
          Size ion_number = ion_nr_string.toUInt();
          s.append("\n");
          // extract peptide ion sequence
          QString aa_ss;
          for (Size j = aa_sequence.size() - 1; j >= aa_sequence.size() - ion_number; --j)
          {
            const Residue & r = aa_sequence.getResidue(j);
            aa_ss.append(r.getOneLetterCode().toQString());
            if (r.getModification() != "")
            {
              aa_ss.append("*");
            }
          }
          s.append(aa_ss);
          Annotation1DItem * item = tv_->getActive1DWidget()->canvas()->addPeakAnnotation(pi, s, Qt::darkRed);
          temporary_annotations_.push_back(item);
        }
        else if (s.at(0) == 'b')
        {
          ion_nr_string.replace("b", "");
          ion_nr_string.replace("+", "");
          UInt ion_number = ion_nr_string.toUInt();
          s.append("\n");
          // extract peptide ion sequence
          AASequence aa_subsequence = aa_sequence.getSubsequence(0, ion_number);
          QString aa_ss = aa_subsequence.toString().toQString();
          // shorten modifications "(MODNAME)" to "*"
          aa_ss.replace(QRegExp("[(].*[)]"), "*");
          // append to label
          s.append(aa_ss);
          Annotation1DItem * item = tv_->getActive1DWidget()->canvas()->addPeakAnnotation(pi, s, Qt::darkGreen);
          // save label for later removal
          temporary_annotations_.push_back(item);
        } else
        {
          s.append("\n");
          Annotation1DItem * item = tv_->getActive1DWidget()->canvas()->addPeakAnnotation(pi, s, Qt::black);
          // save label for later removal
          temporary_annotations_.push_back(item);
        }
      }

      tv_->updateLayerBar();
      tv_->getSpectraIdentificationViewWidget()->ignore_update = false;
    }
  }

  void TOPPViewIdentificationViewBehavior::deactivate1DSpectrum(int spectrum_index)
  {
    LayerData & current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();
    int ms_level = (*current_layer.getPeakData())[spectrum_index].getMSLevel();

    removeTemporaryAnnotations_(spectrum_index);

    if (ms_level == 2)
    {
      removeTheoreticalSpectrumLayer_();
    }

    // the next line is meant to be disabled to allow switching between spectra without loosing the current view range (to compare across spectra)
    // tv_->getActive1DWidget()->canvas()->resetZoom();
  }

  void TOPPViewIdentificationViewBehavior::removeTheoreticalSpectrumLayer_()
  {
    Spectrum1DWidget * spectrum_widget_1D = tv_->getActive1DWidget();
    if (spectrum_widget_1D)
    {
      Spectrum1DCanvas * canvas_1D = spectrum_widget_1D->canvas();

      // Find the automatically generated layer with theoretical spectrum and remove it and the associated alignment.
      // before activating the next normal spectrum
      Size lc = canvas_1D->getLayerCount();
      for (Size i = 0; i != lc; ++i)
      {
        String ln = canvas_1D->getLayerName(i);
        if (ln.hasSubstring("(identification view)"))
        {
          canvas_1D->removeLayer(i);
          canvas_1D->resetAlignment();
          tv_->updateLayerBar();
          break;
        }
      }
    }
  }

  void TOPPViewIdentificationViewBehavior::activateBehavior()
  {
    Spectrum1DWidget* w = tv_->getActive1DWidget();
    if ( w == 0)
    {
      return;
    }
    SpectrumCanvas * current_canvas = w->canvas();
    LayerData & current_layer = current_canvas->getCurrentLayer();
    SpectrumType & current_spectrum = current_layer.getCurrentSpectrum();

    // find first MS2 spectrum with peptide identification and set current spectrum to it
    if (current_spectrum.getMSLevel() == 1)  // no fragment spectrum
    {
      for (Size i = 0; i < current_layer.getPeakData()->size(); ++i)
      {
        UInt ms_level = (*current_layer.getPeakData())[i].getMSLevel();
        const vector<PeptideIdentification> peptide_ids = (*current_layer.getPeakData())[i].getPeptideIdentifications();
        Size peptide_ids_count = peptide_ids.size();

        if (ms_level != 2 || peptide_ids_count == 0)  // skip non ms2 spectra and spectra with no identification
        {
          continue;
        }
        current_layer.setCurrentSpectrumIndex(i);
        break;
      }
    }
  }

  void TOPPViewIdentificationViewBehavior::deactivateBehavior()
  {
    // remove precusor labels, theoretical spectra and trigger repaint
    if (tv_->getActive1DWidget() != 0)
    {
      removeTemporaryAnnotations_(tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentSpectrumIndex());
      removeTheoreticalSpectrumLayer_();
      tv_->getActive1DWidget()->canvas()->repaint();
    }
  }

  void TOPPViewIdentificationViewBehavior::setVisibleArea1D(double l, double h)
  {
    if (tv_->getActive1DWidget() != 0)
    {
      DRange<2> range = tv_->getActive1DWidget()->canvas()->getVisibleArea();
      range.setMinX(l);
      range.setMaxX(h);
      tv_->getActive1DWidget()->canvas()->setVisibleArea(range);
      tv_->getActive1DWidget()->canvas()->repaint();
    }
  }

}
