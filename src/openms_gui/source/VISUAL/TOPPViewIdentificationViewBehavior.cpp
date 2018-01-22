// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Container.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopeDistribution.h>

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
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <QtCore/QString>
#include <QtGui/QMessageBox>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  TOPPViewIdentificationViewBehavior::TOPPViewIdentificationViewBehavior(TOPPViewBase* parent) :
    tv_(parent)
  {
  }

  void TOPPViewIdentificationViewBehavior::showSpectrumAs1D(int index)
  {
    // call without selecting an identification
    showSpectrumAs1D(index, -1, -1);
  }

  void TOPPViewIdentificationViewBehavior::showSpectrumAs1D(int spectrum_index, int peptide_id_index, int peptide_hit_index)
  {
    // basic behavior 1
    const LayerData& layer = tv_->getActiveCanvas()->getCurrentLayer();
    ExperimentSharedPtrType exp_sptr = layer.getPeakData();

    if (layer.type == LayerData::DT_PEAK)
    {
      // open new 1D widget with the current default parameters
      Spectrum1DWidget* w = new Spectrum1DWidget(tv_->getSpectrumParameters(1), (QWidget*)tv_->getWorkspace());
      // add data
      if (!w->canvas()->addLayer(exp_sptr, layer.filename) || (Size)spectrum_index >= w->canvas()->getCurrentLayer().getPeakData()->size())
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

      String caption = layer.name;
      w->canvas()->setLayerName(w->canvas()->activeLayerIndex(), caption);

      tv_->showSpectrumWidgetInWindow(w, caption);

      // special behavior
      if (peptide_id_index == -1 || peptide_hit_index == -1) { return; }

      const vector<PeptideIdentification>& pis = w->canvas()->getCurrentLayer().getCurrentSpectrum().getPeptideIdentifications();
      if (!pis.empty())
      {
        switch (ms_level)
        {
          // mass fingerprint annotation of name etc
          case 1: { addPeakAnnotations_(pis); break; }
          // annotation with stored fragments or synthesized theoretical spectrum
          case 2:
          {
            // check if index in bounds and hits are present
            if (peptide_id_index < static_cast<int>(pis.size()) && peptide_hit_index < static_cast<int>(pis[peptide_id_index].getHits().size()))
            {
              // get hit
              PeptideHit ph = pis[peptide_id_index].getHits()[peptide_hit_index];
              if (ph.getPeakAnnotations().empty())
              {
                // if no fragment annotations are stored, create a theoretical spectrum
                addTheoreticalSpectrumLayer_(ph);
              }
              else
              {
                // otherwise, use stored fragment annotations
                addAnnotationsSpectrumLayer_(ph);
              }
            }
            break;
          }
          default:
            LOG_WARN << "Annotation of MS level > 2 not supported.!" << std::endl;
        }
      }

      tv_->updateLayerBar();
      tv_->updateViewBar();
      tv_->updateFilterBar();
      tv_->updateMenu();
    }
    // else if (layer.type == LayerData::DT_CHROMATOGRAM)
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
      if (Math::getPPMAbs(mz, current_layer.getCurrentSpectrum()[peak_idx].getMZ()) > ppm) continue;

      double peak_int = current_layer.getCurrentSpectrum()[peak_idx].getIntensity();

      Annotation1DCaret* first_dit(nullptr);
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
          String cf = ith->getMetaValue("chemical_formula");
          if (cf.empty()) continue; // skip unannotated "null" peaks
          formula_to_names[cf].push_back(name);
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
        IsotopeDistribution id = ef.getIsotopeDistribution(new CoarseIsotopeDistribution(3)); // three isotopes at most
        double int_factor = peak_int / id.begin()->getIntensity();
        Annotation1DCaret::PositionsType points;
        Size itic(0);
        for (IsotopeDistribution::ConstIterator iti = id.begin(); iti != id.end(); ++iti)
        {
          points.push_back(Annotation1DCaret::PointType(mz + itic*Constants::C13C12_MASSDIFF_U, iti->getIntensity() * int_factor));
          ++itic;
        }
        Annotation1DCaret* ditem = new Annotation1DCaret(points,
                                                         QString(),
                                                         cols[i],
                                                         current_layer.param.getValue("peak_color").toQString());
        ditem->setSelected(false);
        temporary_annotations_.push_back(ditem); // for removal (no ownership)
        current_layer.getCurrentAnnotations().push_front(ditem); // for visualization (ownership)
        if (first_dit==nullptr) first_dit = ditem; // remember first item (we append the text, when ready)

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
      if (first_dit!=nullptr)
      {
        first_dit->setRichText(text.toQString());
      }
    }
  }

  void TOPPViewIdentificationViewBehavior::activate1DSpectrum(int index)
  {
    activate1DSpectrum(index, -1, -1);
  }

  void TOPPViewIdentificationViewBehavior::activate1DSpectrum(int spectrum_index, int peptide_id_index, int peptide_hit_index)
  {
    Spectrum1DWidget* widget_1D = tv_->getActive1DWidget();

    // return if no active 1D widget is present
    if (widget_1D == nullptr) return;

    widget_1D->canvas()->activateSpectrum(spectrum_index);
    LayerData& current_layer = widget_1D->canvas()->getCurrentLayer();
    current_layer.peptide_id_index = peptide_id_index;
    current_layer.peptide_hit_index = peptide_hit_index;

    if (current_layer.type == LayerData::DT_PEAK)
    {
      UInt ms_level = current_layer.getCurrentSpectrum().getMSLevel();

      const vector<PeptideIdentification>& pis = current_layer.getCurrentSpectrum().getPeptideIdentifications();
      switch (ms_level)
      {
        case 1: // mass fingerprint annotation of name etc and precursor labels
        {
          addPeakAnnotations_(pis);
          vector<Precursor> precursors;

          // collect all MS2 spectra precursor till next MS1 spectrum is encountered
          for (Size i = spectrum_index + 1; i < current_layer.getPeakData()->size(); ++i)
          {
            if ((*current_layer.getPeakData())[i].getMSLevel() == 1) break;

            // skip MS2 without precursor
            if ((*current_layer.getPeakData())[i].getPrecursors().empty()) continue;

            // there should be only one precursor per MS2 spectrum.
            vector<Precursor> pcs = (*current_layer.getPeakData())[i].getPrecursors();
            copy(pcs.begin(), pcs.end(), back_inserter(precursors));
          }
          addPrecursorLabels1D_(precursors);
          break;
        }
        case 2: // annotation with stored fragments or synthesized theoretical spectrum
        {
          // check if index in bounds and hits are present
          if (peptide_id_index < static_cast<int>(pis.size()) && peptide_hit_index < static_cast<int>(pis[peptide_id_index].getHits().size()))
          {
            // get selected hit
            PeptideHit ph = pis[peptide_id_index].getHits()[peptide_hit_index];

            if (ph.getPeakAnnotations().empty())
            {
              // if no fragment annotations are stored, create a theoretical spectrum
              addTheoreticalSpectrumLayer_(ph);
            }
            else
            {
              // otherwise, use stored fragment annotations
              addAnnotationsSpectrumLayer_(ph);

              if (ph.metaValueExists("xl_chain")) // if this meta value exists, this should be an XLMS annotation
              {
                String box_text;
                String vert_bar = "&#124;";

                if (ph.metaValueExists("xl_pos2")) // if this meta value exists, this should be the special case of a loop-link
                {
                  String hor_bar = "_";
                  PeptideHit ph_alpha = pis[peptide_id_index].getHits()[0];
                  String seq_alpha = ph.getSequence().toUnmodifiedString();
                  int xl_pos_alpha = String(ph.getMetaValue("xl_pos")).toInt();
                  int xl_pos_beta = String(ph.getMetaValue("xl_pos2")).toInt() - xl_pos_alpha - 1;

                  String alpha_cov;
                  String beta_cov;
                  extractCoverageStrings(ph.getPeakAnnotations(), alpha_cov, beta_cov, seq_alpha.size(), 0);

                  // String formatting
                  box_text += alpha_cov + "<br>" +  seq_alpha +  "<br>" + String(xl_pos_alpha, ' ') +  vert_bar + n_times(xl_pos_beta, hor_bar) + vert_bar;
                  // cut out line: "<br>" + String(xl_pos_alpha, ' ') + vert_bar + String(xl_pos_beta, ' ') + vert_bar +
                }
                else if (pis[peptide_id_index].getHits().size() == 2) // xl_chain exists and 2 PeptideHits: should be a cross-link
                {
                  PeptideHit ph_alpha = pis[peptide_id_index].getHits()[0];
                  PeptideHit ph_beta = pis[peptide_id_index].getHits()[1];
                  String seq_alpha = ph_alpha.getSequence().toUnmodifiedString();
                  String seq_beta = ph_beta.getSequence().toUnmodifiedString();
                  int xl_pos_alpha = String(ph_alpha.getMetaValue("xl_pos")).toInt();
                  int xl_pos_beta = String(ph_beta.getMetaValue("xl_pos")).toInt();


                  // String formatting
                  Size prefix_length = std::max(xl_pos_alpha, xl_pos_beta);
                  //Size suffix_length = std::max(seq_alpha.size() - xl_pos_alpha, seq_beta.size() - xl_pos_beta);
                  Size alpha_space = prefix_length - xl_pos_alpha;
                  Size beta_space = prefix_length - xl_pos_beta;

                  String alpha_cov;
                  String beta_cov;
                  extractCoverageStrings(ph_alpha.getPeakAnnotations(), alpha_cov, beta_cov, seq_alpha.size(), seq_beta.size());

                  box_text += String(alpha_space, ' ') + alpha_cov + "<br>" + String(alpha_space, ' ') + seq_alpha + "<br>" + String(prefix_length, ' ') + vert_bar + "<br>" + String(beta_space, ' ') + seq_beta + "<br>" + String(beta_space, ' ') + beta_cov;
                  // color: <font color=\"green\">&boxur;</font>
                }
                else // no xl_pos2 and no second PeptideHit, should be a mono-link
                {
                  String seq_alpha = ph.getSequence().toUnmodifiedString();
                  int xl_pos_alpha = String(ph.getMetaValue("xl_pos")).toInt();
                  Size prefix_length = xl_pos_alpha;

                  String alpha_cov;
                  String beta_cov;
                  extractCoverageStrings(ph.getPeakAnnotations(), alpha_cov, beta_cov, seq_alpha.size(), 0);

                  box_text += alpha_cov + "<br>" + seq_alpha + "<br>" + String(prefix_length, ' ') + vert_bar;

                }
                box_text =  "<font size=\"5\" style=\"background-color:white;\"><pre>" + box_text + "</pre></font> ";
                widget_1D->canvas()->setTextBox(box_text.toQString());
              }
              else
              {
                String seq = ph.getSequence().toString();
                if (seq.empty()) seq = ph.getMetaValue("label");
                widget_1D->canvas()->setTextBox(seq.toQString());
              }
            }
          }
          break;
        }
        default:
          LOG_WARN << "Annotation of MS level > 2 not supported.!" << std::endl;
      }
    } // end DT_PEAK
    // else if (current_layer.type == LayerData::DT_CHROMATOGRAM)
  }

  // Helper function for text formatting
  String TOPPViewIdentificationViewBehavior::n_times(Size n, String input)
  {
    String result;
    for (Size i = 0; i < n; ++i)
    {
      result.append(input);
    }
    return result;
  }

  // Helper function that collapses a vector of strings into one string
  String TOPPViewIdentificationViewBehavior::collapseStringVector(vector<String> strings)
  {
    String result;
    for (Size i = 0; i < strings.size(); ++i)
    {
      result.append(strings[i]);
    }
    return result;
  }

  // Helper function that turns fragment annotations into coverage strings for visualization with the sequence
  void TOPPViewIdentificationViewBehavior::extractCoverageStrings(vector<PeptideHit::PeakAnnotation> frag_annotations, String& alpha_string, String& beta_string, Size alpha_size, Size beta_size)
  {
    vector<String> alpha_strings(alpha_size, " ");
    vector<String> beta_strings(beta_size, " ");
    // vectors to keep track of assigned symbols, 0 = nothing, -1 = left, 1 = right, 2 = both
    vector<int> alpha_direction(alpha_size, 0);
    vector<int> beta_direction(beta_size, 0);

    for (Size i = 0; i < frag_annotations.size(); ++i)
    {
      bool has_alpha = frag_annotations[i].annotation.hasSubstring(String("alpha|"));
      bool has_beta = frag_annotations[i].annotation.hasSubstring(String("beta|"));
      // if it has both, it is a complex fragment and more difficult to parse
      // those are ignored for the coverage indicator for now
      if ( has_alpha != has_beta )
      {
        vector<String> dol_split;
        frag_annotations[i].annotation.split("$", dol_split);

        vector<String> bar_split;
        dol_split[0].split("|", bar_split);

        bool alpha = bar_split[0] == "[alpha";
        bool ci = bar_split[1] == "ci";

        vector<String> loss_split;
        dol_split[1].split("-", loss_split);
        String pos_string = loss_split[0].suffix(loss_split[0].size()-1);
        int pos;
        if (pos_string.hasSubstring("]"))
        {
          pos = pos_string.prefix(pos_string.size()-1).toInt()-1;
        }
        else
        {
          pos = pos_string.toInt()-1;
        }

        String frag_type = dol_split[1][0];
        //bool left = (frag_type == "a" || frag_type == "b" || frag_type == "c");
        int direction;
        if (frag_type == "a" || frag_type == "b" || frag_type == "c")
        {
          direction = -1;
        }
        else
        {
          direction = 1;
        }

        if (direction == 1)
        {
          if (alpha)
          {
            pos = alpha_size - pos - 1;
          }
          else
          {
            pos = beta_size - pos - 1;
          }
        }

        String arrow;
        if (ci)
        {
          arrow += "<font color=\"green\">";
        }
        else
        {
          arrow += "<font color=\"red\">";
        }

        if (direction == -1)
        {
          arrow += "&#8636;</font>";
        }
        else
        {
          arrow += "&#8641;</font>";
        }

        if (alpha)
        {
          if (alpha_direction[pos] == 0) // no arrow assigned yet
          {
            alpha_strings[pos] = arrow;
            alpha_direction[pos] = direction;
          }
          else if (alpha_direction[pos] != direction && alpha_direction[pos] != 2) // assigned arrow has different direction, make bidirectional arrow
          {
            alpha_strings[pos] = String("<font color=\"blue\">&#8651;</font>");
            alpha_direction[pos] = 2;
          } // otherwise an arrow with the correct direction is already assigned
        }
        else
        {
          if (beta_direction[pos] == 0) // no arrow assigned yet
          {
            beta_strings[pos] = arrow;
            beta_direction[pos] = direction;
          }
          else if (beta_direction[pos] != direction && beta_direction[pos] != 2) // assigned arrow has different direction, make bidirectional arrow
          {
            beta_strings[pos] = String("<font color=\"blue\">&#8651;</font>");
            beta_direction[pos] = 2;
          } // otherwise an arrow with the correct direction is already assigned
        }
      }
    }
    alpha_string = "<font style=\"\">" + collapseStringVector(alpha_strings) + "</font>";
    beta_string = collapseStringVector(beta_strings);
  }

  void TOPPViewIdentificationViewBehavior::addPrecursorLabels1D_(const vector<Precursor>& pcs)
  {
    LayerData& current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();

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
    else if (current_layer.type == LayerData::DT_CHROMATOGRAM)
    {

    }
  }

  void TOPPViewIdentificationViewBehavior::removeTemporaryAnnotations_(Size spectrum_index)
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

  void TOPPViewIdentificationViewBehavior::addTheoreticalSpectrumLayer_(const PeptideHit& ph)
  {
    SpectrumCanvas* current_canvas = tv_->getActive1DWidget()->canvas();
    LayerData& current_layer = current_canvas->getCurrentLayer();
    SpectrumType& current_spectrum = current_layer.getCurrentSpectrum();

    AASequence aa_sequence = ph.getSequence();

    // get measured spectrum indices and spectrum
    Size current_spectrum_layer_index = current_canvas->activeLayerIndex();
    Size current_spectrum_index = current_layer.getCurrentSpectrumIndex();

    const Param& tv_params = tv_->getParameters();

    PeakSpectrum spectrum;
    TheoreticalSpectrumGenerator generator;
    Param p;
    p.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");

    // these two are true by default, initialize to false here and set to true in the loop below
    p.setValue("add_y_ions", "false", "Add peaks of y-ions to the spectrum");
    p.setValue("add_b_ions", "false", "Add peaks of b-ions to the spectrum");

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

    p.setValue("add_a_ions", tv_params.getValue("preferences:idview:show_a_ions"), "Add peaks of a-ions to the spectrum");
    p.setValue("add_b_ions", tv_params.getValue("preferences:idview:show_b_ions"), "Add peaks of b-ions to the spectrum");
    p.setValue("add_c_ions", tv_params.getValue("preferences:idview:show_c_ions"), "Add peaks of c-ions to the spectrum");
    p.setValue("add_x_ions", tv_params.getValue("preferences:idview:show_x_ions"), "Add peaks of x-ions to the spectrum");
    p.setValue("add_y_ions", tv_params.getValue("preferences:idview:show_y_ions"), "Add peaks of y-ions to the spectrum");
    p.setValue("add_z_ions", tv_params.getValue("preferences:idview:show_z_ions"), "Add peaks of z-ions to the spectrum");
    p.setValue("add_precursor_peaks", tv_params.getValue("preferences:idview:show_precursor"), "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");

    try
    {
      Int max_charge = max(1, ph.getCharge()); // at least generate charge 1 if no charge (0) is annotated

      // generate mass ladder for all charge states
      generator.setParameters(p);
      generator.getSpectrum(spectrum, aa_sequence, 1, max_charge);

    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::warning(tv_, "Error", QString("Spectrum generation failed! (") + e.what() + "). Please report this to the developers (specify what input you used)!");
      return;
    }

    PeakMap new_exp;
    new_exp.addSpectrum(spectrum);
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
    if (current_spectrum_layer_index != theoretical_spectrum_layer_index && !spectrum.getStringDataArrays().empty())
    {
      // Ensure theoretical spectrum is drawn as dashed sticks
      tv_->setDrawMode1D(Spectrum1DCanvas::DM_PEAKS);
      tv_->getActive1DWidget()->canvas()->setCurrentLayerPeakPenStyle(Qt::DashLine);

      // Add ion names as annotations to the theoretical spectrum
      PeakSpectrum::StringDataArray sa = spectrum.getStringDataArrays()[0];

      for (Size i = 0; i != spectrum.size(); ++i)
      {
        DPosition<2> position = DPosition<2>(spectrum[i].getMZ(), spectrum[i].getIntensity());
        QString s(sa[i].c_str());

        if (s.at(0) == 'y')
        {
          Annotation1DItem* item = new Annotation1DPeakItem(position, s, Qt::darkRed);
          item->setSelected(false);
          tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
        }
        else if (s.at(0) == 'b')
        {
          Annotation1DItem* item = new Annotation1DPeakItem(position, s, Qt::darkGreen);
          item->setSelected(false);
          tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
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
        QString s(sa[aligned_peak_indices[i].second].c_str());
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
            const Residue& r = aa_sequence.getResidue(j);
            aa_ss.append(r.getOneLetterCode().toQString());
            if (r.isModified())
            {
              aa_ss.append("*");
            }
          }
          s.append(aa_ss);
          Annotation1DItem* item = tv_->getActive1DWidget()->canvas()->addPeakAnnotation(pi, s, Qt::darkRed);
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
          Annotation1DItem* item = tv_->getActive1DWidget()->canvas()->addPeakAnnotation(pi, s, Qt::darkGreen);
          // save label for later removal
          temporary_annotations_.push_back(item);
        }
        else
        {
          s.append("\n");
          Annotation1DItem* item = tv_->getActive1DWidget()->canvas()->addPeakAnnotation(pi, s, Qt::black);
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
    // Retrieve active 1D widget
    Spectrum1DWidget* widget_1D = tv_->getActive1DWidget();

    // Return if none present
    if (widget_1D == nullptr) return;

    LayerData& current_layer = widget_1D->canvas()->getCurrentLayer();

    // Return if no valid peak layer attached
    if (current_layer.getPeakData()->size() == 0 || current_layer.type != LayerData::DT_PEAK) { return; }

    MSSpectrum& spectrum = (*current_layer.getPeakData())[spectrum_index];
    int ms_level = spectrum.getMSLevel();

    if (ms_level == 2)
    {
      // synchronize PeptideHits with the annotations in the spectrum
      current_layer.synchronizePeakAnnotations();

      // remove all graphical peak annotations as these will be recreated from the stored peak annotations
      Annotations1DContainer& las = current_layer.getAnnotations(spectrum_index);
      auto new_end = std::remove_if(las.begin(), las.end(),
                              [](const Annotation1DItem* a)
                              { return dynamic_cast<const Annotation1DPeakItem*>(a) != nullptr; });
      las.erase(new_end, las.end());

      removeTheoreticalSpectrumLayer_();
    }

    removeTemporaryAnnotations_(spectrum_index);

    // reset selected id indices
    current_layer.peptide_id_index = -1;
    current_layer.peptide_hit_index = -1;

    widget_1D->canvas()->setTextBox(QString());
  }

  void TOPPViewIdentificationViewBehavior::addAnnotationsSpectrumLayer_(const PeptideHit& hit, bool align)
  {
    const vector<PeptideHit::PeakAnnotation>& annotations =
      hit.getPeakAnnotations();
    String seq = hit.getSequence().toString();
    if (seq.empty()) seq = hit.getMetaValue("label");

    SpectrumCanvas* current_canvas = tv_->getActive1DWidget()->canvas();
    LayerData& current_layer = current_canvas->getCurrentLayer();
    Size current_spectrum_layer_index = current_canvas->activeLayerIndex();
    Size current_spectrum_index = current_layer.getCurrentSpectrumIndex();

    const MSSpectrum& current_spectrum = current_layer.getCurrentSpectrum();
    if (align)
    {
      if (current_spectrum.empty())
      {
        LOG_WARN << "Spectrum is empty! Nothing to annotate!" << std::endl;
      }
      else if (!current_spectrum.isSorted())
      {
        QMessageBox::warning(tv_, "Error", "The spectrum is not sorted! Aborting!"); // @TODO: improve error message
        return;
      }
    }

    MSSpectrum ann_spectrum;
    vector<String> labels;
    for (const auto& ann : annotations) // NOLINT
    {
      Peak1D peak(ann.mz, ann.intensity);
      if (align) // align to the measured spectrum
      {
        // @TODO: avoid magic constant (m/z tolerance)
        Int peak_idx = current_spectrum.findNearest(ann.mz, 1e-2);
        if (peak_idx == -1) // no match
        {
          LOG_WARN << "Annotation present for missing peak. m/z: " << ann.mz
                   << endl;
          continue;
        }
        peak = current_spectrum[peak_idx];
      }
      ann_spectrum.push_back(peak);

      String label = ann.annotation;
      // write out positive and negative charges with the correct sign at the end of the annotation string
      switch (ann.charge)
      {
      case 0: break;
      case 1: label += "+"; break;
      case 2: label += "++"; break;
      case -1: label += "-"; break;
      case -2: label += "--"; break;
      default: label += ((ann.charge > 0) ? "+" : "") + String(ann.charge);
      }
      labels.push_back(label);
    }
    ann_spectrum.sortByPosition();

    if (ann_spectrum.getMaxInt() <= 1.0) // undo scaling of intensities
    {
      double max_int = current_layer.getCurrentSpectrum().getMaxInt();
      for (auto& peak : ann_spectrum)
      {
        peak.setIntensity(peak.getIntensity() * max_int);
      }
    }

    PeakMap new_exp;
    new_exp.addSpectrum(ann_spectrum);
    ExperimentSharedPtrType new_exp_sptr(new PeakMap(new_exp));
    FeatureMapSharedPtrType f_dummy(new FeatureMapType());
    ConsensusMapSharedPtrType c_dummy(new ConsensusMapType());
    vector<PeptideIdentification> p_dummy;

    // Block update events for identification widget
    tv_->getSpectraIdentificationViewWidget()->ignore_update = true;

    String layer_caption = seq + " (identification view)";
    tv_->addData(f_dummy, c_dummy, p_dummy, new_exp_sptr, LayerData::DT_PEAK, true, false, false, "", layer_caption);

    // get layer index of new layer
    Size theoretical_spectrum_layer_index = tv_->getActive1DWidget()->canvas()->activeLayerIndex();

    // kind of a hack to check whether adding the layer was successful
    if (current_spectrum_layer_index != theoretical_spectrum_layer_index)
    {
      // Ensure theoretical spectrum is drawn as sticks
      tv_->setDrawMode1D(Spectrum1DCanvas::DM_PEAKS);
      // ensure intensities are on the same scale as the measured spectrum:
      tv_->setIntensityMode(SpectrumCanvas::IM_SNAP);

      // Add ion names to the annotations spectrum
      for (Size i = 0; i != ann_spectrum.size(); ++i)
      {
        DPosition<2> position(ann_spectrum[i].getMZ(),
                              ann_spectrum[i].getIntensity());
        const String& label = labels[i];
        QColor color;
        // XL-MS specific coloring of the labels, green for linear fragments and red for cross-linked fragments
        if (label.hasSubstring("[alpha|") || label.hasSubstring("[beta|"))
        {
          if (label.hasSubstring("|ci$"))
          {
            color = Qt::darkGreen;
          }
          else if (label.hasSubstring("|xi$"))
          {
            color = Qt::darkRed;
          }
        }
        else // different colors for left/right fragments (e.g. b/y ions)
        {
          color = (label.at(0) < 'n') ? Qt::darkRed : Qt::darkGreen;
        }

        Annotation1DItem* item = new Annotation1DPeakItem(position, label.toQString(), color);
        item->setSelected(false);
        tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
      }

      tv_->getActive1DWidget()->canvas()->activateLayer(current_spectrum_layer_index);
      tv_->getActive1DWidget()->canvas()->getCurrentLayer().setCurrentSpectrumIndex(current_spectrum_index);

      // zoom visible area to real data range:
      DRange<2> visible_area = tv_->getActive1DWidget()->canvas()->getVisibleArea();
      double min_mz = tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentSpectrum().getMin()[0];
      double max_mz = tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentSpectrum().getMax()[0];
      double delta_mz = max_mz - min_mz;
      visible_area.setMin(min_mz - 0.1 * delta_mz);
      visible_area.setMax(max_mz + 0.1 * delta_mz);
      tv_->getActive1DWidget()->canvas()->setVisibleArea(visible_area);

      tv_->updateLayerBar();
      tv_->getSpectraIdentificationViewWidget()->ignore_update = false;
    }
  }


  void TOPPViewIdentificationViewBehavior::removeTheoreticalSpectrumLayer_()
  {
    Spectrum1DWidget* spectrum_widget_1D = tv_->getActive1DWidget();
    if (spectrum_widget_1D)
    {
      Spectrum1DCanvas* canvas_1D = spectrum_widget_1D->canvas();

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
    if (w == nullptr) return;

    SpectrumCanvas* current_canvas = w->canvas();
    LayerData& current_layer = current_canvas->getCurrentLayer();
    SpectrumType& current_spectrum = current_layer.getCurrentSpectrum();

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
    Spectrum1DWidget* widget_1D = tv_->getActive1DWidget();

    // return if no active 1D widget is present
    if (widget_1D == nullptr) return;

    // clear textbox
    widget_1D->canvas()->setTextBox(QString());

    // remove precusor labels, theoretical spectra and trigger repaint
    LayerData& cl = tv_->getActive1DWidget()->canvas()->getCurrentLayer();
    removeTemporaryAnnotations_(cl.getCurrentSpectrumIndex());
    removeTheoreticalSpectrumLayer_();
    cl.peptide_id_index = -1;
    cl.peptide_hit_index = -1;
    tv_->getActive1DWidget()->canvas()->repaint();
  }

  void TOPPViewIdentificationViewBehavior::setVisibleArea1D(double l, double h)
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
