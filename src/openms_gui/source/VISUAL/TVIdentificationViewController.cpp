// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TVIdentificationViewController.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DCaret.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/SpectraIDViewTab.h>


#include <boost/range/adaptor/reversed.hpp>
#include <QtWidgets/QMessageBox>
#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

//#define DEBUG_IDENTIFICATION_VIEW

namespace OpenMS
{
  TVIdentificationViewController::TVIdentificationViewController(TOPPViewBase* parent, SpectraIDViewTab* spec_id_view) :
    TVControllerBase(parent),
    spec_id_view_(spec_id_view)
  {
  }

  void TVIdentificationViewController::showSpectrumAsNew1D(int index)
  {
    // Show spectrum "index" without selecting an identification
    showSpectrumAsNew1D(index, -1, -1);
  }

  void TVIdentificationViewController::showSpectrumAsNew1D(int spectrum_index, int peptide_id_index, int peptide_hit_index)
  {
    // basic behavior 1
    auto& layer = tv_->getActiveCanvas()->getCurrentLayer();

    if (layer.type == LayerDataBase::DT_PEAK)
    {
      // open new 1D widget with the current default parameters
      auto* w = new Plot1DWidget(tv_->getCanvasParameters(1), DIM::Y, (QWidget*)tv_->getWorkspace());

      // copy data from current layer (keeps the TYPE and underlying data identical)
      auto new_1d = layer.to1DLayer();
      bool index_ok = new_1d->hasIndex(spectrum_index);
      if (!index_ok || !w->canvas()->addLayer(std::move(new_1d)))
      {
        // Behavior if its neither (user may have clicked on an empty tree or a
        // dummy entry as drawn by SpectraTreeTab::updateEntries)
        QMessageBox::critical(w, "Error", "Cannot open data. Aborting!");
        return;
      }

      w->canvas()->activateSpectrum(spectrum_index);

      // set relative (%) view of visible area
      w->canvas()->setIntensityMode(PlotCanvas::IM_SNAP);

      // set visible area to visible area in 2D view
      w->canvas()->setVisibleArea(tv_->getActiveCanvas()->getVisibleArea());
      
      w->canvas()->getCurrentLayer().setName(layer.getName());
      w->canvas()->getCurrentLayer().setNameSuffix(layer.getNameSuffix());

      tv_->showPlotWidgetInWindow(w);

      ///////////////////////////////////////////////////////////////////////////////
      // Visualization of ID data

      // if no peptide identification or peptide hit index provided we can return now
      if (peptide_id_index == -1 || peptide_hit_index == -1)
      { 
        return;
      }

      // get peptide identification
      auto layer_1d_peak = dynamic_cast<const LayerData1DPeak*>(&w->canvas()->getCurrentLayer());
      const vector<PeptideIdentification>& pis = layer_1d_peak->getCurrentSpectrum().getPeptideIdentifications();

      if (!pis.empty())
      {
        switch (layer_1d_peak->getCurrentSpectrum().getMSLevel())
        {
          // mass fingerprint annotation of name etc.
          case 1:
          { 
            addPeakAnnotations_(pis); 
            break;
          }

          // annotation with stored fragments or synthesized theoretical spectrum
          case 2:
          {
            // check if index in bounds and hits are present
            if (peptide_id_index < static_cast<int>(pis.size())
              && peptide_hit_index < static_cast<int>(pis[peptide_id_index].getHits().size()))
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
                addPeakAnnotationsFromID_(ph);
              }
            }
            break;
          }
          default:
            OPENMS_LOG_WARN << "Annotation of MS level > 2 not supported.!" << endl;
        }
      }

      // TODO Why would this need to trigger an update in e.g. the Tab Views??
      tv_->updateLayerBar(); // todo replace
      tv_->updateViewBar();
      tv_->updateFilterBar();
      tv_->updateMenu();
    }
    // else if (layer.type == LayerDataBase::DT_CHROMATOGRAM)
  }

  void TVIdentificationViewController::addPeakAnnotations_(const vector<PeptideIdentification>& ph)
  {
    // called anew for every click on a spectrum
    auto getCurrentLayer = [&]() -> LayerData1DPeak& { return dynamic_cast<LayerData1DPeak&>(tv_->getActive1DWidget()->canvas()->getCurrentLayer()); };

    if (getCurrentLayer().getCurrentSpectrum().empty())
    {
      OPENMS_LOG_WARN << "Spectrum is empty! Nothing to annotate!" << endl;
    }

    // mass precision to match a peak's m/z to a feature m/z
    // m/z values of features are usually an average over multiple scans...
    constexpr double ppm = 0.5;

    array<QColor, 5> cols{ Qt::blue, Qt::green, Qt::red, Qt::gray, Qt::darkYellow };

    if (!getCurrentLayer().getCurrentSpectrum().isSorted())
    {
      QMessageBox::warning(tv_, "Error", "The spectrum is not sorted! Aborting!");
      return;
    }
    for (vector<PeptideIdentification>::const_iterator it = ph.begin();
                                                            it!= ph.end();
                                                            ++it)
    {
      if (!it->hasMZ())
      {
        continue;
      }
      double mz = it->getMZ();
      Size peak_idx = getCurrentLayer().getCurrentSpectrum().findNearest(mz);

      // m/z fits ?
      if (Math::getPPMAbs(mz, getCurrentLayer().getCurrentSpectrum()[peak_idx].getMZ()) > ppm)
      {
        continue;
      }

      double peak_int = getCurrentLayer().getCurrentSpectrum()[peak_idx].getIntensity();

      Annotation1DCaret<Peak1D>* first_dit(nullptr);
      // we could have many hits for different compounds which have the exact same sum formula... so first group by sum formula
      map<String, StringList> formula_to_names;
      for (const PeptideHit& pep : it->getHits())
      {
        if (pep.metaValueExists("identifier") && pep.metaValueExists("chemical_formula"))
        {
          String name = pep.getMetaValue("identifier");
          if (name.length() > 20)
          {
            name = name.substr(0, 17) + "...";
          }
          String cf = pep.getMetaValue("chemical_formula");
          if (cf.empty())
          {
            continue; // skip unannotated "null" peaks
          }
          formula_to_names[cf].push_back(name);
        }
        else
        {
          StringList msg;
          if (!pep.metaValueExists("identifier"))
          {
            msg.push_back("identifier");
          }
          if (!pep.metaValueExists("chemical_formula"))
          {
            msg.push_back("chemical_formula");
          }
          OPENMS_LOG_WARN << "Missing meta-value(s): " << ListUtils::concatenate(msg, ", ") << ". Cannot annotate!\n";
        }
      }

      // assemble annotation (each formula gets a paragraph)
      String text = "<html><body>";
      Size i = 0;
      for (map<String, StringList>::iterator ith = formula_to_names.begin();
           ith!= formula_to_names.end(); ++ith)
      {
        if (++i == cols.size())
        { // at this point, this is the 4th entry.. which we don't show any more...
          text += String("<b><span style=\"color:") + cols[i].name() + "\">..." + Size(distance(formula_to_names.begin(), formula_to_names.end()) - 4 + 1) + " more</span></b><br>";
          break;
        }
        text += String("<b><span style=\"color:") + cols[i].name() + "\">" + ith->first + "</span></b><br>\n";
        // carets for isotope profile
        EmpiricalFormula ef(ith->first);
        IsotopeDistribution id = ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(3)); // three isotopes at most
        double int_factor = peak_int / id.begin()->getIntensity();
        Annotation1DCaret<Peak1D>::PositionsType points;
        Size itic(0);
        for (const Peak1D& iso : id)
        {
          points.push_back(Annotation1DCaret<Peak1D>::PointType(mz + itic * Constants::C13C12_MASSDIFF_U, iso.getIntensity() * int_factor));
          ++itic;
        }
        auto ditem = new Annotation1DCaret<Peak1D>(points,
                                                   QString(),
                                                   cols[i],
                                                   String(getCurrentLayer().param.getValue("peak_color").toString()).toQString());
        ditem->setSelected(false);
        temporary_annotations_.push_back(ditem); // for removal (no ownership)
        getCurrentLayer().getCurrentAnnotations().push_front(ditem); // for visualization (ownership)
        if (first_dit == nullptr) 
        {
          first_dit = ditem; // remember first item (we append the text, when ready)
        }

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

  void TVIdentificationViewController::activate1DSpectrum(int index)
  {
    activate1DSpectrum(index, -1, -1);
  }

  void TVIdentificationViewController::activate1DSpectrum(
    int spectrum_index,
    int peptide_id_index,
    int peptide_hit_index)
  {
    Plot1DWidget* widget_1D = tv_->getActive1DWidget();

    // if no active 1D widget is present
    if (widget_1D == nullptr) 
    { // ... create one
      showSpectrumAsNew1D(spectrum_index, peptide_id_index, peptide_hit_index);
      return;
    }

    // lambda which returns the current layer. This has to be used throughout this function to ensure
    // being up-to-date (no invalidated pointer etc.)
    // even after adding a layer with e.g. addTheoreticalSpectrumLayer_ in L372.
    auto current_layer = [&]() -> LayerData1DPeak& { return dynamic_cast<LayerData1DPeak&>(tv_->getActive1DWidget()->canvas()->getCurrentLayer()); };

    widget_1D->canvas()->activateSpectrum(spectrum_index);
    current_layer().peptide_id_index = peptide_id_index;
    current_layer().peptide_hit_index = peptide_hit_index;

    if (current_layer().type == LayerDataBase::DT_PEAK)
    {
      UInt ms_level = current_layer().getCurrentSpectrum().getMSLevel();

      const vector<PeptideIdentification>& pis = current_layer().getCurrentSpectrum().getPeptideIdentifications();
      switch (ms_level)
      {
        case 1: // mass fingerprint annotation of name etc and precursor labels
        {
          addPeakAnnotations_(pis);
          vector<Precursor> precursors;

          // collect all MS2 spectra precursor till next MS1 spectrum is encountered
          for (Size i = spectrum_index + 1; i < current_layer().getPeakData()->size(); ++i)
          {
            if ((*current_layer().getPeakData())[i].getMSLevel() == 1)
            {
              break;
            }
            // skip MS2 without precursor
            if ((*current_layer().getPeakData())[i].getPrecursors().empty())
            {
              continue;
            }
            // there should be only one precursor per MS2 spectrum.
            vector<Precursor> pcs = (*current_layer().getPeakData())[i].getPrecursors();
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

              // synchronize PeptideHits with the annotations in the spectrum
              current_layer().synchronizePeakAnnotations();
              // remove labels and theoretical spectrum (will be recreated using PH annotations)
              removeGraphicalPeakAnnotations_(spectrum_index);
              removeTheoreticalSpectrumLayer_();

              // return if no active 1D widget is present
              if (widget_1D == nullptr)
              { 
                return;
              }
              // update current PeptideHit with the synchronized one
              widget_1D->canvas()->activateSpectrum(spectrum_index);
              const vector<PeptideIdentification>& pis2 = current_layer().getCurrentSpectrum().getPeptideIdentifications();
              ph = pis2[peptide_id_index].getHits()[peptide_hit_index];

            }
            // use stored fragment annotations
            addPeakAnnotationsFromID_(ph);

            if (ph.metaValueExists(Constants::UserParam::OPENPEPXL_XL_TYPE)) // if this meta value exists, this should be an XL-MS annotation
            {
              String box_text;
              String vert_bar = "&#124;";

              if (ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "loop-link")
              {
                String hor_bar = "_";
                String seq_alpha = ph.getSequence().toUnmodifiedString();
                int xl_pos_alpha = String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1)).toInt();
                int xl_pos_beta = String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2)).toInt() - xl_pos_alpha - 1;

                String alpha_cov;
                String beta_cov;
                extractCoverageStrings(ph.getPeakAnnotations(), alpha_cov, beta_cov, seq_alpha.size(), 0);

                // String formatting
                box_text += alpha_cov + "<br>" +  seq_alpha +  "<br>" + String(xl_pos_alpha, ' ') +  vert_bar + n_times(xl_pos_beta, hor_bar) + vert_bar;
                // cut out line: "<br>" + String(xl_pos_alpha, ' ') + vert_bar + String(xl_pos_beta, ' ') + vert_bar +
              }
              else if (ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "cross-link")
              {
                String seq_alpha = ph.getSequence().toUnmodifiedString();
                String seq_beta = AASequence::fromString(ph.getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE)).toUnmodifiedString();
                int xl_pos_alpha = String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1)).toInt();
                int xl_pos_beta = String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2)).toInt();

                // String formatting
                Size prefix_length = max(xl_pos_alpha, xl_pos_beta);
                //Size suffix_length = max(seq_alpha.size() - xl_pos_alpha, seq_beta.size() - xl_pos_beta);
                Size alpha_space = prefix_length - xl_pos_alpha;
                Size beta_space = prefix_length - xl_pos_beta;

                String alpha_cov;
                String beta_cov;
                extractCoverageStrings(ph.getPeakAnnotations(), alpha_cov, beta_cov, seq_alpha.size(), seq_beta.size());

                box_text += String(alpha_space, ' ') + alpha_cov + "<br>" + String(alpha_space, ' ') + seq_alpha + "<br>" + String(prefix_length, ' ') + vert_bar + "<br>" + String(beta_space, ' ') + seq_beta + "<br>" + String(beta_space, ' ') + beta_cov;
                // color: <font color=\"green\">&boxur;</font>
              }
              else // if (ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE) == "mono-link")
              {
                String seq_alpha = ph.getSequence().toUnmodifiedString();
                int xl_pos_alpha = String(ph.getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1)).toInt();
                Size prefix_length = xl_pos_alpha;

                String alpha_cov;
                String beta_cov;
                extractCoverageStrings(ph.getPeakAnnotations(), alpha_cov, beta_cov, seq_alpha.size(), 0);

                box_text += alpha_cov + "<br>" + seq_alpha + "<br>" + String(prefix_length, ' ') + vert_bar;

              }
              box_text = R"(<font size="5" style="background-color:white;"><pre>)" + box_text + "</pre></font> ";
              widget_1D->canvas()->setTextBox(box_text.toQString());
            }
            else if (ph.getPeakAnnotations().empty()) // only write the sequence
            {
              String seq = ph.getSequence().toString();
              if (seq.empty())
              {
                seq = ph.getMetaValue("label"); // e.g. for RNA sequences
              }
              widget_1D->canvas()->setTextBox(seq.toQString());
            }
            else if (widget_1D->canvas()->isIonLadderVisible())
            {
              if (!ph.getSequence().empty()) // generate sequence diagram for a peptide
              {
                // @TODO: read ion list from the input file (meta value)
                static vector<String> top_ions = ListUtils::create<String>("a,b,c");
                static vector<String> bottom_ions = ListUtils::create<String>("x,y,z");
                String diagram = generateSequenceDiagram_(
                  ph.getSequence(),
                  ph.getPeakAnnotations(),
                  top_ions,
                  bottom_ions);
                widget_1D->canvas()->setTextBox(diagram.toQString());
              }
              else if (ph.metaValueExists("label")) // generate sequence diagram for RNA
              {
                try
                {
                  // @TODO: read ion list from the input file (meta value)
                  NASequence na_seq = NASequence::fromString(ph.getMetaValue("label"));
                  static vector<String> top_ions = ListUtils::create<String>("a-B,a,b,c,d");
                  static vector<String> bottom_ions = ListUtils::create<String>("w,x,y,z");
                  String diagram = generateSequenceDiagram_(na_seq, ph.getPeakAnnotations(),
                                                            top_ions, bottom_ions);
                  widget_1D->canvas()->setTextBox(diagram.toQString());
                }
                catch (Exception::ParseError&) // label doesn't contain have a valid seq.
                {
                }
              }
            }
          }
          break;
        }
        default:
          OPENMS_LOG_WARN << "Annotation of MS level > 2 not supported." << endl;
      }
    } // end DT_PEAK
    // else if (current_layer().type == LayerDataBase::DT_CHROMATOGRAM)
  }

  // Helper function for text formatting
  String TVIdentificationViewController::n_times(Size n, const String& input)
  {
    String result;
    for (Size i = 0; i < n; ++i)
    {
      result.append(input);
    }
    return result;
  }

  // Helper function that collapses a vector of strings into one string
  String TVIdentificationViewController::collapseStringVector(vector<String> strings)
  {
    String result;
    for (Size i = 0; i < strings.size(); ++i)
    {
      result.append(strings[i]);
    }
    return result;
  }

  // Helper function that turns fragment annotations into coverage strings for visualization with the sequence
  void TVIdentificationViewController::extractCoverageStrings(vector<PeptideHit::PeakAnnotation> frag_annotations, String& alpha_string, String& beta_string, Size alpha_size, Size beta_size)
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
        // remove b / y ion type letter (must be at first position of second string after $-split)
        String pos_string = loss_split[0].suffix(loss_split[0].size()-1);
        int pos;
        if (pos_string.hasSubstring("]")) // this means the loss_split with "-" did not split the string
        {
          // remove the "]" and possible charges at its right side
          vector<String> pos_split;
          pos_string.split("]", pos_split);
          pos_string = pos_split[0];
          pos = pos_string.toInt()-1;
        }
        else // loss was found and splitted, so the remaining string is just the position
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


  void TVIdentificationViewController::generateSequenceRow_(const AASequence& seq, vector<String>& row)
  {
    // @TODO: spell out modifications or just use an indicator like "*"?
    // @TODO: support "user defined modifications"?
    if (seq.hasNTerminalModification())
    {
      row[0] =  + "." + seq.getNTerminalModificationName();
    }
    Size col_index = 1;
    for (const auto& aa : seq)
    {
      row[col_index] = "<b>" + aa.getOneLetterCode();
      if (aa.isModified())
      {
        row[col_index] += "(" + aa.getModificationName() + ")";
      }
      row[col_index] += "</b>";
      col_index += 2;
    }
    if (seq.hasCTerminalModification())
    {
      row[row.size() - 1] = "." + seq.getCTerminalModificationName();
    }
  }


  void TVIdentificationViewController::generateSequenceRow_(const NASequence& seq, vector<String>& row)
  {
    if (seq.hasFivePrimeMod())
    {
      const String& code = seq.getFivePrimeMod()->getCode();
      row[0] = (code == "5'-p" ? "p" : code);
    }
    Size col_index = 1;
    for (const auto& ribo : seq)
    {
      row[col_index] = "<b>" + ribo.getCode() + "</b>";
      col_index += 2;
    }
    if (seq.hasThreePrimeMod())
    {
      const String& code = seq.getThreePrimeMod()->getCode();
      row[row.size() - 1] = (code == "3'-p" ? "p" : code);
    }
  }


  template <typename SeqType>
  String TVIdentificationViewController::generateSequenceDiagram_(
    const SeqType& seq,
    const vector<PeptideHit::PeakAnnotation>& annotations,
    const vector<String>& top_ions,
    const vector<String>& bottom_ions)
  {
    #ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "Generating Sequence Diagram: " << endl;
    #endif
    map<String, set<Size>> ion_pos;
    for (const auto& ann : annotations)
    {
      const String& label = ann.annotation;
      #ifdef DEBUG_IDENTIFICATION_VIEW
        cout << "Adding Peak Annotation to Diagram: " << label << endl;
      #endif
      // expected format: [ion][number][...]
      if ((label.size() < 2) || !islower(label[0]) || !isdigit(label[1]))
      {
        continue;
      }
      // cut out the position number:
      Size split = label.find_first_not_of("0123456789", 2);
      String ion = label.prefix(1);
      // special case for RNA: "a[n]-B", where "[n]" is the ion number
      // -> don't forget to add the "-B" back on if it's there:
      String more_ion = label.substr(split);
      if (more_ion == "-B")
      {
        ion += more_ion;
      }
      Size pos = label.substr(1, split - 1).toInt();
      ion_pos[ion].insert(pos);
      #ifdef DEBUG_IDENTIFICATION_VIEW
        cout << "Ion: " << ion << " pos: " << pos << endl;
      #endif
    }

    vector<vector<String>> table; // vector of rows
    table.resize(top_ions.size() + bottom_ions.size() + 3);
    Size n_cols = seq.size() * 2 + 1;
    for (auto& row : table)
    {
      row.resize(n_cols);
    }
    if (!top_ions.empty())
    {
      for (Size i = 1; i < seq.size(); ++i)
      {
        // @TODO: check spacing for i > 9
        table[0][i * 2] = "<small>" + String(i) + "</small>";
      }
    }
    Size row_index = 1;
    // ion annotations above sequence - reverse order to have first ion closest to sequence:
    for (const String& ion : boost::adaptors::reverse(top_ions))
    {
      table[row_index][0] = "<small>" + ion + "</small>";
      for (Size pos : ion_pos[ion])
      {
        #ifdef DEBUG_IDENTIFICATION_VIEW
          cout << "Found ion: " << ion << " pos: " << pos << endl;
        #endif

        Size col_index = 2 * pos;
        if ((row_index == 1) || (table[row_index - 1][col_index].empty()))
        {
          table[row_index][col_index] = "&#9488;"; // box drawing: down and left
        }
        else
        {
          table[row_index][col_index] = "&#9508;"; // box drawing: vertical and left
        }
        table[row_index][col_index - 1] = "&#9590;"; // box drawing: right
      }
      if (row_index > 1)
      {
        for (Size col_index = 2; col_index < n_cols - 2; col_index += 2)
        {
          if (table[row_index][col_index].empty() && !table[row_index - 1][col_index].empty())
          {
            table[row_index][col_index] = "&#9474;"; // box drawing: vertical
          }
        }
      }
      ++row_index;
    }
    // sequence itself:
    generateSequenceRow_(seq, table[row_index]);
    // ion annotations below sequence - iterate over the bottom ions in reverse order (bottom-most first):
    row_index = table.size() - 2;
    for (const String& ion : boost::adaptors::reverse(bottom_ions))
    {
      table[row_index][n_cols - 1] = "<small>" + ion + "<small>";
      for (Size pos : ion_pos[ion])
      {
        Size col_index = n_cols - 2 * pos - 1;
        if ((row_index == table.size() - 1) || (table[row_index + 1][col_index].empty()))
        {
          table[row_index][col_index] = "&#9492;"; // box drawing: up and right
        }
        else
        {
          table[row_index][col_index] = "&#9500;"; // box drawing: vertical and right
        }
        table[row_index][col_index + 1] = "&#9588;"; // box drawing: left
      }
      if (row_index < table.size() - 2)
      {
        for (Size col_index = 2; col_index < n_cols - 2; col_index += 2)
        {
          if (table[row_index][col_index].empty() && !table[row_index + 1][col_index].empty())
          {
            table[row_index][col_index] = "&#9474;"; // box drawing: vertical
          }
        }
      }
      --row_index;
    }
    // "row_index" is again at the sequence row - fill in "split indicators":
    for (Size col_index = 2; col_index < n_cols - 2 ; col_index += 2)
    {
      bool top = !top_ions.empty() && !table[row_index - 1][col_index].empty();
      bool bottom = !bottom_ions.empty() && !table[row_index + 1][col_index].empty();
      if (top && bottom)
      {
        table[row_index][col_index] = "&#9474;"; // box drawing: vertical
      }
      else if (top)
      {
        table[row_index][col_index] = "&#9589;"; // box drawing: up
      }
      else if (bottom)
      {
        table[row_index][col_index] = "&#9591;"; // box drawing: down
      }
    }
    if (!bottom_ions.empty())
    {
      for (Size i = 1; i < seq.size(); ++i)
      {
        // @TODO: check spacing in diagram for i > 9
        table[table.size() - 1][n_cols - 2 * i - 1] = "<small>" + String(i) + "</small>";
      }
    }

    String html = "<table cellspacing=\"0\">";
    for (const auto& row : table)
    {
      html += "<tr>";
      for (const String& cell : row)
      {
        #ifdef DEBUG_IDENTIFICATION_VIEW
          cout << "cell: '" << cell << "'" << endl;
        #endif
        html += "<td align=\"center\">" + cell + "</td>";
      }
      html += "</tr>";
    }
    html += "</table>";

    #ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "Generated html:\n" << html << endl;
    #endif
    return html;
  }


  // add specializations to allow template implementation outside of header file:
  template String TVIdentificationViewController::generateSequenceDiagram_<AASequence>(const AASequence& seq, const vector<PeptideHit::PeakAnnotation>& annotations, const StringList& top_ions, const StringList& bottom_ions);
  template String TVIdentificationViewController::generateSequenceDiagram_<NASequence>(const NASequence& seq, const vector<PeptideHit::PeakAnnotation>& annotations, const StringList& top_ions, const StringList& bottom_ions);


  void TVIdentificationViewController::addPrecursorLabels1D_(const vector<Precursor>& pcs)
  {
    auto& current_layer = dynamic_cast<LayerData1DPeak&>(tv_->getActive1DWidget()->canvas()->getCurrentLayer());

    if (current_layer.type == LayerDataBase::DT_PEAK)
    {
      const SpectrumType& spectrum = current_layer.getCurrentSpectrum();

      for (const Precursor& pre : pcs)
      {
        // determine start and stop of isolation window
        double center_mz = pre.metaValueExists("isolation window target m/z") ?
        double(pre.getMetaValue("isolation window target m/z")) : pre.getMZ();
        double isolation_window_lower_mz = center_mz - pre.getIsolationWindowLowerOffset();
        double isolation_window_upper_mz = center_mz + pre.getIsolationWindowUpperOffset();

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

        Annotation1DDistanceItem* item = new Annotation1DDistanceItem(QString::number(pre.getCharge()), lower_position, upper_position);
        // add additional tick at precursor target position (e.g. to show if isolation window is asymmetric)
        vector<PointXYType> ticks;
        ticks.emplace_back(pre.getMZ(), 0);
        item->setTicks(ticks);
        item->setSelected(false);

        temporary_annotations_.push_back(item); // for removal (no ownership)
        current_layer.getCurrentAnnotations().push_front(item); // for visualization (ownership)
      }
    }
    else if (current_layer.type == LayerDataBase::DT_CHROMATOGRAM)
    {

    }
  }

  void TVIdentificationViewController::removeTemporaryAnnotations_(Size spectrum_index)
  {
#ifdef DEBUG_IDENTIFICATION_VIEW
    cout << "removePrecursorLabels1D_ " << spectrum_index << endl;
#endif
    // Delete annotations added by IdentificationView (but not user added annotations)
    auto& current_layer = tv_->getActive1DWidget()->canvas()->getCurrentLayer();
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

  void TVIdentificationViewController::addTheoreticalSpectrumLayer_(const PeptideHit& ph)
  {
    PlotCanvas* current_canvas = tv_->getActive1DWidget()->canvas();
    auto& current_layer = dynamic_cast<LayerData1DPeak&>(current_canvas->getCurrentLayer());
    const SpectrumType& current_spectrum = current_layer.getCurrentSpectrum();

    const AASequence& aa_sequence = ph.getSequence();

    // get measured spectrum indices and spectrum
    Size current_spectrum_layer_index = current_canvas->getCurrentLayerIndex();
    Size current_spectrum_index = current_layer.getCurrentIndex();

    const Param& tv_params = tv_->getParameters();
    Param tag_params = tv_params.copy("preferences:user:idview:tsg:", true);
    // override: enable metavalues for simulated peaks (needed for annotation)
    assert(tag_params.exists("add_metainfo"));
    tag_params.setValue("add_metainfo", "true");

    PeakSpectrum theo_spectrum;
    try
    {
      TheoreticalSpectrumGenerator generator;
      Int max_charge = max(1, ph.getCharge()); // at least generate charge 1 if no charge (0) is annotated

      // generate mass ladder for all charge states
      generator.setParameters(tag_params);
      generator.getSpectrum(theo_spectrum, aa_sequence, 1, max_charge);

      // scale spectrum to maximum peak intensity of real spectrum
      auto scale_by = current_spectrum.getMaxIntensity();
      std::for_each(theo_spectrum.begin(), theo_spectrum.end(), [scale_by](Peak1D& p) { p.setIntensity(p.getIntensity() * scale_by); });
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::warning(tv_, "Error", QString("Spectrum generation failed! (") + e.what() + "). Please report this to the developers (specify what input you used)!");
      return;
    }

    // Block update events for identification widget
    spec_id_view_->ignore_update = true;
    RAIICleanup cleanup([&]() { spec_id_view_->ignore_update = false; });

    PeakMap new_exp;
    new_exp.addSpectrum(theo_spectrum);
    ExperimentSharedPtrType new_exp_sptr(new PeakMap(new_exp));
    LayerDataBase::ODExperimentSharedPtrType od_dummy(new OnDiscMSExperiment());
    String layer_caption = aa_sequence.toString() + " (identification view)";
    current_canvas->addPeakLayer(new_exp_sptr, od_dummy, layer_caption);

    // get layer index of new layer
    Size theoretical_spectrum_layer_index = tv_->getActive1DWidget()->canvas()->getCurrentLayerIndex();

    // kind of a hack to check whether adding the layer was successful
    if (current_spectrum_layer_index != theoretical_spectrum_layer_index && !theo_spectrum.getStringDataArrays().empty())
    {
      // Ensure theoretical spectrum is drawn as dashed sticks
      tv_->setDrawMode1D(Plot1DCanvas::DM_PEAKS);
      tv_->getActive1DWidget()->canvas()->setCurrentLayerPeakPenStyle(Qt::DashLine);

      // Add ion names as annotations to the theoretical spectrum
      PeakSpectrum::StringDataArray sa = theo_spectrum.getStringDataArrays()[0];

      for (Size i = 0; i != theo_spectrum.size(); ++i)
      {
        Peak1D position(theo_spectrum[i].getMZ(), theo_spectrum[i].getIntensity());
        QString s(sa[i].c_str());

        if (s.at(0) == 'y')
        {
          auto item = new Annotation1DPeakItem<Peak1D>(position, s, Qt::darkRed);
          item->setSelected(false);
          tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
        }
        else if (s.at(0) == 'b')
        {
          auto item = new Annotation1DPeakItem<Peak1D>(position, s, Qt::darkGreen);
          item->setSelected(false);
          tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
        }
      }

      // remove theoretical and activate real data layer and spectrum
      tv_->getActive1DWidget()->canvas()->changeVisibility(theoretical_spectrum_layer_index, false);
      tv_->getActive1DWidget()->canvas()->activateLayer(current_spectrum_layer_index);
      tv_->getActive1DWidget()->canvas()->getCurrentLayer().setCurrentIndex(current_spectrum_index);

      // zoom to maximum visible area in real data (as theoretical might be much larger and therefor squeezes the interesting part)
      auto spec_range = tv_->getActive1DWidget()->canvas()->getCurrentLayer().getRange();
      spec_range.scaleBy(1.2);
      tv_->getActive1DWidget()->canvas()->setVisibleArea(RangeAllType().assign(spec_range));

      // spectra alignment
      Param p_align = tv_params.copy("preferences:user:idview:align", true);
      tv_->getActive1DWidget()->performAlignment(current_spectrum_layer_index, theoretical_spectrum_layer_index, p_align);

      vector<pair<Size, Size> > aligned_peak_indices = tv_->getActive1DWidget()->canvas()->getAlignedPeaksIndices();

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
          aa_ss.replace(QRegularExpression("[(].*[)]"), "*");
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
    }
  }

  void TVIdentificationViewController::removeGraphicalPeakAnnotations_(int spectrum_index)
  {
    auto* widget_1D = tv_->getActive1DWidget();
    auto& current_layer = widget_1D->canvas()->getCurrentLayer();

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
                              return dynamic_cast<const Annotation1DPeakItem<Peak1D>*>(a) != nullptr;
                            });
    las.erase(new_end, las.end());

    return;
  }

  void TVIdentificationViewController::deactivate1DSpectrum(int spectrum_index)
  {
    // Retrieve active 1D widget
    Plot1DWidget* widget_1D = tv_->getActive1DWidget();

    // Return if none present
    if (widget_1D == nullptr)
    {
      return;
    }
    auto& current_layer = widget_1D->canvas()->getCurrentLayer();

    // Return if no valid peak layer attached
    auto* current_layer_ptr = dynamic_cast<LayerData1DPeak*>(&current_layer);
    if (!current_layer_ptr || current_layer_ptr->getPeakData()->empty())
    { 
      return;
    }
    MSSpectrum& spectrum = (*current_layer_ptr->getPeakDataMuteable())[spectrum_index];
    int ms_level = spectrum.getMSLevel();
    if (ms_level == 2)
    {
      // synchronize PeptideHits with the annotations in the spectrum
      current_layer_ptr->synchronizePeakAnnotations();
      removeGraphicalPeakAnnotations_(spectrum_index);
      removeTheoreticalSpectrumLayer_();
    }

    removeTemporaryAnnotations_(spectrum_index);

    // reset selected id indices
    current_layer.peptide_id_index = -1;
    current_layer.peptide_hit_index = -1;

    widget_1D->canvas()->setTextBox(QString());
  }

  void TVIdentificationViewController::addPeakAnnotationsFromID_(const PeptideHit& hit)
  {
    // get annotations and sequence
    const vector<PeptideHit::PeakAnnotation>& annotations =
      hit.getPeakAnnotations();

    String seq = hit.getSequence().toString();
    if (seq.empty())
    { // no sequence information stored? use label
      if (hit.metaValueExists("label"))
      {
        seq = hit.getMetaValue("label");
      }
    }

    auto* current_canvas = tv_->getActive1DWidget()->canvas();
    LayerDataBase& current_layer = current_canvas->getCurrentLayer();
    auto& current_layer2 = dynamic_cast<LayerData1DPeak&>(current_layer);
    const MSSpectrum& current_spectrum = current_layer2.getCurrentSpectrum();

    if (current_spectrum.empty())
    {
      OPENMS_LOG_WARN << "Spectrum is empty! Nothing to annotate!" << endl;
    }
    else if (!current_spectrum.isSorted())
    {
      QMessageBox::warning(tv_, "Error", "The spectrum is not sorted! Aborting!"); // @TODO: improve error message
      return;
    }

    // init all peak colors to black (=no annotation)
    current_layer2.peak_colors_1d.assign(current_spectrum.size(), Qt::black);

    for (const auto& ann : annotations)
    {
      // find matching peak in experimental spectrum
      Int peak_idx = current_spectrum.findNearest(ann.mz, 1e-2);
      if (peak_idx == -1) // no match
      {
        OPENMS_LOG_WARN << "Annotation present for missing peak. m/z: " << ann.mz
                  << endl;
        continue;
      }

      String label = ann.annotation;
      label.trim();

#ifdef DEBUG_IDENTIFICATION_VIEW
      cout << "Adding annotation item based on fragment annotations: " << label << endl;
#endif
      QStringList lines = label.toQString().split(QRegularExpression("[\r\n]"), Qt::SkipEmptyParts);
      if (lines.size() > 1)
      {
        label = String(lines[0]);
      }

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

      QColor color(Qt::black);
      QColor peak_color(Qt::black);

      // XL-MS specific coloring of the labels, green for linear fragments and red for cross-linked fragments
      if (label.hasSubstring("[alpha|")
        || label.hasSubstring("[beta|"))
      {
        if (label.hasSubstring("|ci$"))
        {
          color = Qt::darkGreen;
          peak_color = Qt::green;
        }
        else if (label.hasSubstring("|xi$"))
        {
          color = Qt::darkRed;
          peak_color = Qt::red;
        }
      }
      else
      { // different colors for left/right fragments (e.g. b/y ions)
        color = (label.at(0) < 'n') ? Qt::darkRed : Qt::darkGreen;
        peak_color = (label.at(0) < 'n') ? Qt::red : Qt::green;
      }

      Peak1D position(current_spectrum[peak_idx].getMZ(),
        current_spectrum[peak_idx].getIntensity());

      if (lines.size() > 1)
      {
        label.append("\n").append(String(lines[1]));
      }

      auto item = new Annotation1DPeakItem<Peak1D>(
        position,
        label.toQString(),
        color);

      // set peak color
      current_layer2.peak_colors_1d[peak_idx] = peak_color;

      item->setSelected(false);
      tv_->getActive1DWidget()->canvas()->getCurrentLayer().getCurrentAnnotations().push_front(item);
    }

    // Block update events for identification widget
    // TODO: Why? If it is to avoid a new selection while repainting, why just do it now and not much earlier??
    spec_id_view_->ignore_update = true;
    RAIICleanup cleanup([&]() {spec_id_view_->ignore_update = false; });

    // zoom visible area to real data range:
    auto spec_range = current_layer.getRange();
    spec_range.scaleBy(1.2);
    tv_->getActive1DWidget()->canvas()->setVisibleArea(RangeAllType().assign(spec_range));

    tv_->updateLayerBar();
  }

  void TVIdentificationViewController::removeTheoreticalSpectrumLayer_()
  {
    auto* spectrum_widget_1D = tv_->getActive1DWidget();
    if (spectrum_widget_1D)
    {
      Plot1DCanvas* canvas_1D = spectrum_widget_1D->canvas();

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

  // override
  void TVIdentificationViewController::activateBehavior() 
  {
    Plot1DWidget* w = tv_->getActive1DWidget();
    if (w == nullptr)
    {
      return;
    }
    auto* current_canvas = w->canvas();
    auto& current_layer = dynamic_cast<LayerData1DPeak&>(current_canvas->getCurrentLayer());

    const SpectrumType& current_spectrum = current_layer.getCurrentSpectrum();

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
        current_layer.setCurrentIndex(i);
        break;
      }
    }
  }

  // override
  void TVIdentificationViewController::deactivateBehavior()
  {
    Plot1DWidget* widget_1D = tv_->getActive1DWidget();

    // return if no active 1D widget is present
    if (widget_1D == nullptr)
    {
      return;
    }
    // clear textbox
    widget_1D->canvas()->setTextBox(QString());

    // remove precursor labels, theoretical spectra and trigger repaint
    auto cl = dynamic_cast<LayerData1DPeak*> (&tv_->getActive1DWidget()->canvas()->getCurrentLayer());
    if (!cl) return;
    removeTemporaryAnnotations_(cl->getCurrentIndex());
    removeTheoreticalSpectrumLayer_();
    cl->peptide_id_index = -1;
    cl->peptide_hit_index = -1;
    tv_->getActive1DWidget()->canvas()->repaint();
  }

  void TVIdentificationViewController::setVisibleArea1D(double l, double h)
  {
    Plot1DWidget* widget_1D = tv_->getActive1DWidget();

    // return if no active 1D widget is present
    if (widget_1D == nullptr)
    {
      return;
    }
    DRange<2> range = tv_->getActive1DWidget()->canvas()->getVisibleArea().getAreaXY();
    range.setMinX(l);
    range.setMaxX(h);
    tv_->getActive1DWidget()->canvas()->setVisibleArea(range);
  }
}
