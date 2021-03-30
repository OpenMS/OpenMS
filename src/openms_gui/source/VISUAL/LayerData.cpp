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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerData.h>

#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h> // for AMS annotation
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/DATASTRUCTURES/OSWData.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

//#include <iostream>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>

using namespace std;

namespace OpenMS
{
  const std::string LayerData::NamesOfLabelType[] = {"None", "Index", "Label meta data", "Peptide identification", "All peptide identifications"};

  std::ostream & operator<<(std::ostream & os, const LayerData & rhs)
  {
    os << "--LayerData BEGIN--" << std::endl;
    os << "name: " << rhs.getName() << std::endl;
    os << "visible: " << rhs.visible << std::endl;
    os << "number of peaks: " << rhs.getPeakData()->getSize() << std::endl;
    os << "--LayerData END--" << std::endl;
    return os;
  }


  /// Default constructor

  LayerData::LayerData() :
    flags(),
    visible(true),
    flipped(false),
    type(DT_UNKNOWN),
    name_(),
    filename(),
    peptides(),
    param(),
    gradient(),
    filters(),
    annotations_1d(),
    peak_colors_1d(),
    modifiable(false),
    modified(false),
    label(L_NONE),
    peptide_id_index(-1),
    peptide_hit_index(-1),
    features_(new FeatureMapType()),
    consensus_map_(new ConsensusMapType()),
    peak_map_(new ExperimentType()),
    on_disc_peaks(new OnDiscMSExperiment()),
    chromatogram_map_(new ExperimentType()),
    current_spectrum_(0),
    cached_spectrum_()
  {
    annotations_1d.resize(1);
  }

  const LayerData::ConstExperimentSharedPtrType LayerData::getPeakData() const
  {
    return boost::static_pointer_cast<const ExperimentType>(peak_map_);
  }

  void LayerData::updateRanges()
  {
    peak_map_->updateRanges();
    features_->updateRanges();
    consensus_map_->updateRanges();
    // on_disc_peaks->updateRanges(); // note: this is not going to work since its on disk! We currently don't have a good way to access these ranges
    chromatogram_map_->updateRanges();
    cached_spectrum_.updateRanges();
  }

  /// Returns the minimum intensity of the internal data, depending on type

  float LayerData::getMinIntensity() const
  {
    if (type == LayerData::DT_PEAK || type == LayerData::DT_CHROMATOGRAM)
    {
      return getPeakData()->getMinInt();
    }
    else if (type == LayerData::DT_FEATURE)
    {
      return getFeatureMap()->getMinInt();
    }
    else
    {
      return getConsensusMap()->getMinInt();
    }
  }

  /// Returns the maximum intensity of the internal data, depending on type

  float LayerData::getMaxIntensity() const
  {
    if (type == LayerData::DT_PEAK || type == LayerData::DT_CHROMATOGRAM)
    {
      return getPeakData()->getMaxInt();
    }
    else if (type == LayerData::DT_FEATURE)
    {
      return getFeatureMap()->getMaxInt();
    }
    else
    {
      return getConsensusMap()->getMaxInt();
    }
  }


  /// get name augmented with attributes, e.g. [flipped], or '*' if modified

   String LayerData::getDecoratedName() const
   {
    String n = name_;
    if (flipped)
    {
      n += " [flipped]";
    }
    if (modified)
    {
      n += '*';
    }
    return n;
  }

  void LayerData::updateCache_()
  {
    if (peak_map_->getNrSpectra() > current_spectrum_ && (*peak_map_)[current_spectrum_].size() > 0)
    {
      cached_spectrum_ = (*peak_map_)[current_spectrum_];
    }
    else if (on_disc_peaks->getNrSpectra() > current_spectrum_)
    {
      cached_spectrum_ = on_disc_peaks->getSpectrum(current_spectrum_);
    }
  }


  /// add annotation from an OSW sqlite file.


  /// get annotation (e.g. to build a hierachical ID View)
  /// Not const, because we might have incomplete data, which needs to be loaded from sql source

  LayerData::OSWDataSharedPtrType& LayerData::getChromatogramAnnotation()
  {
    return chrom_annotation_;
  }

  const LayerData::OSWDataSharedPtrType& LayerData::getChromatogramAnnotation() const
  {
    return chrom_annotation_;
  }

  void LayerData::setChromatogramAnnotation(OSWData&& data)
  {
    chrom_annotation_ = OSWDataSharedPtrType(new OSWData(std::move(data)));
  }

  bool LayerData::annotate(const vector<PeptideIdentification>& identifications,
    const vector<ProteinIdentification>& protein_identifications)
  {
    IDMapper mapper;
    if (this->type == LayerData::DT_PEAK)
    {
      Param p = mapper.getDefaults();
      p.setValue("rt_tolerance", 0.1, "RT tolerance (in seconds) for the matching");
      p.setValue("mz_tolerance", 1.0, "m/z tolerance (in ppm or Da) for the matching");
      p.setValue("mz_measure", "Da", "unit of 'mz_tolerance' (ppm or Da)");
      mapper.setParameters(p);
      mapper.annotate(*getPeakDataMuteable(), identifications, protein_identifications, true);
    }
    else if (type == LayerData::DT_FEATURE)
    {
      mapper.annotate(*getFeatureMap(), identifications, protein_identifications);
    }
    else if (type == LayerData::DT_CONSENSUS)
    {
      mapper.annotate(*getConsensusMap(), identifications, protein_identifications);
    }
    else
    {
      return false;
    }

    return true;
  }

  const LayerData::ExperimentType::SpectrumType& LayerData::getCurrentSpectrum() const
  {
    return cached_spectrum_;
  }

  /// Returns a const-copy of the required spectrum which is guaranteed to be populated with raw data

  const LayerData::ExperimentType::SpectrumType LayerData::getSpectrum(Size spectrum_idx) const
  {
    if (spectrum_idx == current_spectrum_) return cached_spectrum_;

    if ((*peak_map_)[spectrum_idx].size() > 0)
    {
      return (*peak_map_)[spectrum_idx];
    }
    else if (!on_disc_peaks->empty())
    {
      return on_disc_peaks->getSpectrum(spectrum_idx);
    }
    return (*peak_map_)[spectrum_idx];
  }

  void LayerData::synchronizePeakAnnotations()
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->empty() || type != LayerData::DT_PEAK) { return; }

    // get mutable access to the spectrum
    MSSpectrum & spectrum = getPeakDataMuteable()->getSpectrum(current_spectrum_);

    int ms_level = spectrum.getMSLevel();

    if (ms_level == 2)
    {
      // store user fragment annotations
      vector<PeptideIdentification>& pep_ids = spectrum.getPeptideIdentifications();

      // no ID selected
      if (peptide_id_index == -1 || peptide_hit_index == -1)
      {
        return;
      }

      if (!pep_ids.empty())
      {
        vector<PeptideHit>& hits = pep_ids[peptide_id_index].getHits();

        if (!hits.empty())
        {
          PeptideHit& hit = hits[peptide_hit_index];
          updatePeptideHitAnnotations_(hit);
        }
        else
        { // no hits? add empty hit
          PeptideHit hit;
          updatePeptideHitAnnotations_(hit);
          hits.push_back(hit);
        }
      }
      else // PeptideIdentifications are empty, create new PepIDs and PeptideHits to store the PeakAnnotations
      {
        // copy user annotations to fragment annotation vector
        const Annotations1DContainer & las = getAnnotations(current_spectrum_);

        // no annotations so we don't need to synchronize
        bool has_peak_annotation(false);
        for (auto& a : las)
        {
          // only store peak annotations
          Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(a);
          if (pa != nullptr) { has_peak_annotation = true; break; }
        }
        if (has_peak_annotation == false) return;

        PeptideIdentification pep_id;
        pep_id.setIdentifier("Unknown");

        // create a dummy ProteinIdentification for all ID-less PeakAnnotations
        vector<ProteinIdentification>& prot_ids = getPeakDataMuteable()->getProteinIdentifications();
        if (prot_ids.empty() || prot_ids.back().getIdentifier() != String("Unknown"))
        {
          ProteinIdentification prot_id;
          prot_id.setIdentifier("Unknown");
          prot_ids.push_back(prot_id);
        }

        PeptideHit hit;
        if (spectrum.getPrecursors().empty() == false)
        {
          pep_id.setMZ(spectrum.getPrecursors()[0].getMZ());
          hit.setCharge(spectrum.getPrecursors()[0].getCharge());
        }
        pep_id.setRT(spectrum.getRT());

        updatePeptideHitAnnotations_(hit);
        std::vector<PeptideHit> hits;
        hits.push_back(hit);
        pep_id.setHits(hits);
        pep_ids.push_back(pep_id);
      }
    }
  }

  void LayerData::updatePeptideHitAnnotations_(PeptideHit& hit)
  {
    // copy user annotations to fragment annotation vector
    const Annotations1DContainer & las = getCurrentAnnotations();

    // initialize with an empty vector
    vector<PeptideHit::PeakAnnotation> fas;

    // do not change PeptideHit annotations, if there are no annotations on the spectrum
    bool annotations_changed(false);

    // regular expression for a charge at the end of the annotation
    QRegExp reg_exp("([\\+|\\-]\\d+)$");

    // for each annotation item on the canvas
    for (auto& a : las)
    {
      // only store peak annotations (skip general lables and distance annotations)
      Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(a);
      if (pa == nullptr) { continue; }

      // add new fragment annotation
      QString peak_anno = pa->getText().trimmed();

      // check for newlines in the label and only continue with the first line for charge determination
      QStringList lines = peak_anno.split(QRegExp("[\r\n]"), QString::SkipEmptyParts);
      if (lines.size() > 1)
      {
        peak_anno = lines[0];
      }

      // read charge and text from annotation item string
      // we support two notations for the charge suffix: '2+' or '++'
      // cut and convert the trailing + or - to a proper charge
      int match_pos = reg_exp.indexIn(peak_anno);
      int tmp_charge(0);
      if (match_pos >= 0)
      {
        tmp_charge = reg_exp.cap(1).toInt();
        peak_anno = peak_anno.left(match_pos);
      }
      else
      {
        // count number of + and - in suffix (e.g., to support "++" as charge 2 anotation)
        int plus(0), minus(0);

        for (int p = (int)peak_anno.size() - 1; p >= 0; --p)
        {
          if (peak_anno[p] == '+')
          {
            ++plus;
            continue;
          }
          else if (peak_anno[p] == '-')
          {
            ++minus;
            continue;
          }
          else // not '+' or '-'?
          {
            if (plus > 0 && minus == 0) // found pluses?
            {
              tmp_charge = plus;
              peak_anno = peak_anno.left(peak_anno.size() - plus);
              break;
            }
            else if (minus > 0 && plus == 0)  // found minuses?
            {
              tmp_charge = -minus;
              peak_anno = peak_anno.left(peak_anno.size() - minus);
              break;
            }
            break;
          }
        }
      }

      PeptideHit::PeakAnnotation fa;
      fa.charge = tmp_charge;
      fa.mz = pa->getPeakPosition()[0];
      fa.intensity = pa->getPeakPosition()[1];
      if (lines.size() > 1)
      {
        peak_anno.append("\n").append(lines[1]);
      }
      fa.annotation = peak_anno;

      fas.push_back(fa);
      annotations_changed = true;
    }

    if (annotations_changed)
    {
      hit.setPeakAnnotations(fas);
    }
  }

  void LayerData::removePeakAnnotationsFromPeptideHit(const std::vector<Annotation1DItem*>& selected_annotations)
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->empty() || type != LayerData::DT_PEAK) { return; }

    // no ID selected
    if (peptide_id_index == -1 || peptide_hit_index == -1) { return; }

    // get mutable access to the spectrum
    MSSpectrum & spectrum = getPeakDataMuteable()->getSpectrum(current_spectrum_);
    int ms_level = spectrum.getMSLevel();

    // wrong MS level
    if (ms_level < 2) { return; }

    // extract PeptideIdentification and PeptideHit if possible.
    // that this function returns prematurely is unlikely,
    // since we are deleting existing annotations,
    // that have to be somewhere, but better make sure
    vector<PeptideIdentification>& pep_ids = spectrum.getPeptideIdentifications();
    if (pep_ids.empty()) { return; }
    vector<PeptideHit>& hits = pep_ids[peptide_id_index].getHits();
    if (hits.empty()) { return; }
    PeptideHit& hit = hits[peptide_hit_index];
    vector<PeptideHit::PeakAnnotation> fas = hit.getPeakAnnotations();
    if (fas.empty()) { return; }

    // all requirements fulfilled, PH in hit and annotations in selected_annotations
    vector<PeptideHit::PeakAnnotation> to_remove;
    bool annotations_changed(false);

    // collect annotations, that have to be removed
    for (auto const& tmp_a : fas)
    {
      for (auto const& it : selected_annotations)
      {
        Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(it);
        // only search for peak annotations
        if (pa == nullptr) { continue; }
        if (fabs(tmp_a.mz - pa->getPeakPosition()[0]) < 1e-6)
        {
          if (String(pa->getText()).hasPrefix(tmp_a.annotation))
          {
            to_remove.push_back(tmp_a);
            annotations_changed = true;
          }
        }
      }
    }
    // remove the collected annotations from the PeptideHit annotations
    for (auto const& tmp_a : to_remove)
    {
      fas.erase(std::remove(fas.begin(), fas.end(), tmp_a), fas.end());
    }
    if (annotations_changed) { hit.setPeakAnnotations(fas); }
  }

  LayerAnnotatorBase::LayerAnnotatorBase(const FileTypes::FileTypeList& supported_types, const String& file_dialog_text, QWidget* gui_lock)
    : supported_types_(supported_types),
      file_dialog_text_(file_dialog_text),
      gui_lock_(gui_lock)
  {
  }

  bool LayerAnnotatorBase::annotateWithFileDialog(LayerData& layer, LogWindow& log, const String& current_path) const
  {
    // warn if hidden layer => wrong layer selected...
    if (!layer.visible)
    {
      log.appendNewHeader(LogWindow::LogState::NOTICE, "The current layer is not visible", "Have you selected the right layer for this action? Aborting.");
      return false;
    }

    // load id data
    QString fname = QFileDialog::getOpenFileName(nullptr,
                                                 file_dialog_text_.toQString(),
                                                 current_path.toQString(),
                                                 supported_types_.toFileDialogFilter(FileTypes::Filter::BOTH, true).toQString());
    
    bool success = annotateWithFilename(layer, log, fname);

    return success;
  }

  bool LayerAnnotatorBase::annotateWithFilename(LayerData& layer, LogWindow& log, const String& fname) const
  {
    if (fname.empty()) return false;

    FileTypes::Type type = FileHandler::getType(fname);

    if (!supported_types_.contains(type))
    {
      log.appendNewHeader(LogWindow::LogState::NOTICE, "Error", String("Filename '" + fname + "' has unsupported file type. No annotation performed.").toQString());
      return false;
    }

    GUIHelpers::GUILock glock(gui_lock_);
    bool success = annotateWorker_(layer, fname, log);

    if (success) log.appendNewHeader(LogWindow::LogState::NOTICE, "Done", "Annotation finished. Open identification view to see results!");
    return success;
  }

  std::unique_ptr<LayerAnnotatorBase> LayerAnnotatorBase::getAnnotatorWhichSupports(const FileTypes::Type& type)
  {
    std::unique_ptr<LayerAnnotatorBase> ptr(nullptr);
    auto assign = [&type, &ptr](std::unique_ptr<LayerAnnotatorBase> other)
    {
      if (other->supported_types_.contains(type))
      {
        if (ptr.get() != nullptr) throw Exception::IllegalSelfOperation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        ptr = std::move(other);
      }
    };
    // hint: add new derived classes here, so they are checked as well
    assign(std::unique_ptr<LayerAnnotatorBase>(new LayerAnnotatorAMS(nullptr)));
    assign(std::unique_ptr<LayerAnnotatorBase>(new LayerAnnotatorPeptideID(nullptr)));
    assign(std::unique_ptr<LayerAnnotatorBase>(new LayerAnnotatorOSW(nullptr)));

    return ptr; // Note: no std::move here needed because of copy elision
  }

  std::unique_ptr<LayerAnnotatorBase> LayerAnnotatorBase::getAnnotatorWhichSupports(const String& filename)
  {
    return getAnnotatorWhichSupports(FileHandler::getType(filename));
  }

  bool LayerAnnotatorPeptideID::annotateWorker_(LayerData& layer, const String& filename, LogWindow& /*log*/) const
  {
    FileTypes::Type type = FileHandler::getType(filename);
    vector<PeptideIdentification> identifications;
    vector<ProteinIdentification> protein_identifications;
    if (type == FileTypes::MZIDENTML)
    {
      MzIdentMLFile().load(filename, protein_identifications, identifications);
    }
    else
    {
      String document_id;
      IdXMLFile().load(filename, protein_identifications, identifications, document_id);
    }

    layer.annotate(identifications, protein_identifications);
    return true;
  }

  bool LayerAnnotatorAMS::annotateWorker_(LayerData& layer, const String& filename, LogWindow& log) const
  {
    FeatureMap fm;
    FeatureXMLFile().load(filename, fm);

    // last protein ID must be from AccurateMassSearch (it gets appended there)
    String engine = "no protein identification section found";
    if (fm.getProteinIdentifications().size() > 0)
    {
      engine = fm.getProteinIdentifications().back().getSearchEngine();
      if (engine == AccurateMassSearchEngine::search_engine_identifier)
      {
        if (layer.type != LayerData::DT_PEAK)
        {
          QMessageBox::warning(nullptr, "Error", "Layer type is not DT_PEAK!");
          return false;
        }
        IDMapper im;
        Param p = im.getParameters();
        p.setValue("rt_tolerance", 30.0);
        im.setParameters(p);
        log.appendNewHeader(LogWindow::LogState::NOTICE, "Note", "Mapping matches with 30 sec tolerance and no m/z limit to spectra...");
        im.annotate((*layer.getPeakDataMuteable()), fm, true, true);

        return true;
      }
    }

    QMessageBox::warning(nullptr, "Error", (String("FeatureXML is currently only supported for files generated by the AccurateMassSearch tool (got '") + engine + "', expected 'AccurateMassSearch'.").toQString());
    return false;
  }

  bool LayerAnnotatorOSW::annotateWorker_(LayerData& layer,
                                          const String& filename,
                                          LogWindow& log) const
  {
    log.appendNewHeader(LogWindow::LogState::NOTICE, "Note", "Reading OSW data ...");
    try
    {
      OSWFile oswf(filename); // this can throw if file does not exist
      OSWData data;
      oswf.readMinimal(data);
      // allow data to map from transition.id (=native.id) to a chromatogram index in MSExperiment
      data.buildNativeIDResolver(*layer.getFullChromData().get());
      layer.setChromatogramAnnotation(std::move(data));
      return true;
    }
    catch (Exception::BaseException& e)
    {
      log.appendText(e.what());
      return false;
    }
  }

} // namespace OpenMS
