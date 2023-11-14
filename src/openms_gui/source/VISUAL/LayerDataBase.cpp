// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerDataBase.h>

#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>// for AMS annotation
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/DATASTRUCTURES/OSWData.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/LayerDataConsensus.h>
#include <OpenMS/VISUAL/LayerDataFeature.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>

using namespace std;

namespace OpenMS
{
  LayerDataDefs::ProjectionData::ProjectionData() = default;
  LayerDataDefs::ProjectionData::ProjectionData(ProjectionData&&) = default;
  LayerDataDefs::ProjectionData::~ProjectionData() = default;

  const std::string LayerDataDefs::NamesOfLabelType[] = {"None", "Index", "Label meta data", "Peptide identification", "All peptide identifications"};

  std::ostream& operator<<(std::ostream& os, const LayerDataBase& rhs)
  {
    os << "--LayerDataBase BEGIN--\n";
    os << "name: " << rhs.getName() << '\n';
    os << "visible: " << rhs.visible << '\n';
    os << "--LayerDataBase END--\n";
    return os;
  }


  /// get name plus name_extra plus optionally augmented with attributes, e.g. '*' if modified
  String LayerDataBase::getDecoratedName() const
  {
    String n = name_ + name_suffix_;
    if (modified)
    {
      n += '*';
    }
    return n;
  }

  bool LayerDataBase::annotate(const vector<PeptideIdentification>& identifications,
                           const vector<ProteinIdentification>& protein_identifications)
  {
    IDMapper mapper;
    if (auto* lp = dynamic_cast<LayerDataPeak*>(this))
    {
      Param p = mapper.getDefaults();
      p.setValue("rt_tolerance", 0.1, "RT tolerance (in seconds) for the matching");
      p.setValue("mz_tolerance", 1.0, "m/z tolerance (in ppm or Da) for the matching");
      p.setValue("mz_measure", "Da", "unit of 'mz_tolerance' (ppm or Da)");
      mapper.setParameters(p);
      mapper.annotate(*lp->getPeakDataMuteable(), identifications, protein_identifications, true);
    }
    if (auto* lp = dynamic_cast<LayerDataFeature*>(this))
    {
      mapper.annotate(*lp->getFeatureMap(), identifications, protein_identifications);
    }
    else if (auto* lp = dynamic_cast<LayerDataConsensus*>(this))
    {
      mapper.annotate(*lp->getConsensusMap(), identifications, protein_identifications);
    }
    else
    {
      return false;
    }

    return false;
  }


  float LayerDataBase::getMinIntensity() const
  {
    return getRange().getMinIntensity();
  }

  float LayerDataBase::getMaxIntensity() const
  {
    return getRange().getMaxIntensity();
  }

/*
<<<<<<< HEAD
  void LayerDataBase::synchronizePeakAnnotations()
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->empty() || type != LayerDataBase::DT_PEAK)
    {
      return;
    }

    // get mutable access to the spectrum
    MSSpectrum& spectrum = getPeakDataMuteable()->getSpectrum(current_spectrum_idx_);

    int ms_level = spectrum.getMSLevel();

    if (ms_level != 2) return;

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
      {// no hits? add empty hit
        PeptideHit hit;
        updatePeptideHitAnnotations_(hit);
        hits.push_back(hit);
      }
    }
    else// PeptideIdentifications are empty, create new PepIDs and PeptideHits to store the PeakAnnotations
    {
      // copy user annotations to fragment annotation vector
      const Annotations1DContainer& las = getAnnotations(current_spectrum_idx_);

      // no annotations so we don't need to synchronize
      bool has_peak_annotation(false);
      for (auto& a : las)
      {
        // only store peak annotations
        Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(a);
        if (pa != nullptr)
        {
          has_peak_annotation = true;
          break;
        }
      }
      if (has_peak_annotation == false)
      {
        return;
      }
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

  void LayerDataBase::updatePeptideHitAnnotations_(PeptideHit& hit)
  {
    // copy user annotations to fragment annotation vector
    const Annotations1DContainer& las = getCurrentAnnotations();


    const vector<PeptideHit::PeakAnnotation>& old_fas = hit.getPeakAnnotations();

    // initialize with an empty vector
    vector<PeptideHit::PeakAnnotation> fas;

    // do not change PeptideHit annotations, if there are no annotations on the spectrum
    bool annotations_changed(false);

    // for each annotation item on the canvas
    for (auto& a : las)
    {
      // only store peak annotations (skip general labels and distance annotations)
      Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(a);
      if (pa == nullptr)
      {
        continue;
      }
      
      auto fa = pa->toPeakAnnotation();

      // check if peptide hit alread contains the annotation 
      // Note: we only match by m/z and label. Intensity is ignored to prevent updates from switching 
      // between relative and absolute mode change the annotation in unpredictable ways.
      auto it = std::find_if(old_fas.begin(), old_fas.end(),
       [&fa](const PeptideHit::PeakAnnotation& old_fa) { return std::abs(old_fa.mz - fa.mz) < 0.0001 && old_fa.annotation == fa.annotation ;});

      if (it == old_fas.end())
      {
        fas.push_back(fa);
        annotations_changed = true;
      }      
    }

    if (annotations_changed)
    {
      hit.setPeakAnnotations(fas);
    }
  }

  void LayerDataBase::removePeakAnnotationsFromPeptideHit(const std::vector<Annotation1DItem*>& selected_annotations)
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->empty() || type != LayerDataBase::DT_PEAK)
    {
      return;
    }

    // no ID selected
    if (peptide_id_index == -1 || peptide_hit_index == -1)
    {
      return;
    }

    // get mutable access to the spectrum
    MSSpectrum& spectrum = getPeakDataMuteable()->getSpectrum(current_spectrum_idx_);
    int ms_level = spectrum.getMSLevel();

    // wrong MS level
    if (ms_level < 2)
    {
      return;
    }

    // extract PeptideIdentification and PeptideHit if possible.
    // that this function returns prematurely is unlikely,
    // since we are deleting existing annotations,
    // that have to be somewhere, but better make sure
    vector<PeptideIdentification>& pep_ids = spectrum.getPeptideIdentifications();
    if (pep_ids.empty())
    {
      return;
    }
    vector<PeptideHit>& hits = pep_ids[peptide_id_index].getHits();
    if (hits.empty())
    {
      return;
    }
    PeptideHit& hit = hits[peptide_hit_index];
    vector<PeptideHit::PeakAnnotation> fas = hit.getPeakAnnotations();
    if (fas.empty())
    {
      return;
    }

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
        if (pa == nullptr)
        {
          continue;
        }
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
    if (annotations_changed)
    {
      hit.setPeakAnnotations(fas);
    }
  }

=======
>>>>>>> origin/develop
*/

  LayerAnnotatorBase::LayerAnnotatorBase(const FileTypeList& supported_types, const String& file_dialog_text, QWidget* gui_lock) :
      supported_types_(supported_types),
      file_dialog_text_(file_dialog_text),
      gui_lock_(gui_lock)
  {
  }

  bool LayerAnnotatorBase::annotateWithFileDialog(LayerDataBase& layer, LogWindow& log, const String& current_path) const
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
                                                 supported_types_.toFileDialogFilter(FilterLayout::BOTH, true).toQString());

    bool success = annotateWithFilename(layer, log, fname);

    return success;
  }

  bool LayerAnnotatorBase::annotateWithFilename(LayerDataBase& layer, LogWindow& log, const String& fname) const
  {
    if (fname.empty())
    {
      return false;
    }
    FileTypes::Type type = FileHandler::getType(fname);

    if (!supported_types_.contains(type))
    {
      log.appendNewHeader(LogWindow::LogState::NOTICE, "Error", String("Filename '" + fname + "' has unsupported file type. No annotation performed.").toQString());
      return false;
    }

    GUIHelpers::GUILock glock(gui_lock_);
    bool success = annotateWorker_(layer, fname, log);

    if (success)
    {
      log.appendNewHeader(LogWindow::LogState::NOTICE, "Done", "Annotation finished. Open the corresponding view to see results!");
    }
    return success;
  }

  std::unique_ptr<LayerAnnotatorBase> LayerAnnotatorBase::getAnnotatorWhichSupports(const FileTypes::Type& type)
  {
    std::unique_ptr<LayerAnnotatorBase> ptr(nullptr);
    auto assign = [&type, &ptr](std::unique_ptr<LayerAnnotatorBase> other) {
      if (other->supported_types_.contains(type))
      {
        if (ptr.get() != nullptr)
        {
          throw Exception::IllegalSelfOperation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
        ptr = std::move(other);
      }
    };
    // hint: add new derived classes here, so they are checked as well
    assign(std::unique_ptr<LayerAnnotatorBase>(new LayerAnnotatorAMS(nullptr)));
    assign(std::unique_ptr<LayerAnnotatorBase>(new LayerAnnotatorPeptideID(nullptr)));
    assign(std::unique_ptr<LayerAnnotatorBase>(new LayerAnnotatorOSW(nullptr)));

    return ptr;// Note: no std::move here needed because of copy elision
  }

  std::unique_ptr<LayerAnnotatorBase> LayerAnnotatorBase::getAnnotatorWhichSupports(const String& filename)
  {
    return getAnnotatorWhichSupports(FileHandler::getType(filename));
  }

  bool LayerAnnotatorPeptideID::annotateWorker_(LayerDataBase& layer, const String& filename, LogWindow& /*log*/) const
  {
    FileTypes::Type type = FileHandler::getType(filename);
    vector<PeptideIdentification> identifications;
    vector<ProteinIdentification> protein_identifications;
    FileHandler().loadIdentifications(filename, protein_identifications, identifications, {type});

    layer.annotate(identifications, protein_identifications);
    return true;
  }

  bool LayerAnnotatorAMS::annotateWorker_(LayerDataBase& layer, const String& filename, LogWindow& log) const
  {
    FeatureMap fm;
    FileHandler().loadFeatures(filename, fm, {FileTypes::FEATUREXML});

    // last protein ID must be from AccurateMassSearch (it gets appended there)
    String engine = "no protein identification section found";
    if (fm.getProteinIdentifications().size() > 0)
    {
      engine = fm.getProteinIdentifications().back().getSearchEngine();
      if (engine == AccurateMassSearchEngine::search_engine_identifier)
      {
        auto* lp = dynamic_cast<LayerDataPeak*>(&layer);
        if (!lp)
        {
          QMessageBox::warning(nullptr, "Error", "Layer type is not DT_PEAK!");
          return false;
        }
        IDMapper im;
        Param p = im.getParameters();
        p.setValue("rt_tolerance", 30.0);
        im.setParameters(p);
        log.appendNewHeader(LogWindow::LogState::NOTICE, "Note", "Mapping matches with 30 sec tolerance and no m/z limit to spectra...");
        im.annotate((*lp->getPeakDataMuteable()), fm, true, true);

        return true;
      }
    }

    QMessageBox::warning(nullptr, "Error", (String("FeatureXML is currently only supported for files generated by the AccurateMassSearch tool (got '") + engine + "', expected 'AccurateMassSearch'.").toQString());
    return false;
  }

  bool LayerAnnotatorOSW::annotateWorker_(LayerDataBase& layer,
                                          const String& filename,
                                          LogWindow& log) const
  {
    log.appendNewHeader(LogWindow::LogState::NOTICE, "Note", "Reading OSW data ...");
    auto* lp = dynamic_cast<LayerDataChrom*>(&layer);
    if (!lp)
    {
      QMessageBox::warning(nullptr, "Error", "Layer type is not DT_CHROM!");
      return false;
    }
    try
    {
      OSWFile oswf(filename);// this can throw if file does not exist
      OSWData data;
      oswf.readMinimal(data);
      // allow data to map from transition.id (=native.id) to a chromatogram index in MSExperiment
      data.buildNativeIDResolver(*lp->getChromatogramData().get());
      lp->setChromatogramAnnotation(std::move(data));
      return true;
    }
    catch (Exception::BaseException& e)
    {
      log.appendText(e.what());
      return false;
    }
  }


} // namespace OpenMS
