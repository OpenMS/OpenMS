// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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


  /**
   * @brief Returns the layer's name with any suffix and an asterisk if modified.
   *
   * The returned string consists of the layer's base name, an optional suffix, and appends '*' if the layer has been modified.
   *
   * @return The decorated name of the layer.
   */
  String LayerDataBase::getDecoratedName() const
  {
    String n = name_ + name_suffix_;
    if (modified)
    {
      n += '*';
    }
    return n;
  }

  /**
   * @brief Annotates the layer with peptide and protein identifications.
   *
   * Applies peptide and protein identification data to the layer's underlying data structure using ID mapping.
   * The annotation method and parameters depend on the specific layer type (peak, feature, or consensus).
   * Returns false if the layer type is unsupported or after annotation.
   *
   * @param identifications Peptide identifications to be mapped to the layer.
   * @param protein_identifications Protein identifications to be mapped to the layer.
   * @return false Always returns false, regardless of annotation outcome.
   */
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

  /**
   * @brief Returns the minimum intensity value from the layer's data range.
   *
   * @return The minimum intensity as a float.
   */
  float LayerDataBase::getMinIntensity() const
  {
    return getRange().getMinIntensity();
  }

  float LayerDataBase::getMaxIntensity() const
  {
    return getRange().getMaxIntensity();
  }

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

  /**
   * @brief Annotates a chromatogram layer with OSW data from a file.
   *
   * Loads OSW data from the specified file, builds a native ID resolver using the layer's chromatogram data, and sets the annotation on the layer. Returns false if the layer is not a chromatogram layer, if the file cannot be read, or if an exception occurs during processing.
   *
   * @param filename Path to the OSW file to load.
   * @return true if annotation was successful, false otherwise.
   */
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
      data.buildNativeIDResolver(lp->getChromatogramData().get()->getMSExperiment());
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
