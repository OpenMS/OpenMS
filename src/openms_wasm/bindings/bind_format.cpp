// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Reading and writing files (FORMAT/). Paths refer to the module's virtual file system
// (OpenMS.FS, WasmFS): write a file there first (FS.writeFile) or, in Node.js, mount a
// directory of the host (mountHostDirectory). MzMLFile also reads and writes strings directly.

#include "binding_utils.h"

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace emscripten;

namespace OpenMS::WASM
{
  using ProteinIdentificationList = std::vector<ProteinIdentification>;

  namespace
  {
    //------------------------------------------------------------------------------
    // Load options of peak files, as plain object:
    // {msLevels: [1, 2], rtRange: [min, max], mzRange: [min, max], intensityRange: [min, max],
    //  metadataOnly, fillData, sortSpectraByMZ, skipXMLChecks}

    DRange<1> rangeFromJS(const val& v)
    {
      const auto bounds = numbersFromJS<double>(v);
      if (bounds.size() != 2)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "A range is an array [min, max].");
      }
      return DRange<1>(DPosition<1>(bounds[0]), DPosition<1>(bounds[1]));
    }

    void applyPeakFileOptions(PeakFileOptions& o, const val& v)
    {
      if (!isNullish(v["msLevels"])) o.setMSLevels(numbersFromJS<Int>(v["msLevels"]));
      if (!isNullish(v["rtRange"])) o.setRTRange(rangeFromJS(v["rtRange"]));
      if (!isNullish(v["mzRange"])) o.setMZRange(rangeFromJS(v["mzRange"]));
      if (!isNullish(v["intensityRange"])) o.setIntensityRange(rangeFromJS(v["intensityRange"]));
      if (!isNullish(v["metadataOnly"])) o.setMetadataOnly(v["metadataOnly"].as<bool>());
      if (!isNullish(v["fillData"])) o.setFillData(v["fillData"].as<bool>());
      if (!isNullish(v["sortSpectraByMZ"])) o.setSortSpectraByMZ(v["sortSpectraByMZ"].as<bool>());
      if (!isNullish(v["skipXMLChecks"])) o.setSkipXMLChecks(v["skipXMLChecks"].as<bool>());
    }

    val rangeToJS(const DRange<1>& r)
    {
      val a = val::array();
      a.call<void>("push", r.minPosition()[0]);
      a.call<void>("push", r.maxPosition()[0]);
      return a;
    }

    JS::PeakFileOptions peakFileOptionsToJS(const PeakFileOptions& o)
    {
      val v = val::object();
      v.set("msLevels", toArray(o.getMSLevels()));
      if (o.hasRTRange()) v.set("rtRange", rangeToJS(o.getRTRange()));
      if (o.hasMZRange()) v.set("mzRange", rangeToJS(o.getMZRange()));
      if (o.hasIntensityRange()) v.set("intensityRange", rangeToJS(o.getIntensityRange()));
      v.set("metadataOnly", o.getMetadataOnly());
      v.set("fillData", o.getFillData());
      v.set("sortSpectraByMZ", o.getSortSpectraByMZ());
      v.set("skipXMLChecks", o.getSkipXMLChecks());
      return v;
    }

    //------------------------------------------------------------------------------
    // MzMLFile

    void mzmlSetOptions(MzMLFile& f, const JS::PeakFileOptions& options)
    {
      PeakFileOptions o = f.getOptions();
      applyPeakFileOptions(o, options);
      f.setOptions(o);
    }
    JS::PeakFileOptions mzmlGetOptions(const MzMLFile& f) { return peakFileOptionsToJS(f.getOptions()); }
    void mzmlLoad(MzMLFile& f, const std::string& path, PeakMap& exp) { f.load(path, exp); }
    void mzmlStore(const MzMLFile& f, const std::string& path, const PeakMap& exp) { f.store(path, exp); }
    void mzmlLoadBuffer(MzMLFile& f, const std::string& buffer, PeakMap& exp) { f.loadBuffer(buffer, exp); }
    std::string mzmlStoreBuffer(const MzMLFile& f, const PeakMap& exp)
    {
      std::string out;
      f.storeBuffer(out, exp);
      return out;
    }

    //------------------------------------------------------------------------------
    // Other formats

    void mzxmlLoad(MzXMLFile& f, const std::string& path, PeakMap& exp) { f.load(path, exp); }
    void mzxmlStore(const MzXMLFile& f, const std::string& path, const PeakMap& exp) { f.store(path, exp); }

    void mgfLoad(MascotGenericFile& f, const std::string& path, PeakMap& exp) { f.load(path, exp); }
    void mgfStore(MascotGenericFile& f, const std::string& path, const PeakMap& exp) { f.store(path, exp); }
    void mgfStoreCompact(MascotGenericFile& f, const std::string& path, const PeakMap& exp, bool compact) { f.store(path, exp, compact); }

    void featureXMLLoad(FeatureXMLFile& f, const std::string& path, FeatureMap& map) { f.load(path, map); }
    void featureXMLStore(FeatureXMLFile& f, const std::string& path, const FeatureMap& map) { f.store(path, map); }

    void consensusXMLLoad(ConsensusXMLFile& f, const std::string& path, ConsensusMap& map) { f.load(path, map); }
    void consensusXMLStore(ConsensusXMLFile& f, const std::string& path, const ConsensusMap& map) { f.store(path, map); }

    void idXMLLoad(IdXMLFile& f, const std::string& path, ProteinIdentificationList& proteins, PeptideIdentificationList& peptides)
    {
      f.load(path, proteins, peptides);
    }
    void idXMLStore(IdXMLFile& f, const std::string& path, const ProteinIdentificationList& proteins, const PeptideIdentificationList& peptides)
    {
      f.store(path, proteins, peptides);
    }

    void mzIdentMLLoad(MzIdentMLFile& f, const std::string& path, ProteinIdentificationList& proteins, PeptideIdentificationList& peptides)
    {
      f.load(path, proteins, peptides);
    }
    void mzIdentMLStore(const MzIdentMLFile& f, const std::string& path, const ProteinIdentificationList& proteins, const PeptideIdentificationList& peptides)
    {
      f.store(path, proteins, peptides);
    }

    /// FASTA entries as plain objects {identifier, description, sequence}.
    JS::FASTAEntries fastaLoad(const FASTAFile& f, const std::string& path)
    {
      std::vector<FASTAFile::FASTAEntry> entries;
      f.load(path, entries);
      return toArray(entries, [](const FASTAFile::FASTAEntry& e) {
        val o = val::object();
        o.set("identifier", e.identifier);
        o.set("description", e.description);
        o.set("sequence", e.sequence);
        return o;
      });
    }

    void fastaStore(const FASTAFile& f, const std::string& path, const JS::FASTAEntries& entries)
    {
      std::vector<FASTAFile::FASTAEntry> result;
      const unsigned n = static_cast<unsigned>(lengthOf(entries));
      for (unsigned i = 0; i < n; ++i)
      {
        const val e = entries[i];
        result.emplace_back(e["identifier"].as<std::string>(),
                            isNullish(e["description"]) ? std::string() : e["description"].as<std::string>(),
                            e["sequence"].as<std::string>());
      }
      f.store(path, result);
    }

    //------------------------------------------------------------------------------
    // FileHandler: any supported format, chosen by file name (or content)

    std::string handlerGetType(const std::string& path) { return FileTypes::typeToName(FileHandler::getType(path)); }
    std::string handlerGetTypeByFileName(const std::string& path) { return FileTypes::typeToName(FileHandler::getTypeByFileName(path)); }
    std::string handlerGetTypeByContent(const std::string& path) { return FileTypes::typeToName(FileHandler::getTypeByContent(path)); }

    void handlerSetOptions(FileHandler& h, const JS::PeakFileOptions& options)
    {
      PeakFileOptions o = h.getOptions();
      applyPeakFileOptions(o, options);
      h.setOptions(o);
    }
    JS::PeakFileOptions handlerGetOptions(const FileHandler& h) { return peakFileOptionsToJS(h.getOptions()); }
    void handlerLoadExperiment(FileHandler& h, const std::string& path, PeakMap& exp) { h.loadExperiment(path, exp); }
    void handlerStoreExperiment(FileHandler& h, const std::string& path, const PeakMap& exp) { h.storeExperiment(path, exp); }
    void handlerLoadFeatures(FileHandler& h, const std::string& path, FeatureMap& map) { h.loadFeatures(path, map); }
    void handlerStoreFeatures(FileHandler& h, const std::string& path, const FeatureMap& map) { h.storeFeatures(path, map); }
    void handlerLoadConsensusFeatures(FileHandler& h, const std::string& path, ConsensusMap& map) { h.loadConsensusFeatures(path, map); }
    void handlerStoreConsensusFeatures(FileHandler& h, const std::string& path, const ConsensusMap& map) { h.storeConsensusFeatures(path, map); }
    void handlerLoadIdentifications(FileHandler& h, const std::string& path, ProteinIdentificationList& proteins, PeptideIdentificationList& peptides)
    {
      h.loadIdentifications(path, proteins, peptides);
    }
    void handlerStoreIdentifications(FileHandler& h, const std::string& path, const ProteinIdentificationList& proteins, const PeptideIdentificationList& peptides)
    {
      h.storeIdentifications(path, proteins, peptides);
    }
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS;
using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_format)
{
  class_<MzMLFile>("MzMLFile")
    .constructor<>()
    .function("load", OMS_FN(&mzmlLoad))
    .function("store", OMS_FN(&mzmlStore))
    .function("loadBuffer", OMS_FN(&mzmlLoadBuffer))
    .function("storeBuffer", OMS_FN(&mzmlStoreBuffer))
    .function("getOptions", OMS_FN(&mzmlGetOptions))
    .function("setOptions", OMS_FN(&mzmlSetOptions))
    ;

  class_<MzXMLFile>("MzXMLFile")
    .constructor<>()
    .function("load", OMS_FN(&mzxmlLoad))
    .function("store", OMS_FN(&mzxmlStore))
    ;

  class_<MascotGenericFile, base<DefaultParamHandler>>("MascotGenericFile")
    .constructor<>()
    .function("load", OMS_FN(&mgfLoad))
    .function("store", OMS_FN(&mgfStore))
    .function("store", OMS_FN(&mgfStoreCompact))
    ;

  class_<FeatureXMLFile>("FeatureXMLFile")
    .constructor<>()
    .function("load", OMS_FN(&featureXMLLoad))
    .function("store", OMS_FN(&featureXMLStore))
    ;

  class_<ConsensusXMLFile>("ConsensusXMLFile")
    .constructor<>()
    .function("load", OMS_FN(&consensusXMLLoad))
    .function("store", OMS_FN(&consensusXMLStore))
    ;

  class_<IdXMLFile>("IdXMLFile")
    .constructor<>()
    .function("load", OMS_FN(&idXMLLoad))
    .function("store", OMS_FN(&idXMLStore))
    ;

  class_<MzIdentMLFile>("MzIdentMLFile")
    .constructor<>()
    .function("load", OMS_FN(&mzIdentMLLoad))
    .function("store", OMS_FN(&mzIdentMLStore))
    ;

  class_<FASTAFile>("FASTAFile")
    .constructor<>()
    .function("load", OMS_FN(&fastaLoad))
    .function("store", OMS_FN(&fastaStore))
    ;

  class_<FileHandler>("FileHandler")
    .constructor<>()
    .class_function("getType", OMS_FN(&handlerGetType))
    .class_function("getTypeByFileName", OMS_FN(&handlerGetTypeByFileName))
    .class_function("getTypeByContent", OMS_FN(&handlerGetTypeByContent))
    .function("getOptions", OMS_FN(&handlerGetOptions))
    .function("setOptions", OMS_FN(&handlerSetOptions))
    .function("loadExperiment", OMS_FN(&handlerLoadExperiment))
    .function("storeExperiment", OMS_FN(&handlerStoreExperiment))
    .function("loadFeatures", OMS_FN(&handlerLoadFeatures))
    .function("storeFeatures", OMS_FN(&handlerStoreFeatures))
    .function("loadConsensusFeatures", OMS_FN(&handlerLoadConsensusFeatures))
    .function("storeConsensusFeatures", OMS_FN(&handlerStoreConsensusFeatures))
    .function("loadIdentifications", OMS_FN(&handlerLoadIdentifications))
    .function("storeIdentifications", OMS_FN(&handlerStoreIdentifications))
    ;
}
