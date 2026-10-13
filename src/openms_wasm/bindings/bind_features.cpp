// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Features and consensus features (KERNEL/), and feature detection of small molecules
// (FEATUREFINDER/: mass traces, elution peaks, FeatureFindingMetabo).

#include "binding_utils.h"

#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace emscripten;

namespace OpenMS::WASM
{
  using MassTraceList = std::vector<MassTrace>;

  namespace
  {
    template <typename T>
    val metaValue(const T& object, const std::string& name)
    {
      return object.metaValueExists(name) ? toJS(object.getMetaValue(name)) : val::undefined();
    }

    template <typename T>
    void setMetaValue(T& object, const std::string& name, const val& value)
    {
      object.setMetaValue(name, dataValueFromJS(value));
    }

    template <typename Container>
    void checkIndex(const Container& c, Size i)
    {
      if (i >= c.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, i, c.size());
    }

    val boundingBoxToJS(const DBoundingBox<2>& bb)
    {
      val o = val::object();
      o.set("minRT", bb.minPosition()[0]);
      o.set("maxRT", bb.maxPosition()[0]);
      o.set("minMZ", bb.minPosition()[1]);
      o.set("maxMZ", bb.maxPosition()[1]);
      return o;
    }

    //------------------------------------------------------------------------------
    // Feature

    std::string featureUniqueId(const Feature& f) { return std::to_string(f.getUniqueId()); }
    void featureSetCharge(Feature& f, int charge) { f.setCharge(charge); }
    int featureGetCharge(const Feature& f) { return f.getCharge(); }
    double featureGetIntensity(const Feature& f) { return f.getIntensity(); }
    void featureSetIntensity(Feature& f, double i) { f.setIntensity(static_cast<float>(i)); }
    double featureGetQuality(const Feature& f) { return f.getOverallQuality(); }
    void featureSetQuality(Feature& f, double q) { f.setOverallQuality(static_cast<float>(q)); }
    double featureGetWidth(const Feature& f) { return f.getWidth(); }
    void featureSetWidth(Feature& f, double w) { f.setWidth(static_cast<float>(w)); }

    val featureBoundingBox(const Feature& f)
    {
      if (f.getConvexHulls().empty()) return val::null();
      return boundingBoxToJS(f.getConvexHull().getBoundingBox());
    }

    /// Bounding boxes of the mass traces (one convex hull per isotope trace).
    val featureMassTraceBoxes(const Feature& f)
    {
      return toArray(f.getConvexHulls(), [](const ConvexHull2D& h) { return boundingBoxToJS(h.getBoundingBox()); });
    }

    val featureHullPoints(const Feature& f)
    {
      val result = val::array();
      if (f.getConvexHulls().empty()) return result;
      for (const auto& p : f.getConvexHull().getHullPoints())
      {
        val point = val::array();
        point.call<void>("push", p[0]);
        point.call<void>("push", p[1]);
        result.call<void>("push", point);
      }
      return result;
    }

    val featureToJS(const Feature& f)
    {
      val o = val::object();
      o.set("uniqueId", std::to_string(f.getUniqueId()));
      o.set("rt", f.getRT());
      o.set("mz", f.getMZ());
      o.set("intensity", f.getIntensity());
      o.set("charge", f.getCharge());
      o.set("quality", f.getOverallQuality());
      o.set("width", f.getWidth());
      o.set("boundingBox", featureBoundingBox(f));
      o.set("metaValues", metaValuesToObject(f));
      return o;
    }

    //------------------------------------------------------------------------------
    // FeatureMap

    double mapSize(const FeatureMap& m) { return static_cast<double>(m.size()); }
    Feature mapGet(const FeatureMap& m, Size i) { checkIndex(m, i); return m[i]; }
    void mapSet(FeatureMap& m, Size i, const Feature& f) { checkIndex(m, i); m[i] = f; }
    void mapPush(FeatureMap& m, const Feature& f) { m.push_back(f); }
    void mapClear(FeatureMap& m, bool clear_meta_data) { m.clear(clear_meta_data); }
    void mapSortByIntensity(FeatureMap& m, bool reverse) { m.sortByIntensity(reverse); }
    void mapSortByQuality(FeatureMap& m, bool reverse) { m.sortByOverallQuality(reverse); }
    double mapMinRT(const FeatureMap& m) { return m.getMinRT(); }
    double mapMaxRT(const FeatureMap& m) { return m.getMaxRT(); }
    double mapMinMZ(const FeatureMap& m) { return m.getMinMZ(); }
    double mapMaxMZ(const FeatureMap& m) { return m.getMaxMZ(); }
    double mapMinIntensity(const FeatureMap& m) { return m.getMinIntensity(); }
    double mapMaxIntensity(const FeatureMap& m) { return m.getMaxIntensity(); }

    /// Assigns unique ids to the map and to features without one (needed before storing).
    void mapEnsureUniqueIds(FeatureMap& m)
    {
      m.ensureUniqueId();
      for (auto& f : m) f.ensureUniqueId();
    }

    val mapToArray(const FeatureMap& m) { return toArray(m, featureToJS); }

    /// The features as columns of typed arrays (for plots and tables).
    val mapFeatureTable(const FeatureMap& m)
    {
      val o = val::object();
      o.set("rt", toTypedArray<double>(m, [](const Feature& f) { return f.getRT(); }));
      o.set("mz", toTypedArray<double>(m, [](const Feature& f) { return f.getMZ(); }));
      o.set("intensity", toTypedArray<float>(m, [](const Feature& f) { return f.getIntensity(); }));
      o.set("charge", toTypedArray<int32_t>(m, [](const Feature& f) { return f.getCharge(); }));
      o.set("quality", toTypedArray<float>(m, [](const Feature& f) { return f.getOverallQuality(); }));
      o.set("width", toTypedArray<float>(m, [](const Feature& f) { return f.getWidth(); }));
      return o;
    }

    std::vector<ProteinIdentification> mapGetProteinIds(const FeatureMap& m) { return m.getProteinIdentifications(); }
    void mapSetProteinIds(FeatureMap& m, const std::vector<ProteinIdentification>& ids) { m.setProteinIdentifications(ids); }
    PeptideIdentificationList mapGetUnassignedPeptideIds(const FeatureMap& m) { return m.getUnassignedPeptideIdentifications(); }
    void mapSetUnassignedPeptideIds(FeatureMap& m, const PeptideIdentificationList& ids) { m.setUnassignedPeptideIdentifications(ids); }

    //------------------------------------------------------------------------------
    // ConsensusMap

    double consensusSize(const ConsensusMap& m) { return static_cast<double>(m.size()); }

    val consensusFeatureToJS(const ConsensusFeature& f)
    {
      val o = val::object();
      o.set("uniqueId", std::to_string(f.getUniqueId()));
      o.set("rt", f.getRT());
      o.set("mz", f.getMZ());
      o.set("intensity", f.getIntensity());
      o.set("charge", f.getCharge());
      o.set("quality", f.getQuality());
      o.set("elements", toArray(f.getFeatures(), [](const FeatureHandle& h) {
        val e = val::object();
        e.set("mapIndex", static_cast<double>(h.getMapIndex()));
        e.set("uniqueId", std::to_string(h.getUniqueId()));
        e.set("rt", h.getRT());
        e.set("mz", h.getMZ());
        e.set("intensity", h.getIntensity());
        e.set("charge", h.getCharge());
        return e;
      }));
      o.set("metaValues", metaValuesToObject(f));
      return o;
    }

    val consensusGet(const ConsensusMap& m, Size i) { checkIndex(m, i); return consensusFeatureToJS(m[i]); }
    val consensusToArray(const ConsensusMap& m) { return toArray(m, consensusFeatureToJS); }

    val consensusColumnHeaders(const ConsensusMap& m)
    {
      val result = val::array();
      for (const auto& [index, header] : m.getColumnHeaders())
      {
        val o = val::object();
        o.set("index", static_cast<double>(index));
        o.set("filename", header.filename);
        o.set("label", header.label);
        o.set("size", static_cast<double>(header.size));
        result.call<void>("push", o);
      }
      return result;
    }

    std::string consensusExperimentType(const ConsensusMap& m) { return m.getExperimentType(); }

    //------------------------------------------------------------------------------
    // Feature finding for small molecules

    double tracesSize(const MassTraceList& l) { return static_cast<double>(l.size()); }
    void tracesClear(MassTraceList& l) { l.clear(); }

    val tracesToArray(const MassTraceList& l)
    {
      return toArray(l, [](const MassTrace& t) {
        val o = val::object();
        o.set("label", t.getLabel());
        o.set("mz", t.getCentroidMZ());
        o.set("rt", t.getCentroidRT());
        o.set("mzSD", t.getCentroidSD());
        o.set("intensity", t.getIntensity(false));
        o.set("maxIntensity", t.getMaxIntensity(false));
        o.set("fwhm", t.getFWHM());
        o.set("traceLength", t.getTraceLength());
        o.set("size", static_cast<double>(t.getSize()));
        return o;
      });
    }

    void mtdRun(MassTraceDetection& d, const PeakMap& exp, MassTraceList& traces) { d.run(exp, traces); }
    void mtdRunMax(MassTraceDetection& d, const PeakMap& exp, MassTraceList& traces, Size max_traces) { d.run(exp, traces, max_traces); }
    void epdDetectPeaks(ElutionPeakDetection& d, MassTraceList& traces, MassTraceList& split) { d.detectPeaks(traces, split); }

    void ffmRun(FeatureFindingMetabo& f, MassTraceList& traces, FeatureMap& features)
    {
      std::vector<std::vector<MSChromatogram>> chromatograms;
      f.run(traces, features, chromatograms);
    }
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS;
using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_features)
{
  class_<Feature>("Feature")
    .constructor<>()
    .constructor<const Feature&>()
    .function("getUniqueId", OMS_FN(&featureUniqueId))
    .function("getRT", OMS_METHOD(Feature, &Feature::getRT))
    .function("setRT", OMS_METHOD(Feature, &Feature::setRT))
    .function("getMZ", OMS_METHOD(Feature, &Feature::getMZ))
    .function("setMZ", OMS_METHOD(Feature, &Feature::setMZ))
    .function("getIntensity", OMS_FN(&featureGetIntensity))
    .function("setIntensity", OMS_FN(&featureSetIntensity))
    .function("getCharge", OMS_FN(&featureGetCharge))
    .function("setCharge", OMS_FN(&featureSetCharge))
    .function("getOverallQuality", OMS_FN(&featureGetQuality))
    .function("setOverallQuality", OMS_FN(&featureSetQuality))
    .function("getWidth", OMS_FN(&featureGetWidth))
    .function("setWidth", OMS_FN(&featureSetWidth))
    .function("getBoundingBox", OMS_FN(&featureBoundingBox))
    .function("getMassTraceBoundingBoxes", OMS_FN(&featureMassTraceBoxes))
    .function("getConvexHullPoints", OMS_FN(&featureHullPoints))
    .function("encloses", OMS_METHOD(Feature, &Feature::encloses))
    .function("getMetaValue", OMS_FN(&metaValue<Feature>))
    .function("setMetaValue", OMS_FN(&setMetaValue<Feature>))
    .function("getMetaValues", OMS_FN(&metaValuesToObject<Feature>))
    .function("toObject", OMS_FN(&featureToJS))
    ;

  class_<FeatureMap>("FeatureMap")
    .constructor<>()
    .constructor<const FeatureMap&>()
    .function("size", OMS_FN(&mapSize))
    .function("get", OMS_FN(&mapGet))
    .function("set", OMS_FN(&mapSet))
    .function("push_back", OMS_FN(&mapPush))
    .function("clear", OMS_FN(&mapClear))
    .function("sortByRT", OMS_METHOD(FeatureMap, &FeatureMap::sortByRT))
    .function("sortByMZ", OMS_METHOD(FeatureMap, &FeatureMap::sortByMZ))
    .function("sortByPosition", OMS_METHOD(FeatureMap, &FeatureMap::sortByPosition))
    .function("sortByIntensity", OMS_FN(&mapSortByIntensity))
    .function("sortByOverallQuality", OMS_FN(&mapSortByQuality))
    .function("updateRanges", OMS_METHOD(FeatureMap, &FeatureMap::updateRanges))
    .function("getMinRT", OMS_FN(&mapMinRT))
    .function("getMaxRT", OMS_FN(&mapMaxRT))
    .function("getMinMZ", OMS_FN(&mapMinMZ))
    .function("getMaxMZ", OMS_FN(&mapMaxMZ))
    .function("getMinIntensity", OMS_FN(&mapMinIntensity))
    .function("getMaxIntensity", OMS_FN(&mapMaxIntensity))
    .function("ensureUniqueIds", OMS_FN(&mapEnsureUniqueIds))
    .function("toArray", OMS_FN(&mapToArray))
    .function("getFeatureTable", OMS_FN(&mapFeatureTable))
    .function("getProteinIdentifications", OMS_FN(&mapGetProteinIds))
    .function("setProteinIdentifications", OMS_FN(&mapSetProteinIds))
    .function("getUnassignedPeptideIdentifications", OMS_FN(&mapGetUnassignedPeptideIds))
    .function("setUnassignedPeptideIdentifications", OMS_FN(&mapSetUnassignedPeptideIds))
    .function("getMetaValue", OMS_FN(&metaValue<FeatureMap>))
    .function("setMetaValue", OMS_FN(&setMetaValue<FeatureMap>))
    .function("getMetaValues", OMS_FN(&metaValuesToObject<FeatureMap>))
    ;

  class_<ConsensusMap>("ConsensusMap")
    .constructor<>()
    .constructor<const ConsensusMap&>()
    .function("size", OMS_FN(&consensusSize))
    .function("get", OMS_FN(&consensusGet))
    .function("toArray", OMS_FN(&consensusToArray))
    .function("getColumnHeaders", OMS_FN(&consensusColumnHeaders))
    .function("getExperimentType", OMS_FN(&consensusExperimentType))
    .function("updateRanges", OMS_METHOD(ConsensusMap, &ConsensusMap::updateRanges))
    .function("sortByRT", OMS_METHOD(ConsensusMap, &ConsensusMap::sortByRT))
    .function("sortByMZ", OMS_METHOD(ConsensusMap, &ConsensusMap::sortByMZ))
    ;

  class_<MassTraceList>("MassTraceList")
    .constructor<>()
    .function("size", OMS_FN(&tracesSize))
    .function("clear", OMS_FN(&tracesClear))
    .function("toArray", OMS_FN(&tracesToArray))
    ;

  class_<MassTraceDetection, base<DefaultParamHandler>>("MassTraceDetection")
    .constructor<>()
    .function("run", OMS_FN(&mtdRun))
    .function("run", OMS_FN(&mtdRunMax))
    ;

  class_<ElutionPeakDetection, base<DefaultParamHandler>>("ElutionPeakDetection")
    .constructor<>()
    .function("detectPeaks", OMS_FN(&epdDetectPeaks))
    ;

  class_<FeatureFindingMetabo, base<DefaultParamHandler>>("FeatureFindingMetabo")
    .constructor<>()
    .function("run", OMS_FN(&ffmRun))
    ;
}
