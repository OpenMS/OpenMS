// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Parameters of algorithms (Param, DefaultParamHandler), signal processing (PROCESSING/)
// and spectrum comparison (COMPARISON/).

#include "binding_utils.h"

#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SpectrumAlignmentScore.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>

using namespace emscripten;

namespace OpenMS::WASM
{
  namespace
  {
    //------------------------------------------------------------------------------
    // Param

    JS::ParamValue paramGetValue(const Param& p, const std::string& key)
    {
      const Param::ParamEntry& entry = p.getEntry(key); // throws ElementNotFound
      return toJS(entry.value, entry.isBool());
    }

    void paramSetValue(Param& p, const std::string& key, const JS::ParamValue& value)
    {
      if (p.exists(key))
      {
        const Param::ParamEntry& entry = p.getEntry(key);
        p.setValue(key, paramValueFromJS(value, entry.value.valueType()), entry.description,
                   std::vector<std::string>(entry.tags.begin(), entry.tags.end()));
      }
      else
      {
        p.setValue(key, paramValueFromJS(value));
      }
    }

    void paramSetValueDescribed(Param& p, const std::string& key, const JS::ParamValue& value, const std::string& description)
    {
      p.setValue(key, paramValueFromJS(value), description);
    }

    double paramSize(const Param& p) { return static_cast<double>(p.size()); }

    JS::StringArray paramKeys(const Param& p)
    {
      val result = val::array();
      for (auto it = p.begin(); it != p.end(); ++it) result.call<void>("push", it.getName());
      return result;
    }

    std::string valueTypeName(ParamValue::ValueType t)
    {
      switch (t)
      {
        case ParamValue::STRING_VALUE: return "string";
        case ParamValue::INT_VALUE: return "int";
        case ParamValue::DOUBLE_VALUE: return "double";
        case ParamValue::STRING_LIST: return "string[]";
        case ParamValue::INT_LIST: return "int[]";
        case ParamValue::DOUBLE_LIST: return "double[]";
        default: return "empty";
      }
    }

    JS::ParamEntry entryToJS(const std::string& name, const Param::ParamEntry& e)
    {
      val o = val::object();
      o.set("name", name);
      o.set("value", toJS(e.value, e.isBool()));
      o.set("type", e.isBool() ? std::string("bool") : valueTypeName(e.value.valueType()));
      o.set("description", e.description);
      o.set("tags", toArray(e.tags));
      if (!e.valid_strings.empty()) o.set("validStrings", toArray(e.valid_strings));
      if (e.value.valueType() == ParamValue::INT_VALUE || e.value.valueType() == ParamValue::INT_LIST)
      {
        if (e.min_int != -std::numeric_limits<int>::max()) o.set("min", e.min_int);
        if (e.max_int != std::numeric_limits<int>::max()) o.set("max", e.max_int);
      }
      if (e.value.valueType() == ParamValue::DOUBLE_VALUE || e.value.valueType() == ParamValue::DOUBLE_LIST)
      {
        if (e.min_float != -std::numeric_limits<double>::max()) o.set("min", e.min_float);
        if (e.max_float != std::numeric_limits<double>::max()) o.set("max", e.max_float);
      }
      return o;
    }

    /// All entries with their documentation: what a user interface needs to build a parameter form.
    JS::ParamEntries paramEntries(const Param& p)
    {
      val result = val::array();
      for (auto it = p.begin(); it != p.end(); ++it) result.call<void>("push", entryToJS(it.getName(), *it));
      return result;
    }

    JS::ParamEntry paramGetEntry(const Param& p, const std::string& key) { return entryToJS(key, p.getEntry(key)); }
    std::string paramGetDescription(const Param& p, const std::string& key) { return p.getDescription(key); }
    JS::StringArray paramGetTags(const Param& p, const std::string& key) { return toArray(p.getTags(key)); }
    Param paramCopy(const Param& p, const std::string& prefix, bool remove_prefix) { return p.copy(prefix, remove_prefix); }
    void paramInsert(Param& p, const std::string& prefix, const Param& other) { p.insert(prefix, other); }

    //------------------------------------------------------------------------------
    // DefaultParamHandler

    Param handlerGetParameters(const DefaultParamHandler& h) { return h.getParameters(); }
    Param handlerGetDefaults(const DefaultParamHandler& h) { return h.getDefaults(); }
    std::string handlerGetName(const DefaultParamHandler& h) { return h.getName(); }
    JS::ParamValues handlerGetParameterValues(const DefaultParamHandler& h) { return paramToObject(h.getParameters()); }

    /// Sets the parameters named in @p values ({name: value}), keeping the others.
    void handlerSetParameterValues(DefaultParamHandler& h, const JS::ParamValues& values)
    {
      Param p = h.getParameters();
      updateParamFromObject(p, values);
      h.setParameters(p);
    }

    //------------------------------------------------------------------------------
    // Algorithms

    void pickerPickChromatogram(const PeakPickerHiRes& p, const MSChromatogram& in, MSChromatogram& out) { p.pick(in, out); }
    void pickerPickSpectrum(const PeakPickerHiRes& p, const MSSpectrum& in, MSSpectrum& out) { p.pick(in, out); }
    void pickerPickExperiment(const PeakPickerHiRes& p, const PeakMap& in, PeakMap& out) { p.pickExperiment(in, out); }
    void pickerPickExperimentChecked(const PeakPickerHiRes& p, const PeakMap& in, PeakMap& out, bool check_spectrum_type) { p.pickExperiment(in, out, check_spectrum_type); }

    template <typename Filter> void filterSpectrum(Filter& f, MSSpectrum& s) { f.filter(s); }
    template <typename Filter> void filterChromatogram(Filter& f, MSChromatogram& c) { f.filter(c); }
    template <typename Filter> void filterExperiment(Filter& f, PeakMap& e) { f.filterExperiment(e); }

    template <typename Filter> void filterPeakSpectrum(Filter& f, MSSpectrum& s) { f.filterPeakSpectrum(s); }
    template <typename Filter> void filterPeakMap(Filter& f, PeakMap& e) { f.filterPeakMap(e); }
    void normalizerSpectrum(const Normalizer& f, MSSpectrum& s) { f.filterPeakSpectrum(s); }
    void normalizerPeakMap(const Normalizer& f, PeakMap& e) { f.filterPeakMap(e); }

    JS::Alignment alignmentGet(const SpectrumAlignment& a, const MSSpectrum& s1, const MSSpectrum& s2)
    {
      std::vector<std::pair<Size, Size>> alignment;
      a.getSpectrumAlignment(alignment, s1, s2);
      return toArray(alignment, [](const std::pair<Size, Size>& p) {
        val pair = val::array();
        pair.call<void>("push", static_cast<double>(p.first));
        pair.call<void>("push", static_cast<double>(p.second));
        return pair;
      });
    }

    double alignmentScore(const SpectrumAlignmentScore& s, const MSSpectrum& s1, const MSSpectrum& s2) { return s(s1, s2); }
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS;
using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_processing)
{
  class_<Param>("Param")
    .constructor<>()
    .constructor<const Param&>()
    .function("getValue", OMS_FN(&paramGetValue))
    .function("setValue", OMS_FN(&paramSetValue))
    .function("setValue", OMS_FN(&paramSetValueDescribed))
    .function("exists", OMS_METHOD(Param, &Param::exists))
    .function("remove", OMS_METHOD(Param, &Param::remove))
    .function("size", OMS_FN(&paramSize))
    .function("empty", OMS_METHOD(Param, &Param::empty))
    .function("clear", OMS_METHOD(Param, &Param::clear))
    .function("keys", OMS_FN(&paramKeys))
    .function("getEntry", OMS_FN(&paramGetEntry))
    .function("getEntries", OMS_FN(&paramEntries))
    .function("getDescription", OMS_FN(&paramGetDescription))
    .function("getTags", OMS_FN(&paramGetTags))
    .function("copy", OMS_FN(&paramCopy))
    .function("insert", OMS_FN(&paramInsert))
    .function("toObject", OMS_FN(&paramToObject))
    .function("update", OMS_FN(&updateParamFromObject))
    ;

  class_<DefaultParamHandler>("DefaultParamHandler")
    .function("getName", OMS_FN(&handlerGetName))
    .function("getParameters", OMS_FN(&handlerGetParameters))
    .function("setParameters", OMS_METHOD(DefaultParamHandler, &DefaultParamHandler::setParameters))
    .function("getDefaults", OMS_FN(&handlerGetDefaults))
    .function("getParameterValues", OMS_FN(&handlerGetParameterValues))
    .function("setParameterValues", OMS_FN(&handlerSetParameterValues))
    ;

  class_<PeakPickerHiRes, base<DefaultParamHandler>>("PeakPickerHiRes")
    .constructor<>()
    .function("pick", OMS_FN(&pickerPickSpectrum))
    .function("pickChromatogram", OMS_FN(&pickerPickChromatogram))
    .function("pickExperiment", OMS_FN(&pickerPickExperiment))
    .function("pickExperiment", OMS_FN(&pickerPickExperimentChecked))
    ;

  class_<GaussFilter, base<DefaultParamHandler>>("GaussFilter")
    .constructor<>()
    .function("filter", OMS_FN(&filterSpectrum<GaussFilter>))
    .function("filterChromatogram", OMS_FN(&filterChromatogram<GaussFilter>))
    .function("filterExperiment", OMS_FN(&filterExperiment<GaussFilter>))
    ;

  class_<SavitzkyGolayFilter, base<DefaultParamHandler>>("SavitzkyGolayFilter")
    .constructor<>()
    .function("filter", OMS_FN(&filterSpectrum<SavitzkyGolayFilter>))
    .function("filterChromatogram", OMS_FN(&filterChromatogram<SavitzkyGolayFilter>))
    .function("filterExperiment", OMS_FN(&filterExperiment<SavitzkyGolayFilter>))
    ;

  class_<Normalizer, base<DefaultParamHandler>>("Normalizer")
    .constructor<>()
    .function("filterPeakSpectrum", OMS_FN(&normalizerSpectrum))
    .function("filterPeakMap", OMS_FN(&normalizerPeakMap))
    ;

  class_<ThresholdMower, base<DefaultParamHandler>>("ThresholdMower")
    .constructor<>()
    .function("filterPeakSpectrum", OMS_FN(&filterPeakSpectrum<ThresholdMower>))
    .function("filterPeakMap", OMS_FN(&filterPeakMap<ThresholdMower>))
    ;

  class_<WindowMower, base<DefaultParamHandler>>("WindowMower")
    .constructor<>()
    .function("filterPeakSpectrum", OMS_FN(&filterPeakSpectrum<WindowMower>))
    .function("filterPeakMap", OMS_FN(&filterPeakMap<WindowMower>))
    ;

  class_<NLargest, base<DefaultParamHandler>>("NLargest")
    .constructor<>()
    .constructor<UInt>()
    .function("filterPeakSpectrum", OMS_FN(&filterPeakSpectrum<NLargest>))
    .function("filterPeakMap", OMS_FN(&filterPeakMap<NLargest>))
    ;

  class_<SpectrumAlignment, base<DefaultParamHandler>>("SpectrumAlignment")
    .constructor<>()
    .function("getSpectrumAlignment", OMS_FN(&alignmentGet))
    ;

  class_<SpectrumAlignmentScore, base<DefaultParamHandler>>("SpectrumAlignmentScore")
    .constructor<>()
    .function("score", OMS_FN(&alignmentScore))
    ;
}
