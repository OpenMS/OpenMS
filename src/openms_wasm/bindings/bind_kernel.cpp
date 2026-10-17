// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Spectra, chromatograms and experiments (KERNEL/).

#include "binding_utils.h"

#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Product.h>

#include <algorithm>

using namespace emscripten;

namespace OpenMS::WASM
{
  namespace
  {
    //------------------------------------------------------------------------------
    // Plain JavaScript objects for small records (no .delete() needed)

    JS::PrecursorObject precursorToJS(const Precursor& p)
    {
      val o = val::object();
      o.set("mz", p.getMZ());
      o.set("intensity", p.getIntensity());
      o.set("charge", p.getCharge());
      o.set("isolationWindowLowerOffset", p.getIsolationWindowLowerOffset());
      o.set("isolationWindowUpperOffset", p.getIsolationWindowUpperOffset());
      o.set("activationEnergy", p.getActivationEnergy());
      o.set("activationMethods", toArray(p.getActivationMethodsAsShortString()));
      o.set("driftTime", p.getDriftTime());
      return o;
    }

    Precursor precursorFromJS(const val& o)
    {
      Precursor p;
      const auto number = [&o](const char* key, auto setter) {
        const val v = o[key];
        if (!isNullish(v)) setter(v.as<double>());
      };
      number("mz", [&p](double v) { p.setMZ(v); });
      number("intensity", [&p](double v) { p.setIntensity(static_cast<float>(v)); });
      number("charge", [&p](double v) { p.setCharge(static_cast<Int>(v)); });
      number("isolationWindowLowerOffset", [&p](double v) { p.setIsolationWindowLowerOffset(v); });
      number("isolationWindowUpperOffset", [&p](double v) { p.setIsolationWindowUpperOffset(v); });
      number("activationEnergy", [&p](double v) { p.setActivationEnergy(v); });
      number("driftTime", [&p](double v) { p.setDriftTime(v); });
      const val methods = o["activationMethods"];
      if (!isNullish(methods))
      {
        std::set<Precursor::ActivationMethod> am;
        for (const auto& name : stringsFromJS(methods))
        {
          am.insert(Precursor::toActivationMethod(name)); // throws for unknown names
        }
        p.setActivationMethods(am);
      }
      return p;
    }

    JS::ProductObject productToJS(const Product& p)
    {
      val o = val::object();
      o.set("mz", p.getMZ());
      o.set("isolationWindowLowerOffset", p.getIsolationWindowLowerOffset());
      o.set("isolationWindowUpperOffset", p.getIsolationWindowUpperOffset());
      return o;
    }

    template <typename DataArrays>
    val dataArraysToJS(const DataArrays& arrays)
    {
      val result = val::array();
      for (const auto& a : arrays)
      {
        val o = val::object();
        o.set("name", a.getName());
        using Element = typename std::decay_t<decltype(a)>::value_type;
        if constexpr (std::is_same_v<Element, std::string>)
        {
          o.set("data", toArray(a));
        }
        else
        {
          o.set("data", toTypedArray(a.data(), a.size()));
        }
        result.call<void>("push", o);
      }
      return result;
    }

    //------------------------------------------------------------------------------
    // MSSpectrum

    JS::SpectrumPeaks spectrumGetPeaks(const MSSpectrum& s)
    {
      val o = val::object();
      o.set("mz", toTypedArray<double>(s, [](const Peak1D& p) { return p.getMZ(); }));
      o.set("intensity", toTypedArray<float>(s, [](const Peak1D& p) { return p.getIntensity(); }));
      return o;
    }

    void spectrumSetPeaks(MSSpectrum& s, const JS::NumberInput& mz, const JS::NumberInput& intensity)
    {
      const auto mzs = numbersFromJS<double>(mz);
      const auto ints = numbersFromJS<float>(intensity);
      if (mzs.size() != ints.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "m/z and intensity arrays differ in length (" + std::to_string(mzs.size()) + " vs. " + std::to_string(ints.size()) + ")");
      }
      s.resize(mzs.size());
      for (Size i = 0; i < mzs.size(); ++i)
      {
        s[i].setMZ(mzs[i]);
        s[i].setIntensity(ints[i]);
      }
      // the old data arrays no longer match the peaks
      s.getFloatDataArrays().clear();
      s.getIntegerDataArrays().clear();
      s.getStringDataArrays().clear();
    }

    JS::PeakObject spectrumGetPeak(const MSSpectrum& s, Size i)
    {
      if (i >= s.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, i, s.size());
      val o = val::object();
      o.set("mz", s[i].getMZ());
      o.set("intensity", s[i].getIntensity());
      return o;
    }

    void spectrumAddPeak(MSSpectrum& s, double mz, double intensity)
    {
      s.emplace_back(mz, static_cast<float>(intensity));
    }

    JS::PrecursorList spectrumGetPrecursors(const MSSpectrum& s) { return toArray(s.getPrecursors(), precursorToJS); }

    void spectrumSetPrecursors(MSSpectrum& s, const JS::PrecursorInputList& precursors)
    {
      std::vector<Precursor> result;
      const unsigned n = precursors["length"].as<unsigned>();
      for (unsigned i = 0; i < n; ++i)
      {
        result.push_back(precursorFromJS(precursors[i]));
      }
      s.setPrecursors(result);
    }

    JS::ProductList spectrumGetProducts(const MSSpectrum& s) { return toArray(s.getProducts(), productToJS); }

    std::string spectrumGetType(const MSSpectrum& s) { return SpectrumSettings::spectrumTypeToString(s.getType()); }
    void spectrumSetType(MSSpectrum& s, const std::string& type) { s.setType(SpectrumSettings::toSpectrumType(type)); }
    std::string spectrumGuessType(const MSSpectrum& s) { return SpectrumSettings::spectrumTypeToString(s.getType(true)); }

    std::string spectrumGetPolarity(const MSSpectrum& s) { return IonSource::polarityToString(s.getInstrumentSettings().getPolarity()); }
    void spectrumSetPolarity(MSSpectrum& s, const std::string& p) { s.getInstrumentSettings().setPolarity(IonSource::toPolarity(p)); }

    void spectrumClear(MSSpectrum& s, bool clear_meta_data) { s.clear(clear_meta_data); }
    // (size() and empty() come from MSSpectrum's private std::vector base)
    double spectrumSize(const MSSpectrum& s) { return static_cast<double>(s.size()); }
    bool spectrumEmpty(const MSSpectrum& s) { return s.empty(); }
    void spectrumSortByIntensity(MSSpectrum& s, bool reverse) { s.sortByIntensity(reverse); }
    bool spectrumIsSorted(const MSSpectrum& s) { return s.isSorted(); }
    double spectrumFindNearest(const MSSpectrum& s, double mz)
    {
      if (s.empty()) return -1;
      return static_cast<double>(s.findNearest(mz));
    }
    int spectrumFindNearestTol(const MSSpectrum& s, double mz, double tolerance) { return s.findNearest(mz, tolerance); }
    double spectrumTIC(const MSSpectrum& s) { return s.calculateTIC(); }

    JS::BasePeak spectrumBasePeak(const MSSpectrum& s)
    {
      if (s.empty()) return val::null();
      const auto it = s.getBasePeak();
      val o = val::object();
      o.set("index", static_cast<double>(it - s.begin()));
      o.set("mz", it->getMZ());
      o.set("intensity", it->getIntensity());
      return o;
    }

    JS::FloatDataArrays spectrumFloatDataArrays(const MSSpectrum& s) { return dataArraysToJS(s.getFloatDataArrays()); }
    JS::IntegerDataArrays spectrumIntegerDataArrays(const MSSpectrum& s) { return dataArraysToJS(s.getIntegerDataArrays()); }
    JS::StringDataArrays spectrumStringDataArrays(const MSSpectrum& s) { return dataArraysToJS(s.getStringDataArrays()); }

    //------------------------------------------------------------------------------
    // MSChromatogram

    JS::ChromatogramPeaks chromatogramGetPeaks(const MSChromatogram& c)
    {
      val o = val::object();
      o.set("rt", toTypedArray<double>(c, [](const ChromatogramPeak& p) { return p.getRT(); }));
      o.set("intensity", toTypedArray<float>(c, [](const ChromatogramPeak& p) { return p.getIntensity(); }));
      return o;
    }

    void chromatogramSetPeaks(MSChromatogram& c, const JS::NumberInput& rt, const JS::NumberInput& intensity)
    {
      const auto rts = numbersFromJS<double>(rt);
      const auto ints = numbersFromJS<float>(intensity);
      if (rts.size() != ints.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "RT and intensity arrays differ in length (" + std::to_string(rts.size()) + " vs. " + std::to_string(ints.size()) + ")");
      }
      c.resize(rts.size());
      for (Size i = 0; i < rts.size(); ++i)
      {
        c[i].setRT(rts[i]);
        c[i].setIntensity(ints[i]);
      }
      c.getFloatDataArrays().clear();
      c.getIntegerDataArrays().clear();
      c.getStringDataArrays().clear();
    }

    double chromatogramSize(const MSChromatogram& c) { return static_cast<double>(c.size()); }
    bool chromatogramEmpty(const MSChromatogram& c) { return c.empty(); }

    void chromatogramAddPeak(MSChromatogram& c, double rt, double intensity)
    {
      c.emplace_back(rt, static_cast<float>(intensity));
    }

    JS::PrecursorObject chromatogramGetPrecursor(const MSChromatogram& c) { return precursorToJS(c.getPrecursor()); }
    void chromatogramSetPrecursor(MSChromatogram& c, const JS::PrecursorInput& p) { c.setPrecursor(precursorFromJS(p)); }
    JS::ProductObject chromatogramGetProduct(const MSChromatogram& c) { return productToJS(c.getProduct()); }
    std::string chromatogramGetType(const MSChromatogram& c) { return ChromatogramSettings::ChromatogramNames[static_cast<int>(c.getChromatogramType())]; }
    void chromatogramClear(MSChromatogram& c, bool clear_meta_data) { c.clear(clear_meta_data); }
    double chromatogramTIC(const MSChromatogram& c)
    {
      double sum = 0;
      for (const auto& p : c) sum += p.getIntensity();
      return sum;
    }
    JS::FloatDataArrays chromatogramFloatDataArrays(const MSChromatogram& c) { return dataArraysToJS(c.getFloatDataArrays()); }

    //------------------------------------------------------------------------------
    // MSExperiment

    void checkSpectrumIndex(const MSExperiment& e, Size i)
    {
      if (i >= e.getNrSpectra()) throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, i, e.getNrSpectra());
    }

    void checkChromatogramIndex(const MSExperiment& e, Size i)
    {
      if (i >= e.getNrChromatograms()) throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, i, e.getNrChromatograms());
    }

    MSSpectrum experimentGetSpectrum(const MSExperiment& e, Size i)
    {
      checkSpectrumIndex(e, i);
      return e.getSpectrum(i);
    }

    void experimentSetSpectrum(MSExperiment& e, Size i, const MSSpectrum& s)
    {
      checkSpectrumIndex(e, i);
      e.getSpectrum(i) = s;
    }

    MSChromatogram experimentGetChromatogram(const MSExperiment& e, Size i)
    {
      checkChromatogramIndex(e, i);
      return e.getChromatogram(i);
    }

    void experimentSetChromatogram(MSExperiment& e, Size i, const MSChromatogram& c)
    {
      checkChromatogramIndex(e, i);
      e.getChromatogram(i) = c;
    }

    void experimentAddSpectrum(MSExperiment& e, const MSSpectrum& s) { e.addSpectrum(s); }
    void experimentAddChromatogram(MSExperiment& e, const MSChromatogram& c) { e.addChromatogram(c); }
    JS::NumberArray experimentGetMSLevels(const MSExperiment& e) { return toArray(e.getMSLevels()); }
    void experimentSortSpectra(MSExperiment& e, bool sort_mz) { e.sortSpectra(sort_mz); }
    void experimentSortChromatograms(MSExperiment& e, bool sort_rt) { e.sortChromatograms(sort_rt); }
    bool experimentIsSorted(const MSExperiment& e, bool check_mz) { return e.isSorted(check_mz); }
    void experimentClear(MSExperiment& e, bool clear_meta_data) { e.clear(clear_meta_data); }
    MSChromatogram experimentTIC(const MSExperiment& e, double rt_bin_size, UInt ms_level) { return e.calculateTIC(static_cast<float>(rt_bin_size), ms_level); }
    double experimentPeakCount(const MSExperiment& e) { return static_cast<double>(e.getSize()); }

    double experimentClosestSpectrum(const MSExperiment& e, double rt, UInt ms_level)
    {
      if (e.empty()) return -1;
      const auto it = (ms_level == 0) ? e.getClosestSpectrumInRT(rt) : e.getClosestSpectrumInRT(rt, ms_level);
      return it == e.end() ? -1 : static_cast<double>(it - e.begin());
    }

    int experimentPrecursorSpectrum(const MSExperiment& e, int index) { return e.getPrecursorSpectrum(index); }
    int experimentFirstProductSpectrum(const MSExperiment& e, int index) { return e.getFirstProductSpectrum(index); }

    /// One row per spectrum: what a spectrum table of a user interface shows, without copying peaks.
    JS::SpectrumSummary experimentSpectrumSummary(const MSExperiment& e)
    {
      val rows = val::array();
      for (Size i = 0; i < e.size(); ++i)
      {
        const MSSpectrum& s = e[i];
        val o = val::object();
        o.set("index", static_cast<double>(i));
        o.set("nativeID", s.getNativeID());
        o.set("rt", s.getRT());
        o.set("msLevel", s.getMSLevel());
        o.set("size", static_cast<double>(s.size()));
        o.set("tic", s.calculateTIC());
        if (!s.empty())
        {
          const auto bp = s.getBasePeak();
          o.set("basePeakMZ", bp->getMZ());
          o.set("basePeakIntensity", bp->getIntensity());
          // not front()/back(): spectra need not be sorted by m/z (setPeaks(), addPeak(), MGF files)
          const auto [min_peak, max_peak] = std::minmax_element(s.begin(), s.end(), Peak1D::MZLess());
          o.set("minMZ", min_peak->getMZ());
          o.set("maxMZ", max_peak->getMZ());
        }
        if (!s.getPrecursors().empty())
        {
          o.set("precursorMZ", s.getPrecursors()[0].getMZ());
          o.set("precursorCharge", s.getPrecursors()[0].getCharge());
        }
        rows.call<void>("push", o);
      }
      return rows;
    }

    /// All peaks of @p ms_level (0: all levels) in an RT/m/z window as three columns (for heat maps etc.).
    JS::PeakData2D experimentGet2DPeakData(const MSExperiment& e, double min_rt, double max_rt, double min_mz, double max_mz, UInt ms_level)
    {
      std::vector<double> rt, mz;
      std::vector<float> intensity;
      for (auto it = e.RTBegin(min_rt); it != e.end() && it->getRT() <= max_rt; ++it)
      {
        if (ms_level != 0 && it->getMSLevel() != ms_level) continue;
        for (auto p = it->MZBegin(min_mz); p != it->end() && p->getMZ() <= max_mz; ++p)
        {
          rt.push_back(it->getRT());
          mz.push_back(p->getMZ());
          intensity.push_back(p->getIntensity());
        }
      }
      val o = val::object();
      o.set("rt", toTypedArray(rt));
      o.set("mz", toTypedArray(mz));
      o.set("intensity", toTypedArray(intensity));
      return o;
    }
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS;
using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_kernel)
{
  class_<MSSpectrum>("MSSpectrum")
    .constructor<>()
    .constructor<const MSSpectrum&>()
    .function("size", OMS_FN(&spectrumSize))
    .function("empty", OMS_FN(&spectrumEmpty))
    .function("clear", OMS_FN(&spectrumClear))
    .function("getRT", OMS_METHOD(MSSpectrum, &MSSpectrum::getRT))
    .function("setRT", OMS_METHOD(MSSpectrum, &MSSpectrum::setRT))
    .function("getMSLevel", OMS_METHOD(MSSpectrum, &MSSpectrum::getMSLevel))
    .function("setMSLevel", OMS_METHOD(MSSpectrum, &MSSpectrum::setMSLevel))
    .function("getName", OMS_METHOD(MSSpectrum, &MSSpectrum::getName))
    .function("setName", OMS_METHOD(MSSpectrum, &MSSpectrum::setName))
    .function("getNativeID", OMS_METHOD(MSSpectrum, &SpectrumSettings::getNativeID))
    .function("setNativeID", OMS_METHOD(MSSpectrum, &SpectrumSettings::setNativeID))
    .function("getDriftTime", OMS_METHOD(MSSpectrum, &MSSpectrum::getDriftTime))
    .function("setDriftTime", OMS_METHOD(MSSpectrum, &MSSpectrum::setDriftTime))
    .function("getType", OMS_FN(&spectrumGetType))
    .function("setType", OMS_FN(&spectrumSetType))
    .function("guessType", OMS_FN(&spectrumGuessType))
    .function("getPolarity", OMS_FN(&spectrumGetPolarity))
    .function("setPolarity", OMS_FN(&spectrumSetPolarity))
    .function("getPeaks", OMS_FN(&spectrumGetPeaks))
    .function("setPeaks", OMS_FN(&spectrumSetPeaks))
    .function("getPeak", OMS_FN(&spectrumGetPeak))
    .function("addPeak", OMS_FN(&spectrumAddPeak))
    .function("getPrecursors", OMS_FN(&spectrumGetPrecursors))
    .function("setPrecursors", OMS_FN(&spectrumSetPrecursors))
    .function("getProducts", OMS_FN(&spectrumGetProducts))
    .function("sortByPosition", OMS_METHOD(MSSpectrum, &MSSpectrum::sortByPosition))
    .function("sortByIntensity", OMS_FN(&spectrumSortByIntensity))
    .function("isSorted", OMS_FN(&spectrumIsSorted))
    .function("findNearest", OMS_FN(&spectrumFindNearest))
    .function("findNearest", OMS_FN(&spectrumFindNearestTol))
    .function("calculateTIC", OMS_FN(&spectrumTIC))
    .function("getBasePeak", OMS_FN(&spectrumBasePeak))
    .function("getFloatDataArrays", OMS_FN(&spectrumFloatDataArrays))
    .function("getIntegerDataArrays", OMS_FN(&spectrumIntegerDataArrays))
    .function("getStringDataArrays", OMS_FN(&spectrumStringDataArrays))
    .function("getMetaValue", OMS_FN(&metaValue<MSSpectrum>))
    .function("setMetaValue", OMS_FN(&setMetaValue<MSSpectrum>))
    .function("metaValueExists", OMS_FN(&metaValueExists<MSSpectrum>))
    .function("getMetaValueKeys", OMS_FN(&metaValueKeys<MSSpectrum>))
    .function("getMetaValues", OMS_FN(&metaValuesToObject<MSSpectrum>))
    ;

  class_<MSChromatogram>("MSChromatogram")
    .constructor<>()
    .constructor<const MSChromatogram&>()
    .function("size", OMS_FN(&chromatogramSize))
    .function("empty", OMS_FN(&chromatogramEmpty))
    .function("clear", OMS_FN(&chromatogramClear))
    .function("getName", OMS_METHOD(MSChromatogram, &MSChromatogram::getName))
    .function("setName", OMS_METHOD(MSChromatogram, &MSChromatogram::setName))
    .function("getNativeID", OMS_METHOD(MSChromatogram, &ChromatogramSettings::getNativeID))
    .function("setNativeID", OMS_METHOD(MSChromatogram, &ChromatogramSettings::setNativeID))
    .function("getMZ", OMS_METHOD(MSChromatogram, &MSChromatogram::getMZ))
    .function("getChromatogramType", OMS_FN(&chromatogramGetType))
    .function("getPrecursor", OMS_FN(&chromatogramGetPrecursor))
    .function("setPrecursor", OMS_FN(&chromatogramSetPrecursor))
    .function("getProduct", OMS_FN(&chromatogramGetProduct))
    .function("getPeaks", OMS_FN(&chromatogramGetPeaks))
    .function("setPeaks", OMS_FN(&chromatogramSetPeaks))
    .function("addPeak", OMS_FN(&chromatogramAddPeak))
    .function("sortByPosition", OMS_METHOD(MSChromatogram, &MSChromatogram::sortByPosition))
    .function("isSorted", OMS_METHOD(MSChromatogram, &MSChromatogram::isSorted))
    .function("calculateTIC", OMS_FN(&chromatogramTIC))
    .function("getFloatDataArrays", OMS_FN(&chromatogramFloatDataArrays))
    .function("getMetaValue", OMS_FN(&metaValue<MSChromatogram>))
    .function("setMetaValue", OMS_FN(&setMetaValue<MSChromatogram>))
    .function("getMetaValueKeys", OMS_FN(&metaValueKeys<MSChromatogram>))
    .function("getMetaValues", OMS_FN(&metaValuesToObject<MSChromatogram>))
    ;

  class_<MSExperiment>("MSExperiment")
    .constructor<>()
    .constructor<const MSExperiment&>()
    .function("size", OMS_METHOD(MSExperiment, &MSExperiment::size))
    .function("empty", OMS_METHOD(MSExperiment, &MSExperiment::empty))
    .function("getNrSpectra", OMS_METHOD(MSExperiment, &MSExperiment::getNrSpectra))
    .function("getNrChromatograms", OMS_METHOD(MSExperiment, &MSExperiment::getNrChromatograms))
    .function("getNrPeaks", OMS_FN(&experimentPeakCount))
    .function("getSpectrum", OMS_FN(&experimentGetSpectrum))
    .function("setSpectrum", OMS_FN(&experimentSetSpectrum))
    .function("addSpectrum", OMS_FN(&experimentAddSpectrum))
    .function("getChromatogram", OMS_FN(&experimentGetChromatogram))
    .function("setChromatogram", OMS_FN(&experimentSetChromatogram))
    .function("addChromatogram", OMS_FN(&experimentAddChromatogram))
    .function("getMSLevels", OMS_FN(&experimentGetMSLevels))
    .function("updateRanges", OMS_METHOD(MSExperiment, &MSExperiment::updateRanges))
    .function("getMinRT", OMS_METHOD(MSExperiment, &MSExperiment::getMinRT))
    .function("getMaxRT", OMS_METHOD(MSExperiment, &MSExperiment::getMaxRT))
    .function("getMinMZ", OMS_METHOD(MSExperiment, &MSExperiment::getMinMZ))
    .function("getMaxMZ", OMS_METHOD(MSExperiment, &MSExperiment::getMaxMZ))
    .function("getMinIntensity", OMS_METHOD(MSExperiment, &MSExperiment::getMinIntensity))
    .function("getMaxIntensity", OMS_METHOD(MSExperiment, &MSExperiment::getMaxIntensity))
    .function("sortSpectra", OMS_FN(&experimentSortSpectra))
    .function("sortChromatograms", OMS_FN(&experimentSortChromatograms))
    .function("isSorted", OMS_FN(&experimentIsSorted))
    .function("clear", OMS_FN(&experimentClear))
    .function("calculateTIC", OMS_FN(&experimentTIC))
    .function("getClosestSpectrumInRT", OMS_FN(&experimentClosestSpectrum))
    .function("getPrecursorSpectrum", OMS_FN(&experimentPrecursorSpectrum))
    .function("getFirstProductSpectrum", OMS_FN(&experimentFirstProductSpectrum))
    .function("getSpectrumSummary", OMS_FN(&experimentSpectrumSummary))
    .function("get2DPeakData", OMS_FN(&experimentGet2DPeakData))
    ;
}
