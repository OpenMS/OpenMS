// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: OpenMS-DIAtracer project $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <functional>
#include <cstdio>
#include <fstream>
#include <atomic>
#include <thread>
#include <map>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_DIAPseudoSpectra DIAPseudoSpectra

  @brief Extracts pseudo-MS/MS ("pseudo-DDA") spectra from diaPASEF (ion-mobility DIA) data.

  Reconstructs DDA-like MS2 spectra from ion-mobility DIA (timsTOF / diaPASEF) data by
  exploiting the fact that, in PASEF, fragments are recorded at the ion-mobility (1/K0)
  elution coordinate of their precursor. Precursor and fragment mass traces are detected
  IM-aware (@ref MassTraceDetection with an ion-mobility tolerance), then fragments are
  assigned to a precursor when they co-localize in ion mobility, retention time, and
  elution-profile (Pearson) correlation. The result is written as searchable mzML.

  This is the OpenMS-native, BSD-3 analogue of Nesvilab diaTracer. It is a discovery
  front-end (enables open / semi-tryptic / PTM searches); it does not perform the database
  search, FDR, library building or quantification.

  See PLAN.md / REVIEW.md in the project workspace for the full design and adversarial review.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_DIAPseudoSpectra.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_DIAPseudoSpectra.html
*/

/// @cond TOPPCLASSES

namespace
{
  /// A detected mass trace reduced to the quantities we correlate/gate on.
  /// Current + peak resident set, read from /proc/self/status (Linux). Logged at phase boundaries so
  /// the memory peak is MEASURED rather than estimated: it tells us whether the peak sits in the
  /// raw-load phase or in the parallel per-window phase. Returns "" where unsupported. [mem]
  /// seconds since tool start, for per-phase bottleneck attribution [perf]
  static double phase_clock_()
  {
    static const auto t0 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  static std::string clk_() { char b[32]; snprintf(b, sizeof b, " [t=%.0fs]", phase_clock_()); return b; }

  static std::string rss_()
  {
    std::ifstream st("/proc/self/status");
    if (!st) return "";
    std::string line, cur = "?", peak = "?";
    while (std::getline(st, line))
    {
      auto val = [&line]() {
        std::string v = line.substr(line.find(':') + 1);
        size_t a = v.find_first_not_of(" \t");
        return a == std::string::npos ? std::string("?") : v.substr(a);
      };
      if (line.rfind("VmRSS:", 0) == 0) cur = val();
      else if (line.rfind("VmHWM:", 0) == 0) peak = val();
    }
    return " [RSS " + cur + ", peak " + peak + "]";
  }

  /// [dyn-mem] Bytes we may still allocate. MemAvailable is the kernel's own estimate of what
  /// can be handed out without swapping (it accounts for reclaimable page cache), which is the
  /// right quantity on a SHARED node -- MemFree would ignore cache and MemTotal would ignore
  /// the other users. Add our own RSS back: memory we already hold is part of our budget, not
  /// somebody else's. Returns 0 if unreadable, which the caller treats as "admit one at a time".
  static size_t availableBytes_()
  {
    std::ifstream mi("/proc/meminfo");
    if (!mi) return 0;
    std::string line;
    size_t avail_kb = 0;
    while (std::getline(mi, line))
      if (line.rfind("MemAvailable:", 0) == 0)
      { avail_kb = strtoull(line.c_str() + 13, nullptr, 10); break; }
    if (!avail_kb) return 0;
    std::ifstream st("/proc/self/status");
    size_t rss_kb = 0;
    while (st && std::getline(st, line))
      if (line.rfind("VmRSS:", 0) == 0)
      { rss_kb = strtoull(line.c_str() + 6, nullptr, 10); break; }
    return (avail_kb + rss_kb) * 1024ull;
  }

  struct Trace
  {
    double mz = 0.0;
    double rt = 0.0;   ///< apex/centroid RT
    double im = 0.0;   ///< centroid ion mobility (1/K0)
    double intensity = 0.0;
    // elution profile as (rt, intensity) points, sorted by rt
    vector<pair<double, double>> xic;
  };

  /// A precursor hypothesis after isotope/charge inference.
  struct Precursor_
  {
    double mono_mz = 0.0;
    int charge = 0;          ///< 0 == unknown
    double rt = 0.0;
    double im = 0.0;
    size_t trace_idx = 0;    ///< index into ms1 traces
    bool guessed = false;    ///< [way-4] charge was DEFAULTED, not isotope-supported
  };

  /// MassTraceDetection locates the per-peak IM array by the exact name
  /// Constants::UserParam::ION_MOBILITY ("Ion Mobility"). Converters name it via CV term and
  /// PeakPickerIM names it "Ion Mobility Centroid"; rename whichever exists so tracing is IM-aware.
  void ensureIMArrayName(MSSpectrum& s)
  {
    auto& fdas = s.getFloatDataArrays();
    if (s.containsIMData())
    {
      fdas[s.getIMData().first].setName(Constants::UserParam::ION_MOBILITY);
      return;
    }
    for (auto& fda : fdas)
    {
      if (fda.getName() == Constants::UserParam::ION_MOBILITY_CENTROID)
      {
        fda.setName(Constants::UserParam::ION_MOBILITY);
        return;
      }
    }
  }

  Trace toTrace(const MassTrace& mt)
  {
    Trace t;
    t.mz = mt.getCentroidMZ();
    t.rt = mt.getCentroidRT();
    t.im = mt.getCentroidIM();
    t.intensity = mt.getMaxIntensity(false);
    t.xic.reserve(mt.getSize());
    for (Size i = 0; i < mt.getSize(); ++i)
    {
      t.xic.emplace_back(mt[i].getRT(), mt[i].getIntensity());
    }
    sort(t.xic.begin(), t.xic.end());
    return t;
  }

  /// Per-window correlation grid: all fragment XICs are resampled onto ONE common RT axis
  /// (the sorted union of fragment-trace RTs) so a precursor can be correlated against every
  /// candidate fragment as a precomputed dot product instead of an O(P log F) per-pair search.
  /// This is the hot-path optimization: precursor XIC is prepared once per precursor (not per
  /// fragment), fragment means/norms are precomputed once per window. [perf]
  /// CSR layout: one `vector<vector<...>>` costs 24 B of header PER FRAGMENT before a single
  /// support point is stored -- ~200 MB/window at 8.5 M fragments, and scattered allocations.
  /// A flat array + offsets costs 4 B/fragment and is contiguous for the correlation loop. [compact]
  struct FragGrid
  {
    vector<double> rt;                          ///< sorted unique grid RTs (size G)
    vector<uint32_t> off;                       ///< size F+1; fragment i spans [off[i], off[i+1])
    vector<pair<int, float>> flat;              ///< (grid index, intensity), concatenated
    vector<double> mean, invnorm;               ///< per fragment: mean and 1/norm over the full grid
    size_t supBegin(size_t i) const { return off[i]; }
    size_t supEnd(size_t i)   const { return off[i + 1]; }
  };


  /// [cross-frame] Sliding-window aggregation of adjacent RT frames within ONE isolation window.
  ///
  /// diaTracer sums "adjacent neighbor RT frames" of the same isolation window into a composite
  /// (m/z x 1/K0) matrix BEFORE peak picking, which raises per-point S/N by ~sqrt(N) and lets weak
  /// fragments form traces at all (Nat Commun 16:95 Methods). Because a diaPASEF isolation window is
  /// sampled exactly ONCE per cycle, "adjacent frames" == adjacent CYCLES.
  ///
  /// SLIDING (stride 1), not block: block summing would divide our RT point count by N, and peaks here
  /// span only 4-22 points, which would leave too few points to correlate at all. The cost of sliding
  /// is that adjacent output points share input frames and are therefore AUTOCORRELATED, which inflates
  /// Pearson r for signal and noise alike - the reason this must be benchmarked, not assumed.
  ///
  /// Output spectrum i = intensity-weighted (m/z, IM) clustering of the peaks of spectra
  /// [i-N/2, i+N/2], carrying the ORIGINAL RT of spectrum i so the RT axis is unchanged.
  void aggregateFrames_(PeakMap& wmap, int n_frames, double ppm, double im_tol)
  {
    if (n_frames <= 1 || wmap.empty()) return;
    const int half = n_frames / 2;
    const Size N = wmap.size();

    struct Pk { double mz; double im; double in; };
    PeakMap out;
    out.reserve(N);

    vector<Pk> buf;
    for (Size i = 0; i < N; ++i)
    {
      const Size lo = (i > (Size)half) ? i - half : 0;
      const Size hi = std::min(N - 1, i + (Size)half);

      buf.clear();
      for (Size k = lo; k <= hi; ++k)
      {
        const MSSpectrum& sp = wmap[k];
        if (!sp.containsIMData()) continue;
        const auto& imarr = sp.getFloatDataArrays()[sp.getIMData().first];
        for (Size p = 0; p < sp.size(); ++p) buf.push_back({sp[p].getMZ(), (double)imarr[p], sp[p].getIntensity()});
      }
      sort(buf.begin(), buf.end(), [](const Pk& a, const Pk& b) {
        if (a.mz != b.mz) return a.mz < b.mz;
        return a.im < b.im;
      });

      MSSpectrum agg;
      agg.setRT(wmap[i].getRT());                 // RT axis preserved (sliding, stride 1)
      agg.setMSLevel(wmap[i].getMSLevel());
      agg.setPrecursors(wmap[i].getPrecursors());
      OpenMS::DataArrays::FloatDataArray im_out;
      im_out.setName(Constants::UserParam::ION_MOBILITY);   // keep tracing IM-aware

      // greedy cluster in m/z (ppm) then IM: same ion across adjacent cycles collapses, intensities SUM
      for (size_t a = 0; a < buf.size(); )
      {
        const double mz0 = buf[a].mz;
        const double tol = mz0 * ppm * 1e-6;
        size_t b = a;
        double sum = 0.0, wmz = 0.0, wim = 0.0;
        while (b < buf.size() && buf[b].mz - mz0 <= tol)
        {
          if (fabs(buf[b].im - buf[a].im) <= im_tol)
          {
            sum += buf[b].in; wmz += buf[b].mz * buf[b].in; wim += buf[b].im * buf[b].in;
            buf[b].in = -1.0;                    // consumed
          }
          ++b;
        }
        if (sum > 0.0)
        {
          Peak1D pk; pk.setMZ(wmz / sum); pk.setIntensity((float)sum);
          agg.push_back(pk);
          im_out.push_back((float)(wim / sum));
        }
        // advance to next unconsumed peak
        while (a < buf.size() && buf[a].in < 0.0) ++a;
        if (a >= buf.size()) break;
      }
      agg.getFloatDataArrays().push_back(im_out);
      agg.sortByPosition();
      out.addSpectrum(std::move(agg));
    }
    wmap = std::move(out);
  }

  /// [compact] Quantised frame store for the per-window buffer.
  ///
  /// Holding all 24 isolation windows as full OpenMS PeakMaps costs ~89 GB resident before any window
  /// is consumed, which is what forces `max_concurrent_windows` down and caps our parallelism at ~7
  /// cores (diaTracer uses ~43). A frame only needs (m/z, intensity, 1/K0) per peak plus one RT, so we
  /// store structure-of-arrays with quantised keys and materialise a PeakMap ONLY inside the worker
  /// that is about to trace it, freeing the compact copy immediately.
  ///
  /// Precision (deliberately far finer than every downstream tolerance, so tracing is unaffected):
  ///   m/z : uint32 at 1e-5 Da  -> 0.007 ppm at m/z 1400, vs a 20 ppm trace tolerance (2800x margin)
  ///   1/K0: uint16 over [0.4,1.8] -> 2.1e-5, vs a 0.01 IM tolerance (470x margin)
  ///   intensity: float32, unchanged (Peak1D already stores float)
  /// Cost: 10 B/peak (4+4+2) vs ~20 B/peak for Peak1D+FloatDataArray, plus no per-MSSpectrum overhead.
  constexpr double MZ_Q = 1e5;      ///< m/z quantum = 1e-5 Da
  constexpr double IM_LO = 0.4, IM_HI = 1.8;
  constexpr double IM_Q = 65535.0 / (IM_HI - IM_LO);

  struct CompactFrame
  {
    double rt = 0.0;
    vector<uint32_t> mzq;
    vector<float> inten;
    vector<uint16_t> imq;
    size_t bytes() const { return mzq.size() * 4 + inten.size() * 4 + imq.size() * 2; }
  };

  /// Peaks that cannot be represented exactly are DROPPED and COUNTED, never silently altered.
  /// [adv-fix] The first version CLAMPED out-of-range 1/K0 to IM_LO/IM_HI, which would pin outliers
  /// from many frames onto two artificial, perfectly-aligned mobility values — manufacturing spurious
  /// traces out of nothing. NaN also passed straight through clamp into llround (UB on unsigned
  /// conversion). Dropping + counting is the only defensible behaviour.
  struct CompactStats { size_t no_im_array = 0, size_mismatch = 0, bad_mz = 0, bad_im = 0, kept = 0; };

  /// Convert one picked MSSpectrum into the compact form (drops everything tracing does not read).
  CompactFrame compactify(const MSSpectrum& s, CompactStats& st)
  {
    CompactFrame f;
    f.rt = s.getRT();
    if (!s.containsIMData()) { st.no_im_array += s.size(); return f; }
    const auto& im = s.getFloatDataArrays()[s.getIMData().first];
    const Size n = s.size();
    if (im.size() != n) { st.size_mismatch += n; return f; }   // [adv-fix] no out-of-bounds read
    f.mzq.reserve(n); f.inten.reserve(n); f.imq.reserve(n);
    for (Size i = 0; i < n; ++i)
    {
      const double mz = s[i].getMZ();
      const double imv = (double)im[i];
      // [adv-fix] validate BEFORE the unsigned conversion (NaN/inf/negative -> UB otherwise)
      if (!std::isfinite(mz) || mz <= 0.0 || mz * MZ_Q > 4.2e9) { ++st.bad_mz; continue; }
      if (!std::isfinite(imv) || imv < IM_LO || imv > IM_HI)    { ++st.bad_im; continue; }
      f.mzq.push_back((uint32_t)llround(mz * MZ_Q));
      f.inten.push_back(s[i].getIntensity());
      f.imq.push_back((uint16_t)llround((imv - IM_LO) * IM_Q));
      ++st.kept;
    }
    return f;
  }

  /// Rebuild a PeakMap for ONE window from its compact frames, for MassTraceDetection.
  /// Frames are consumed (cleared) as they are materialised so the compact store shrinks in step.
  PeakMap materializeWindow(vector<CompactFrame>& frames)
  {
    PeakMap m;
    m.reserve(frames.size());
    for (auto& f : frames)
    {
      MSSpectrum s;
      s.setRT(f.rt);
      s.setMSLevel(1); // MassTraceDetection traces MS1-level only
      s.reserve(f.mzq.size());
      OpenMS::DataArrays::FloatDataArray ima;
      ima.setName(Constants::UserParam::ION_MOBILITY);
      ima.reserve(f.mzq.size());
      for (size_t i = 0; i < f.mzq.size(); ++i)
      {
        Peak1D p;
        p.setMZ((double)f.mzq[i] / MZ_Q);
        p.setIntensity(f.inten[i]);
        s.push_back(p);
        ima.push_back((float)(IM_LO + (double)f.imq[i] / IM_Q));
      }
      s.getFloatDataArrays().push_back(std::move(ima));
      m.addSpectrum(std::move(s));
      // release this frame's compact storage immediately
      vector<uint32_t>().swap(f.mzq);
      vector<float>().swap(f.inten);
      vector<uint16_t>().swap(f.imq);
    }
    vector<CompactFrame>().swap(frames);
    return m;
  }

  /// [stream] Pick-and-compact consumer: removes the one-shot-read memory floor.
  ///
  /// MEASURED: `loadExperiment` reaches 90.3 GB RSS / 99.0 GB peak at t=0s, BEFORE any compaction
  /// exists — so the compact store (10 B/peak, 20.9 GB for 2.19e9 peaks) can never lower the peak.
  /// The floor is the whole run being resident at once. BrukerTimsFile::loadDIAStreaming hands us one
  /// frame at a time, so we peak-pick and compact each frame on arrival and never hold the raw run.
  ///
  /// Windows are keyed by the SAME (lo,hi) isolation-window key used downstream, taken from the
  /// spectrum's own precursor, so routing is identical to the old split - not by swath_nr, whose
  /// ordering is the reader's business and need not match ours.
  class PickCompactConsumer : public FullSwathFileConsumer
  {
  public:
    PickCompactConsumer(PeakPickerIM& picker, map<pair<int, int>, vector<CompactFrame>>& windows,
                        PeakMap& ms1, CompactStats& stats,
                        std::function<pair<int, int>(double, double)> keyfn)
      : picker_(picker), windows_(windows), ms1_(ms1), stats_(stats), keyfn_(std::move(keyfn)) {}

    size_t frames_seen = 0;

  protected:
    void consumeMS1Spectrum_(MapType::SpectrumType& s) override
    {
      // MS1 stays a PeakMap: it is only ~4% of frames (1,343 of 33,553) and MassTraceDetection
      // needs it whole for the run-level precursor pass.
      ensureIMArrayName(s);
      if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) picker_.pickIMCluster(s);
      ensureIMArrayName(s);
      ms1_.addSpectrum(std::move(s));
      ++frames_seen;
    }

    void consumeSwathSpectrum_(MapType::SpectrumType& s, size_t /*swath_nr*/) override
    {
      const auto& prec = s.getPrecursors();
      if (prec.empty()) return;                      // no isolation window -> not routable
      const double c = prec[0].getMZ();
      const double lo = c - prec[0].getIsolationWindowLowerOffset();
      const double hi = c + prec[0].getIsolationWindowUpperOffset();
      ensureIMArrayName(s);
      if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) picker_.pickIMCluster(s);
      ensureIMArrayName(s);
      windows_[keyfn_(lo, hi)].push_back(compactify(s, stats_));   // raw frame dies here
      ++frames_seen;
    }

    void ensureMapsAreFilled_() override {}          // we own the containers; nothing to finalise

  private:
    PeakPickerIM& picker_;
    map<pair<int, int>, vector<CompactFrame>>& windows_;
    PeakMap& ms1_;
    CompactStats& stats_;
    std::function<pair<int, int>(double, double)> keyfn_;
  };

  /// Build the per-window FragGrid from IM-sorted fragment traces.
  FragGrid buildFragGrid(const vector<Trace>& frags)
  {
    FragGrid g;
    // union of all fragment RTs -> the common grid
    for (const auto& f : frags) for (const auto& p : f.xic) g.rt.push_back(p.first);
    sort(g.rt.begin(), g.rt.end());
    g.rt.erase(unique(g.rt.begin(), g.rt.end()), g.rt.end());
    const double G = (double)g.rt.size();
    g.mean.resize(frags.size());
    g.invnorm.resize(frags.size());
    g.off.resize(frags.size() + 1, 0);
    size_t total = 0;
    for (const auto& f : frags) total += f.xic.size();
    g.flat.reserve(total);
    for (size_t i = 0; i < frags.size(); ++i)
    {
      double sum = 0, sumsq = 0;
      g.off[i] = (uint32_t)g.flat.size();
      for (const auto& p : frags[i].xic)
      {
        int gi = (int)(lower_bound(g.rt.begin(), g.rt.end(), p.first) - g.rt.begin());
        g.flat.emplace_back(gi, (float)p.second);
        sum += p.second; sumsq += p.second * p.second;
      }
      g.mean[i] = sum / G;
      double var = sumsq - G * g.mean[i] * g.mean[i];
      g.invnorm[i] = var > 0 ? 1.0 / sqrt(var) : 0.0; // 0 => constant/degenerate -> never correlates
    }
    g.off[frags.size()] = (uint32_t)g.flat.size();   // CSR sentinel
    return g;
  }
}

class TOPPDIAPseudoSpectra : public TOPPBase
{
public:
  TOPPDIAPseudoSpectra() :
    TOPPBase("DIAPseudoSpectra", "Extracts pseudo-MS/MS spectra from diaPASEF (ion-mobility DIA) data.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input diaPASEF data (ion-mobility DIA; 1/K0 / VSSC): mzML, mzPeak, or a Bruker .d directory.");
    setValidFormats_("in", {"mzML", "mzpeak", "d"});
    registerOutputFile_("out", "<file>", "", "Output pseudo-MS2 mzML (searchable by a DDA engine).");
    setValidFormats_("out", {"mzML"});

    registerTOPPSubsection_("gate", "Precursor-to-fragment co-localization gate");
    // Conservative, documented defaults; these are placeholders pending a real-data sweep. [H-3]
    registerDoubleOption_("gate:delta_im", "<1/K0>", 0.01, "Max |precursor - fragment| ion mobility difference (1/K0). (diaTracer deltaApexIM default 0.01)", false);
    registerDoubleOption_("gate:delta_rt", "<sec>", 3.0, "Max |precursor - fragment| apex RT difference (seconds).", false);
    registerDoubleOption_("gate:min_correlation", "<0..1>", 0.3, "Min Pearson correlation of elution profiles. (diaTracer ms1MS2Corr default 0.3; raise for cleaner but sparser spectra)", false);
    registerIntOption_("gate:min_correlation_points", "<n>", 3, "Min overlapping XIC points required to accept a correlation.", false, true);

    registerTOPPSubsection_("assembly", "Pseudo-spectrum assembly");
    registerIntOption_("assembly:min_fragments", "<n>", 3, "Emit a pseudo-spectrum only if it has at least this many fragments.", false);
    registerIntOption_("assembly:max_fragments", "<n>", 500, "Keep at most this many (top-ranked) fragments per pseudo-spectrum.", false);
    registerFlag_("assembly:competitive", "Assign each fragment to only its single best-correlating precursor (de-chimerize). Shorthand for rp_max=1; the harder endpoint of the rp_max dial. [bench]");
    registerDoubleOption_("assembly:apportion", "<p>", 0.0, "[route-4] Split a shared fragment's INTENSITY across the precursors that claim it, weight w_i = corr_i^p / sum(corr^p), instead of copying it at full intensity into all of them (mean fan-out 6.45). Fragment COUNT is unchanged - only the weighting - so it sidesteps every falsified count knob (rank-pruning, competitive, min_corr, min_corr_pts). 0 = OFF. Try 1 (linear) or 2 (sharper); p->inf degenerates to competitive, which IS falsified. Share-all path only.", false);
    registerIntOption_("assembly:rp_max", "<n>", 0, "Soft RP rank-pruning (DIA-Umpire style): keep each fragment in only its top-N best-correlating precursors. 0 = unlimited (share-all, the default and current best); 1 = competitive (winner-take-all). Small N (2-8) is the useful range given typical fan-out ~6. RFmax (per-precursor fragment cap) is 'max_fragments'. Runs emit a fan-out histogram so you can see whether a given rp_max actually bites. [bench]", false);
    registerFlag_("assembly:open_search_safe", "[way-4] ANNOTATION-SAFE output for BLIND/OPEN search. In open search a wrong charge is catastrophic in a way closed search hides: assigning z' instead of z shifts the neutral mass by dM=(z'-z)(m/z-1.0073), so ONE charge step at m/z 500 fabricates a ~499 Da 'modification'. Closed search barely notices (charge unset cost only 1.3%: 8,123->8,019) because the engine enumerates charge; an open engine trusts the annotation. This flag therefore: (a) emits charge UNSET for every precursor whose charge was GUESSED rather than isotope-supported, so the engine enumerates instead of trusting a default; (b) records the isotope-offset hypothesis as a meta value so delta-mass inference can correct BEFORE assigning a modification, rather than reporting phantom +-1.00335/+-2.0067 Da shifts. Use for open/variant search; leave off for closed-search peptide counting.");
    registerIntOption_("assembly:default_charge", "<z>", 2, "Charge for precursors with no isotope support: >0 assigns that charge (diaTracer-style; 2 = most common); 0 = leave charge UNSET so the search engine searches its charge range and picks the right one (recovers 3+/4+ that a forced 2 would lose).", false);

    registerTOPPSubsection_("trace", "Mass-trace detection (see MassTraceDetection)");
    registerDoubleOption_("trace:mass_error_ppm", "<ppm>", 15.0, "m/z tolerance for trace detection.", false);
    registerDoubleOption_("trace:noise_threshold_int", "<int>", 100.0, "MS1 (precursor) noise intensity threshold.", false);
    registerDoubleOption_("trace:ms2_noise_threshold_int", "<int>", 100.0, "MS2 (fragment) noise threshold. Lower it (e.g. 50, 20) to recover weak fragments toward diaTracer-like density; the grid-Pearson optimization keeps that affordable.", false);
    registerDoubleOption_("trace:min_length_sec", "<sec>", 3.0, "Minimum MS1 mass-trace length (seconds).", false, true);
    // diaTracer runs MS2 with deliberately RELAXED constraints vs MS1 (Nat Commun 16:95). The OpenMS
    // apex gate is chrom_peak_snr*noise_threshold_int, so the stock snr=3.0 silently TRIPLES the
    // fragment threshold (a "noise 30" run really gated apices at 90). Decouple both for MS2. [recall]
    registerDoubleOption_("trace:ms1_chrom_peak_snr", "<x>", 3.0, "Apex signal-to-noise multiplier for MS1 traces (effective apex threshold = this x noise_threshold_int). OpenMS default 3.0.", false, true);
    registerDoubleOption_("trace:ms2_chrom_peak_snr", "<x>", 1.0, "Apex signal-to-noise multiplier for MS2 (fragment) traces. 1.0 = the ms2_noise_threshold_int you set IS the apex threshold; 3.0 would reproduce the old hidden 3x. Lower = more weak fragments recovered.", false);
    registerDoubleOption_("trace:ms2_min_length_sec", "<sec>", 0.0, "Minimum MS2 (fragment) mass-trace length (seconds). 0 = no length filter (diaTracer-style relaxed MS2); MS1 keeps trace:min_length_sec.", false);
    registerDoubleOption_("trace:ms2_split_valleys", "<chrom_fwhm>", 0.0, "[way-2] Split MS2 mass traces at chromatographic local minima via ElutionPeakDetection (which MassTraceDetection never does - it only terminates on outliers, so two peptides 25 s apart within the m/z tolerance MERGE). Value is chrom_fwhm in seconds; it sets both the Savitzky-Golay window and the local-extrema half-window. 0 = OFF. Try 6-8 (peaks here are 5-30 s). width_filtering is forced off so this only SPLITS, never deletes.", false);
    registerDoubleOption_("trace:ms1_split_valleys", "<chrom_fwhm>", 0.0, "[way-2] Same, for MS1 traces. Merged MS1 traces produce merged precursors, so this is the side that can change WHICH precursors exist. 0 = OFF.", false);
    registerDoubleOption_("trace:max_trace_length_sec", "<sec>", -1.0, "DANGER - NOT a guard. MassTraceDetection::isTraceValid_ RETURNS FALSE on over-length traces and peaks are marked visited only for VALID traces, so an over-length blob is rejected, re-seeded by the next apex, and rejected again - emitting NOTHING at that m/z. Setting this DELETES signal (the falsified direction). Kept only for diagnostics; leave at -1. Splitting requires ElutionPeakDetection, not this. [merged-trace]", false, true);
registerIntOption_("trace:frame_aggregation_ms1_n", "<n>", 1, "[cross-frame] Same as frame_aggregation_n but for the MS1 (PRECURSOR) map. This is the one that can add DEPTH: measurement shows we identify MORE spectra than diaTracer (3.7%% vs 2.0%% of spectra) but cover FEWER distinct peptides (4.6 vs 1.46 PSMs/peptide), i.e. we miss low-abundance PRECURSORS. Aggregating MS2 alone cannot fix that. 1 = off.", false);
    registerIntOption_("trace:frame_aggregation_n", "<n>", 1, "[cross-frame] Number of ADJACENT RT frames (= adjacent diaPASEF CYCLES, since each isolation window is sampled once per cycle) summed into a composite (m/z x 1/K0) frame BEFORE MS2 peak/trace detection, diaTracer-style. 1 = off (current behaviour). 5 = user request. Sliding window, stride 1, so the RT axis and point count are preserved; note adjacent output points then share input frames and are AUTOCORRELATED, which inflates Pearson r - benchmark, do not assume.", false);
    registerDoubleOption_("trace:ms2_min_sample_rate", "<r>", -1.0, "MS2 min fraction of visited scans that must contain a peak (MassTraceDetection min_sample_rate). OpenMS default 0.5 hard-deletes gappy fragment traces; MEASURED: 0.3 and 0.1 both LOSE peptides (8,865 -> 8,735 at 0.3), so -1 (= OpenMS 0.5) is the default. Lowering it produces longer, more-merged traces with wrong apexes.", false);
    registerDoubleOption_("trace:ms1_min_sample_rate", "<r>", -1.0, "MS1 min_sample_rate. -1 = OpenMS default 0.5 (strict, MS1 stays strict).", false, true);

    // NB: "threads" is reserved by TOPPBase (the standard -threads option), so this lives under "perf".
    // [route-1] feature-level consolidation: collapse spectra that are the SAME chromatographic
    // feature emitted more than once (RT slices / IM sub-ranges / charge hypotheses). Targets the
    // measured 4.6 vs 1.46 PSMs-per-peptide gap. OFF by default (delta_rt 0) - it is a deletion.
    registerTOPPSubsection_("consolidate", "Feature-level de-duplication of emitted spectra");
    registerDoubleOption_("consolidate:delta_rt", "<sec>", 0.0, "Consolidate spectra whose precursors agree within this RT window (and mass_ppm/delta_im) into ONE feature, keeping the richest spectrum. 0 = OFF. Try 5-10. Attacks emission multiplicity (ours 4.6 PSMs/peptide vs diaTracer 1.46) without changing fragment content.", false);
    registerDoubleOption_("consolidate:mass_ppm", "<ppm>", 20.0, "Precursor m/z tolerance for feature consolidation.", false, true);
    registerDoubleOption_("consolidate:delta_im", "<1/K0>", 0.02, "Precursor IM tolerance for feature consolidation. Set larger than the 12 m/z x 2 IM sub-range split to merge a feature emitted from BOTH IM sub-ranges.", false, true);
    registerFlag_("consolidate:same_charge_only", "Only merge spectra that also agree in charge (safer; leaves genuine multi-charge observations intact).");

    registerFlag_("perf:stream_load", "[stream] Read the .d frame-by-frame via BrukerTimsFile::loadDIAStreaming, peak-picking and compacting each frame on arrival, instead of loadExperiment() which holds the entire run (MEASURED 90.3 GB RSS at t=0s - that one-shot read IS the memory floor and no downstream compaction can lower it). Bruker .d only.");
    registerIntOption_("trace:native_ms2_neighbors", "<n>", 0, "[stream] NATIVE cross-frame aggregation in the READER (BrukerTimsFile Config dia_ms2_n_neighbors): adjacent frames each side, 0=off, 1=3-frame sum, 2=5-frame sum. This sums RAW frames BEFORE centroiding, which is what diaTracer does; the tool-level frame_aggregation_n sums AFTER centroiding and lost 26% of peptides. Requires perf:stream_load.", false);
    registerIntOption_("trace:native_ms1_neighbors", "<n>", 0, "[stream] Same, for MS1 frames (BrukerTimsFile Config ms1_n_neighbors). Requires perf:stream_load.", false);

    registerTOPPSubsection_("perf", "Concurrency / memory tradeoff");
    registerIntOption_("perf:trace_bands", "<n>", 0, "[parallel] Split each window's fragment m/z range into N bands traced CONCURRENTLY. The diaPASEF method fixes the window count (24 here), which caps window-level parallelism far below the core count; banding lifts the ceiling to threads. Traces cannot span more than trace:mass_error_ppm in m/z, so banding with a halo is exact. 0 = auto (windows x bands >= threads), 1 = off.", false);
    registerDoubleOption_("perf:mem_fraction", "<f>", 0.75, "[dyn-mem] Fraction of currently-FREE RAM (/proc/meminfo MemAvailable + our own RSS) the window loop may commit. Concurrency is re-decided at every window admission, so the run adapts to other jobs on a shared node. One window is always admitted regardless, so a single oversized window cannot deadlock.", false);
    registerIntOption_("perf:max_concurrent_windows", "<n>", 0, "UPPER bound on isolation windows processed concurrently (memory admission may run fewer). Peak RAM scales with windows IN FLIGHT (each holds its frames + traces + grid), not with the total window count, so lowering this trades wall time for RAM. 0 = unlimited (use all threads).", false);

    registerTOPPSubsection_("diag", "Diagnostics ('debug' is reserved by TOPP)");
    registerStringOption_("diag:dump_ms1_tsv", "<prefix>", "", "[ms1-funnel] Write MS1 traces and inferred precursors to <prefix>.traces.tsv / <prefix>.precursors.tsv so precursor loss can be attributed to a stage (no trace / no hypothesis / no spectrum) rather than inferred from the final count. Empty = off.", false, true);

    registerIntOption_("max_charge", "<n>", 5, "Maximum precursor charge considered during isotope inference.", false);
    // NOTE: the subsection MUST be registered or printUsage_() throws ElementNotFound — which turns
    // any user parameter typo into a FATAL crash instead of a usage message. [bugfix]
    registerTOPPSubsection_("charge", "Isotope / charge-state inference");
    // [route-1] Only 25.7%% of diaTracer's precursors are matched by us with the RIGHT m/z AND charge;
    // 17.2%% have the wrong charge and 22.9%% (decoy-corrected) the wrong monoisotope. A wrong charge or
    // mono gives a wrong NEUTRAL MASS, so the peptide can never be identified no matter how good the
    // fragments are. Partner-COUNTING cannot tell a real envelope from a coincidence; averagine shape
    // x isotope co-elution can.
    registerStringOption_("charge:scoring", "<mode>", "envelope", "Charge/monoisotope inference: 'envelope' scores each (charge, mono) hypothesis by averagine-shape cosine x isotope-XIC co-elution and allows GAPPED envelopes; 'count' is the legacy partner-count heuristic (ties favoured LOW charge, break-on-first-miss).", false);
    setValidStrings_("charge:scoring", {"envelope", "count"});
    registerDoubleOption_("charge:ambiguity_margin", "<x>", 0.0, "[route-3] If a different-charge hypothesis scores within this margin of the best, emit it too (diaTracer retains multiple charges when ambiguous). 0 = single charge per precursor.", false);
    registerFlag_("assembly:dedup_precursors", "[route-4] Collapse precursor hypotheses that are the same species (same charge, m/z within ppm, co-eluting RT, co-located IM). Measured redundancy: 4.6 identified spectra per unique peptide vs diaTracer's 1.46.");
    registerDoubleOption_("charge:iso_im_tolerance", "<1/K0>", 0.05, "IM tolerance for matching isotope partners. [adv-fix] Was 0.05 = 5x the fragment gate, which is indefensible: isotopes of ONE species share essentially identical mobility, so 0.05 admits cross-species partners that corrupt both the envelope cosine and the co-elution score. 0.01 matches the fragment gate.", false);
    // ponytail: mass-defect filter default OFF (would silently drop PTM/nonspecific — the discovery cases). [H-7]
  }

  /// Configure an IM-aware MassTraceDetection and run it on @p map.
  /// MassTraceDetection wrapper. NOTE the apex gate inside MassTraceDetection is
  ///   peak_int > chrom_peak_snr * noise_threshold_int
  /// so the EFFECTIVE apex threshold is snr*noise, while trace MEMBERSHIP only needs > noise.
  /// Leaving snr at the OpenMS default 3.0 silently triples the fragment apex threshold.
  /// diaTracer deliberately runs MS2 with "more relaxed constraints" than MS1 (Nat Commun 16:95
  /// Methods) — so MS1 and MS2 get INDEPENDENT snr / min-length here. [recall]
  /// `map` is taken by NON-const ref and CLEARED as soon as MassTraceDetection has consumed it:
  /// it is dead from that point on, but at ~20 B/peak (Peak1D 16 + IM float 4) it is ~1.1 GB per
  /// window that would otherwise stay resident through the valley-splitting peak. [compact]
  vector<Trace> detectTraces_(PeakMap& map, double delta_im, double noise,
                              double snr, double min_len, double min_sample_rate, double max_len,
                              double split_valleys, String* span_log = nullptr, int bands = 1)
  {
    MassTraceDetection mtd;
    mtd.setLogType(ProgressLogger::NONE); // quiet + thread-safe (called from the parallel window loop)
    Param p = mtd.getParameters();
    p.setValue("mass_error_ppm", getDoubleOption_("trace:mass_error_ppm"));
    p.setValue("noise_threshold_int", noise);
    p.setValue("chrom_peak_snr", snr);
    p.setValue("ion_mobility_tolerance", delta_im);        // IM-aware tracing [C-1]
    p.setValue("min_trace_length", min_len);
    // min_sample_rate (OpenMS default 0.5) requires a peak in >=50% of visited scans AND terminates
    // trace extension when occupancy drops below it — an intensity-INDEPENDENT deletion of gappy but
    // real fragment traces in sparse diaPASEF data. Relaxing it for MS2 is the sole surviving MS2
    // trace-validity lever. [recall]
    if (min_sample_rate >= 0.0) p.setValue("min_sample_rate", min_sample_rate);
    if (max_len > 0.0) p.setValue("max_trace_length", max_len);
    mtd.setParameters(p);

    vector<MassTrace> mts;
    if (bands > 1)
    {
      // [bands] The diaPASEF method fixes the window count at 24, so parallelising over
      // windows caps us at 24 threads however many cores exist. Fragment traces at
      // different m/z are INDEPENDENT, so a window's peaks can be partitioned by m/z and
      // traced concurrently. A trace can never span more than mass_error_ppm in m/z, so a
      // halo of 20x that tolerance makes the partition exact; each trace is kept only by the
      // band whose CORE contains its centroid, so halo duplicates are dropped rather than merged.
      const double ppm = getDoubleOption_("trace:mass_error_ppm");
      // Band edges from m/z QUANTILES, not uniform m/z: peak density varies by orders of
      // magnitude across a window, and uniform edges would leave one band holding most of the work.
      vector<double> sample;
      for (const auto& s : map)
        for (Size i = 0; i < s.size(); i += 97) sample.push_back(s[i].getMZ());
      sort(sample.begin(), sample.end());
      vector<double> edge(bands + 1);
      if (sample.empty()) { edge.assign(bands + 1, 0.0); }
      else
        for (int b = 0; b <= bands; ++b)
          edge[b] = sample[std::min(sample.size() - 1, (size_t)((double)b / bands * sample.size()))];
      edge.front() = 0.0; edge.back() = std::numeric_limits<double>::max();
      vector<double> ().swap(sample);

      // Partition ONCE (peaks are copied into exactly one core band plus any halo it falls in),
      // then release the source map so the halo duplicates are the only overhead.
      vector<PeakMap> sub(bands);
      for (const auto& s : map)
      {
        const auto* ima = s.getFloatDataArrays().empty() ? nullptr : &s.getFloatDataArrays()[0];
        for (int b = 0; b < bands; ++b)
        {
          const double h = edge[b] * ppm * 1e-6 * 20.0;
          const double lo = edge[b] - h, hi = edge[b + 1] + h;
          MSSpectrum t; t.setRT(s.getRT()); t.setMSLevel(1);
          OpenMS::DataArrays::FloatDataArray ia;
          ia.setName(Constants::UserParam::ION_MOBILITY);
          for (Size i = 0; i < s.size(); ++i)
            if (s[i].getMZ() >= lo && s[i].getMZ() < hi)
            { t.push_back(s[i]); if (ima && i < ima->size()) ia.push_back((*ima)[i]); }
          if (!t.empty()) { t.getFloatDataArrays().push_back(std::move(ia)); sub[b].addSpectrum(std::move(t)); }
        }
      }
      map.clear(true);
      vector<vector<MassTrace>> per(bands);
      #pragma omp parallel for schedule(dynamic) num_threads(bands)
      for (int b = 0; b < bands; ++b)
      {
        MassTraceDetection m2; m2.setLogType(ProgressLogger::NONE); m2.setParameters(p);
        vector<MassTrace> t;
        sub[b].sortSpectra();
        m2.run(sub[b], t);
        sub[b].clear(true);
        // keep only traces whose centroid lies in this band's CORE -> no duplicates from the halo
        for (auto& mt : t)
        {
          const double c = mt.getCentroidMZ();
          if (c >= edge[b] && c < edge[b + 1]) per[b].push_back(std::move(mt));
        }
      }
      size_t tot = 0; for (const auto& v : per) tot += v.size();
      mts.reserve(tot);
      for (auto& v : per) { for (auto& t : v) mts.push_back(std::move(t)); vector<MassTrace>().swap(v); }
    }
    else
    {
      mtd.run(map, mts);
      map.clear(true); // [compact] dead from here; ~1.1 GB/window freed before the EPD peak
    }

    // [way-2] SPLIT MASS TRACES AT CHROMATOGRAPHIC VALLEYS.
    // MassTraceDetection never splits at a local minimum: it terminates only on consecutive-miss
    // outliers or sample_rate, and max_trace_length only DISCARDS. So two peptides eluting 25 s apart
    // within the m/z tolerance merge into ONE multi-modal trace, whose XIC then correlates with
    // everything under it. ElutionPeakDetection is the OpenMS stage that actually smooths and splits
    // at local minima - we never ran it. width_filtering is forced to "off" because "fixed" is a
    // DELETION filter ([min_fwhm,max_fwhm]) and blanket deletion is on the falsified list.
    if (split_valleys > 0.0 && !mts.empty())
    {
      ElutionPeakDetection epd;
      Param ep = epd.getParameters();
      ep.setValue("chrom_fwhm", split_valleys);      // sets SG window AND the local-extrema half-window
      ep.setValue("width_filtering", "off");         // never delete on width - only split
      ep.setValue("masstrace_snr_filtering", "false");
      epd.setParameters(ep);
      // MEMORY: detectPeaks(mts, split) materialises the ENTIRE output while the entire input
      // is still live -- measured as ~2.9 GB of a ~5.8 GB per-window peak, the single largest
      // term. Feed it in batches instead and free each batch as it is consumed, so the
      // transient is one batch rather than a full second copy. [compact]
      const size_t BATCH = 250000;
      const size_t n_in = mts.size();
      vector<MassTrace> out;
      out.reserve(n_in);                              // split only ever grows the count
      for (size_t b = 0; b < n_in; b += BATCH)
      {
        const size_t e = std::min(n_in, b + BATCH);
        vector<MassTrace> chunk(std::make_move_iterator(mts.begin() + b),
                                std::make_move_iterator(mts.begin() + e));
        // release the consumed slice's heap NOW, not at end of loop (MassTrace has no swap())
        for (size_t i = b; i < e; ++i) mts[i] = MassTrace();
        vector<MassTrace> part;
        epd.detectPeaks(chunk, part);
        vector<MassTrace>().swap(chunk);
        for (auto& t : part) out.push_back(std::move(t));
      }
      if (!out.empty())
      {
        if (span_log) *span_log += " | valley-split " + String(n_in) + " -> " + String(out.size())
                                 + " traces (" + String((int)(1000.0 * out.size() / n_in) / 10.0) + "%)";
        mts.swap(out);
      }
    }
    // [merged-trace] instrumentation: is trace merging actually happening? Log the RT-span
    // distribution. Peaks are 5-30 s here, so a large p90/max means multi-modal merged traces.
    if (!mts.empty())
    {
      vector<double> spans;
      spans.reserve(mts.size());
      for (const auto& mt : mts)
      {
        double lo = 1e30, hi = -1e30;
        for (Size i = 0; i < mt.getSize(); ++i) { double r = mt[i].getRT(); lo = min(lo, r); hi = max(hi, r); }
        spans.push_back(hi - lo);
      }
      sort(spans.begin(), spans.end());
      size_t n = spans.size();
      size_t over45 = spans.end() - lower_bound(spans.begin(), spans.end(), 45.0);
      if (span_log) *span_log = "traces=" + String(n) + " span_s median=" + String(spans[n/2])
                      + " p90=" + String(spans[(size_t)(0.9*n)]) + " max=" + String(spans.back())
                      + " frac>45s=" + String((int)(1000.0*over45/n)/10.0) + "%";
    }
    vector<Trace> out;
    out.reserve(mts.size());
    for (const auto& mt : mts) out.push_back(toTrace(mt));
    return out;
  }

  /// [route-1] Averagine isotope envelope, Poisson approximation: for neutral mass M the expected
  /// number of 13C is lambda ~ 0.000594*M, giving relative intensities [1, l, l^2/2, l^3/6, ...].
  /// This is what lets us pick the MONOISOTOPE: the mono is NOT always the most intense peak
  /// (above ~1800 Da M+1 exceeds M), so intensity-ranking alone mis-assigns it.
  static void averagineRatios(double neutral_mass, int n, vector<double>& out)
  {
    const double lam = std::max(0.10, 0.000594 * neutral_mass);
    out.assign(n, 0.0);
    double term = 1.0;
    for (int k = 0; k < n; ++k) { out[k] = term; term *= lam / (double)(k + 1); }
    double mx = 0.0; for (double v : out) mx = std::max(mx, v);
    if (mx > 0) for (double& v : out) v /= mx;
  }

  static double cosineSim(const vector<double>& a, const vector<double>& b)
  {
    double dot = 0, na = 0, nb = 0;
    for (size_t i = 0; i < a.size() && i < b.size(); ++i) { dot += a[i]*b[i]; na += a[i]*a[i]; nb += b[i]*b[i]; }
    return (na > 0 && nb > 0) ? dot / (sqrt(na) * sqrt(nb)) : 0.0;
  }

  /// [route-1] Pearson of two XICs over their OVERLAPPING RT points, properly mean-centred on the
  /// overlap (NOT zero-padded over a global grid). Isotopes of one peptide must co-elute; a
  /// coincidental peak at the right m/z spacing usually does not.
  static double xicCorr(const vector<pair<double,double>>& a, const vector<pair<double,double>>& b, double rt_tol)
  {
    vector<double> x, y;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size())
    {
      double d = a[i].first - b[j].first;
      if (fabs(d) <= rt_tol) { x.push_back(a[i].second); y.push_back(b[j].second); ++i; ++j; }
      else if (d < 0) ++i; else ++j;
    }
    const size_t n = x.size();
    if (n < 5) return 0.0;  // [adv-fix] n=3 gives P(r>0.9)~14.4% under the null
    double mx = 0, my = 0;
    for (size_t k = 0; k < n; ++k) { mx += x[k]; my += y[k]; }
    mx /= n; my /= n;
    double num = 0, sx = 0, sy = 0;
    for (size_t k = 0; k < n; ++k) { double dx = x[k]-mx, dy = y[k]-my; num += dx*dy; sx += dx*dx; sy += dy*dy; }
    return (sx > 0 && sy > 0) ? num / sqrt(sx*sy) : 0.0;
  }

  /// [route-1] Isotope/charge inference on MS1 traces. Scores each (charge, monoisotope) hypothesis
  /// by averagine-shape agreement x isotope-XIC co-elution, instead of counting partners.
  vector<Precursor_> inferPrecursors_(const vector<Trace>& ms1, int max_charge,
                                      double delta_rt, double iso_im_tol, double mass_ppm,
                                      bool envelope_scoring, double ambig_margin)
  {
    const double ISO = 1.0033548;
    const double PROTON_ = 1.00727646;
    const size_t N = ms1.size();

    // m/z-sorted index for O(log N) isotope-partner lookup (avoids O(N^2) over all MS1 traces). [H-8]
    vector<size_t> by_mz(N);
    for (size_t i = 0; i < N; ++i) by_mz[i] = i;
    sort(by_mz.begin(), by_mz.end(), [&](size_t a, size_t b) { return ms1[a].mz < ms1[b].mz; });
    vector<double> sorted_mz(N);
    for (size_t i = 0; i < N; ++i) sorted_mz[i] = ms1[by_mz[i]].mz;

    vector<size_t> order(N);
    for (size_t i = 0; i < N; ++i) order[i] = i;
    sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      if (ms1[a].intensity != ms1[b].intensity) return ms1[a].intensity > ms1[b].intensity;
      if (ms1[a].mz != ms1[b].mz) return ms1[a].mz < ms1[b].mz;
      return a < b; // total order so greedy isotope assignment is deterministic [code-review]
    });

    vector<bool> used(N, false);
    vector<Precursor_> out;
    for (size_t seed : order)
    {
      if (used[seed]) continue;
      const Trace& s = ms1[seed];
      const double tol = s.mz * mass_ppm * 1e-6;

      // Find an unused isotope partner near `target` m/z, co-localized in RT and IM (binary-search
      // the m/z index, scan only the ppm-window candidates); -1 if none. [H-8]
      auto findPartner = [&](double target) -> long {
        size_t klo = lower_bound(sorted_mz.begin(), sorted_mz.end(), target - tol) - sorted_mz.begin();
        size_t khi = upper_bound(sorted_mz.begin(), sorted_mz.end(), target + tol) - sorted_mz.begin();
        for (size_t k = klo; k < khi; ++k)
        {
          size_t j = by_mz[k];
          if (used[j] || j == seed) continue;
          if (fabs(ms1[j].rt - s.rt) <= delta_rt && fabs(ms1[j].im - s.im) <= iso_im_tol)
            return (long)j;
        }
        return -1;
      };

      // [baseline-fix v2] When envelope scoring is OFF, run the ORIGINAL algorithm end-to-end, not
      // just the original scoring. The previous attempt reverted only the score while keeping the new
      // candidate generation (gapped walk + lead 0..3), producing a THIRD algorithm: 2.52M MS1 traces
      // -> 1.21M precursors, 967,956 spectra, 7,308 peptides (baseline: 1.76M -> 956k, 741k, 9,103).
      // The original is: contiguous walk with BREAK-ON-FIRST-MISS, seed is always the starting point
      // (no lead offsets), monoisotope = furthest-LEFT partner reached, pick z by partner COUNT with
      // ties favouring the LOWEST z.
      if (!envelope_scoring)
      {
        int best_z = 0, best_n = 1;
        double best_mono = s.mz;
        vector<size_t> best_partners;
        for (int z = 1; z <= max_charge; ++z)
        {
          vector<size_t> partners;
          double mono = s.mz;
          for (int k = 1; k <= 5; ++k)          // lighter partners define the monoisotope
          {
            long j = findPartner(s.mz - k * ISO / z);
            if (j < 0) break;                   // BREAK on first miss (contiguous only)
            partners.push_back((size_t)j); mono = ms1[j].mz;
          }
          for (int k = 1; k <= 5; ++k)          // heavier partners add confidence
          {
            long j = findPartner(s.mz + k * ISO / z);
            if (j < 0) break;
            partners.push_back((size_t)j);
          }
          const int n = 1 + (int)partners.size();
          if (n > best_n) { best_n = n; best_z = z; best_mono = mono; best_partners = partners; }
        }
        Precursor_ pc0;
        pc0.trace_idx = seed;
        pc0.rt = s.rt; pc0.im = s.im;
        pc0.mono_mz = best_mono;
        pc0.charge = best_z;                    // 0 => unknown, default_charge applied later
        out.push_back(pc0);
        used[seed] = true;
        for (size_t j : best_partners) used[j] = true;   // de-isotope
        continue;
      }

      struct Hyp { int z; double mono; double score; vector<size_t> partners; };
      vector<Hyp> hyps;
      vector<double> theo, obs(4), obsn(4);

      for (int z = 1; z <= max_charge; ++z)
      {
        // [route-1] GAPPED envelopes: no break-on-first-miss (a sub-threshold M+1 with a present
        // M+2 previously killed the whole walk). Try the seed AND up to 3 isotopes below it as the
        // monoisotope candidate; averagine decides which one really is the mono.
        for (int lead = 0; lead <= 3; ++lead)
        {
          const double mono_mz = s.mz - lead * ISO / z;
          long mono_idx = (lead == 0) ? (long)seed : findPartner(mono_mz);
          if (lead > 0 && mono_idx < 0) continue;
          const double neutral = mono_mz * z - z * PROTON_;
          if (neutral < 400.0 || neutral > 6000.0) continue;

          vector<size_t> members;
          vector<long> iso_idx(4, -1);
          obs.assign(4, 0.0);
          obs[0] = ms1[mono_idx].intensity;
          iso_idx[0] = mono_idx;
          if ((size_t)mono_idx != seed) members.push_back((size_t)mono_idx);
          int found = 1;
          for (int k = 1; k <= 3; ++k)
          {
            long j = findPartner(mono_mz + k * ISO / z);
            if (j >= 0) { obs[k] = ms1[j].intensity; iso_idx[k] = j; ++found;
                          if ((size_t)j != seed) members.push_back((size_t)j); }
          }
          if (found < 2) continue;   // need >=2 isotopes before a shape can be scored

          averagineRatios(neutral, 4, theo);
          double mx = 0.0; for (double v : obs) mx = std::max(mx, v);
          if (mx <= 0) continue;
          for (int k = 0; k < 4; ++k) obsn[k] = obs[k] / mx;
          const double cos_sim = cosineSim(obsn, theo);

          double csum = 0.0; int cn = 0;
          for (int k = 1; k <= 3; ++k)
            if (iso_idx[k] >= 0) { csum += xicCorr(ms1[iso_idx[0]].xic, ms1[iso_idx[k]].xic, delta_rt); ++cn; }
          const double mean_corr = cn ? std::max(0.0, csum / cn) : 0.0;

          // Composite: shape AND co-elution must BOTH hold, WEIGHTED BY EVIDENCE.
          // [adv-fix] The naive product cos*corr + 0.01*n is statistically indefensible: a 2-isotope
          // hypothesis (cos 0.99 x corr 0.95 = 0.9605) BEATS 4 genuine isotopes (0.95 x 0.85 = 0.8475),
          // because the cosine of a 2-vector is nearly always high and the +0.01*n bonus (0.02 between
          // 2 and 4 isotopes) sits far below correlation noise. That biases us toward exactly the
          // sparse, weakly-supported envelopes we are trying to eliminate.
          // Fix: multiply by an evidence weight (n-1)/3 capped at 1 => 4 isotopes 1.00, 3 -> 0.67,
          // 2 -> 0.33. The counterexample now scores 0.808 vs 0.313 - correctly ordered.
          const double evidence = std::min(1.0, (found - 1) / 3.0);
          hyps.push_back({z, mono_mz, cos_sim * mean_corr * evidence, members});
        }
      }

      Precursor_ pc;
      pc.trace_idx = seed;
      pc.rt = s.rt; pc.im = s.im;

      if (hyps.empty())
      {
        pc.mono_mz = s.mz; pc.charge = 0;   // no envelope evidence at all -> unknown charge
        out.push_back(pc);
        used[seed] = true;
        continue;
      }

      // [baseline-fix] When envelope scoring is OFF we must still DE-ISOTOPE, otherwise every isotope
      // peak seeds its own precursor and the spectrum count explodes (measured: 741k -> 1,209,295
      // spectra, 9,103 -> 7,830 peptides). The off-branch therefore falls back to the ORIGINAL
      // partner-COUNT criterion with the original low-z tie-break, and still marks members used.
      if (!envelope_scoring)
      {
        for (auto& h : hyps) h.score = (double)h.partners.size();  // count, not shape
        sort(hyps.begin(), hyps.end(), [](const Hyp& a, const Hyp& b) {
          if (a.score != b.score) return a.score > b.score;
          if (a.z != b.z) return a.z < b.z;                        // original: ties favour LOW z
          return a.mono < b.mono;
        });
      }
      else
      sort(hyps.begin(), hyps.end(), [](const Hyp& a, const Hyp& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.z != b.z) return a.z > b.z;   // [route-1] ties now favour HIGHER z (old code favoured low z)
        return a.mono < b.mono;
      });
      const Hyp& best = hyps[0];
      pc.mono_mz = best.mono; pc.charge = best.z;
      out.push_back(pc);
      used[seed] = true;
      for (size_t j : best.partners) used[j] = true;

      // [route-3] Ambiguous charge -> retain a second hypothesis (diaTracer: "when the charge state
      // cannot be confidently determined, multiple values are retained"). 0 disables.
      if (ambig_margin > 0.0)
      {
        for (size_t h = 1; h < hyps.size(); ++h)
        {
          if (hyps[h].z == best.z) continue;
          if (best.score - hyps[h].score > ambig_margin) break;
          Precursor_ alt = pc;
          alt.mono_mz = hyps[h].mono; alt.charge = hyps[h].z;
          out.push_back(alt);
          break;                            // at most one alternative
        }
      }
    }
    return out;
  }

  /// [route-4] Collapse precursor hypotheses that are the SAME species (same charge, m/z within ppm,
  /// co-eluting in RT, co-located in IM). Measured redundancy is 4.6 identified spectra per unique
  /// peptide vs diaTracer's 1.46, so duplicates cost spectra, RAM and FDR budget without adding
  /// peptides. Keeps the most intense representative. Deterministic (total-order sort).
  static size_t dedupPrecursors_(vector<Precursor_>& pcs, const vector<Trace>& ms1,
                                 double mass_ppm, double rt_tol, double im_tol)
  {
    if (pcs.empty()) return 0;
    vector<size_t> idx(pcs.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
      if (pcs[a].charge != pcs[b].charge) return pcs[a].charge < pcs[b].charge;
      if (pcs[a].mono_mz != pcs[b].mono_mz) return pcs[a].mono_mz < pcs[b].mono_mz;
      if (pcs[a].rt != pcs[b].rt) return pcs[a].rt < pcs[b].rt;
      return a < b;
    });
    vector<bool> drop(pcs.size(), false);
    for (size_t a = 0; a < idx.size(); ++a)
    {
      size_t i = idx[a];
      if (drop[i]) continue;
      const double tol = pcs[i].mono_mz * mass_ppm * 1e-6;
      for (size_t b = a + 1; b < idx.size(); ++b)
      {
        size_t j = idx[b];
        if (pcs[j].charge != pcs[i].charge) break;
        if (pcs[j].mono_mz - pcs[i].mono_mz > tol) break;
        if (drop[j]) continue;
        if (fabs(pcs[j].rt - pcs[i].rt) <= rt_tol && fabs(pcs[j].im - pcs[i].im) <= im_tol)
        {
          if (ms1[pcs[j].trace_idx].intensity > ms1[pcs[i].trace_idx].intensity) { drop[i] = true; break; }
          drop[j] = true;
        }
      }
    }
    vector<Precursor_> keep;
    keep.reserve(pcs.size());
    size_t removed = 0;
    for (size_t i = 0; i < pcs.size(); ++i) { if (drop[i]) ++removed; else keep.push_back(pcs[i]); }
    pcs.swap(keep);
    return removed;
  }

  /// Resample the precursor XIC onto the window grid (once), then for every IM+RT+overlap-passing
  /// candidate fragment compute the Pearson correlation and, if >= min_corr, call emit(fi, corr).
  /// Shared by the per-precursor path and the competitive (best-precursor) assignment. [bench]
  template <class Emit>
  void scoreCandidates_(const Precursor_& pc, const vector<Trace>& frag_traces,
                        const vector<double>& frag_im, const vector<Trace>& ms1_traces,
                        const FragGrid& fg, double delta_im, double delta_rt, double min_corr,
                        int min_corr_pts, vector<float>& pdense, Emit&& emit) const
  {
    const double G = (double)fg.rt.size();
    const vector<pair<double, double>>& p_xic = ms1_traces[pc.trace_idx].xic;
    vector<int> touched;
    touched.reserve(p_xic.size());
    for (const auto& p : p_xic)
    {
      auto it = lower_bound(fg.rt.begin(), fg.rt.end(), p.first);
      int gi = -1; double bd = delta_rt;
      for (auto c : {it == fg.rt.begin() ? it : prev(it), it})
      {
        if (c == fg.rt.end()) continue;
        double d = fabs(*c - p.first);
        if (d <= bd) { bd = d; gi = (int)(c - fg.rt.begin()); }
      }
      if (gi >= 0) { if (pdense[gi] == 0.0f) touched.push_back(gi); pdense[gi] += (float)p.second; }
    }
    double psum = 0, psumsq = 0;
    for (int gi : touched) { double v = pdense[gi]; psum += v; psumsq += v * v; }
    double pmean = psum / G, pvar = psumsq - G * pmean * pmean;
    double pinv = pvar > 0 ? 1.0 / sqrt(pvar) : 0.0;
    if (pinv > 0.0)
    {
      size_t lo = lower_bound(frag_im.begin(), frag_im.end(), pc.im - delta_im) - frag_im.begin();
      size_t hi = upper_bound(frag_im.begin(), frag_im.end(), pc.im + delta_im) - frag_im.begin();
      for (size_t fi = lo; fi < hi; ++fi)
      {
        const Trace& f = frag_traces[fi];
        if (fabs(f.rt - pc.rt) > delta_rt) continue;                   // RT gate
        if (fg.invnorm[fi] == 0.0) continue;                           // degenerate fragment
        if (fabs(f.mz - pc.mono_mz) < 0.01) continue;                  // exclude precursor peak [M-7]
        double dot = 0; int overlap = 0;
        for (size_t k = fg.supBegin(fi), ke = fg.supEnd(fi); k < ke; ++k)
        { const auto& gv = fg.flat[k]; float pv = pdense[gv.first]; if (pv != 0.0f) { dot += (double)pv * gv.second; ++overlap; } }
        if (overlap < min_corr_pts) continue;                          // overlap guard [H-4]
        double c = (dot - G * pmean * fg.mean[fi]) * pinv * fg.invnorm[fi]; // Pearson
        if (c < min_corr) continue;                                    // correlation gate
        emit(fi, c);
      }
    }
    for (int gi : touched) pdense[gi] = 0.0f;                          // reset scratch
  }

  /// Build a pseudo-MS2 from a fragment list (m/z,intensity) + ranking scores: keeps top
  /// max_frags by score, drops if < min_frags, annotates the synthetic precursor. Deterministic.
  void assembleFromList_(const Precursor_& pc, double win_lo, double win_hi,
                         vector<pair<double, double>>& frags, vector<double>& frag_scores,
                         Size min_frags, Size max_frags, MSSpectrum& out) const
  {
    if (frags.size() < min_frags) return;
    if (frags.size() > max_frags)
    {
      vector<size_t> idx(frags.size());
      for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
      partial_sort(idx.begin(), idx.begin() + max_frags, idx.end(), [&](size_t a, size_t b) {
        if (frag_scores[a] != frag_scores[b]) return frag_scores[a] > frag_scores[b];
        if (frags[a].first != frags[b].first) return frags[a].first < frags[b].first;
        return a < b;
      });
      vector<pair<double, double>> keep; keep.reserve(max_frags);
      for (size_t i = 0; i < max_frags; ++i) keep.push_back(frags[idx[i]]);
      frags.swap(keep);
    }
    out.setMSLevel(2);
    out.setType(SpectrumSettings::SpectrumType::CENTROID);
    out.setRT(pc.rt);
    for (const auto& fr : frags) out.emplace_back(fr.first, fr.second);
    out.sortByPosition();
    OpenMS::Precursor prec;
    prec.setMZ(pc.mono_mz);
    if (pc.charge > 0) prec.setCharge(pc.charge);
    prec.setIsolationWindowLowerOffset(max(0.0, pc.mono_mz - win_lo));   // [C-5]
    prec.setIsolationWindowUpperOffset(max(0.0, win_hi - pc.mono_mz));
    prec.setDriftTime(pc.im);
    prec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    out.setPrecursors({prec});
  }

  /// Per-precursor assembly: this precursor claims every fragment passing its gate (a fragment may
  /// be shared across precursors -> chimeric). Leaves @p out empty on failure.
  void assembleOne_(const Precursor_& pc, double win_lo, double win_hi,
                    const vector<Trace>& frag_traces, const vector<double>& frag_im,
                    const vector<Trace>& ms1_traces, const FragGrid& fg,
                    double delta_im, double delta_rt, double min_corr, int min_corr_pts,
                    Size min_frags, Size max_frags, vector<float>& pdense, MSSpectrum& out) const
  {
    vector<pair<double, double>> frags;
    vector<double> frag_scores;
    scoreCandidates_(pc, frag_traces, frag_im, ms1_traces, fg, delta_im, delta_rt, min_corr,
                     min_corr_pts, pdense, [&](size_t fi, double c) {
      frags.emplace_back(frag_traces[fi].mz, frag_traces[fi].intensity);
      frag_scores.push_back(c * frag_traces[fi].intensity);
    });
    assembleFromList_(pc, win_lo, win_hi, frags, frag_scores, min_frags, max_frags, out);
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const double delta_im = getDoubleOption_("gate:delta_im");
    const double delta_rt = getDoubleOption_("gate:delta_rt");
    const double min_corr = getDoubleOption_("gate:min_correlation");
    const int min_corr_pts = getIntOption_("gate:min_correlation_points");
    const Size min_frags = (Size)getIntOption_("assembly:min_fragments");
    const Size max_frags = (Size)getIntOption_("assembly:max_fragments");
    const int max_charge = getIntOption_("max_charge");
    const double mass_ppm = getDoubleOption_("trace:mass_error_ppm");

    PeakMap exp;
    // Accept mzML, mzPeak, or a Bruker .d directory; FileHandler auto-detects the format
    // (mzPeak -> MzPeakFile; .d -> BrukerTimsFile, requires WITH_OPENTIMS).
    PeakMap ms1_map;
    map<pair<int, int>, vector<CompactFrame>> ms2_by_window;
    CompactStats cstat;
    auto winKey = [](double lo, double hi) {
      return make_pair((int)llround(lo * 100.0), (int)llround(hi * 100.0));
    };
    const bool stream_load = getFlag_("perf:stream_load");
    bool streamed = false;

    if (stream_load)
    {
      // [stream] frame-by-frame: pick + compact on arrival, never hold the whole run.
      // Also the ONLY place the reader's NATIVE pre-centroid frame aggregation is reachable.
      BrukerTimsFile::Config cfg;
      cfg.dia_ms2_n_neighbors = getIntOption_("trace:native_ms2_neighbors");
      cfg.ms1_n_neighbors     = getIntOption_("trace:native_ms1_neighbors");
      PeakPickerIM spicker;
      Param sp = spicker.getParameters();
      sp.setValue("pickIMCluster:im_tolerance_cluster", getDoubleOption_("gate:delta_im"));
      sp.setValue("pickIMCluster:ppm_tolerance_cluster", getDoubleOption_("trace:mass_error_ppm"));
      spicker.setParameters(sp);
      PickCompactConsumer consumer(spicker, ms2_by_window, ms1_map, cstat, winKey);
      try
      {
        BrukerTimsFile().loadDIAStreaming(in, consumer, cfg);
        streamed = true;
        writeLogInfo_("[stream] frame-by-frame load done: " + String(consumer.frames_seen)
                      + " frames, MS1 " + String(ms1_map.size()) + ", windows "
                      + String(ms2_by_window.size()) + " (native ms2_neighbors="
                      + String(cfg.dia_ms2_n_neighbors) + ")" + clk_() + rss_());
      }
      catch (const Exception::BaseException& e)
      {
        writeLogWarn_(String("[stream] streaming load failed (") + e.getName()
                      + "); falling back to loadExperiment. " + e.getMessage());
      }
    }

    if (!streamed)
      FileHandler().loadExperiment(in, exp, {FileTypes::MZML, FileTypes::MZPEAK, FileTypes::BRUKER_TDF}, log_type_);

    //-------------------------------------------------------------
    // Input normalization / validation [C-4]
    //-------------------------------------------------------------
    IMFormat imf = IMTypes::determineIMFormat(exp, 1);
    if (imf == IMFormat::NONE)
    {
      writeLogError_("Error: input has no ion-mobility data at MS1. DIAPseudoSpectra requires diaPASEF (IM DIA) input.");
      return ILLEGAL_PARAMETERS;
    }
    // ponytail: assumes 1/K0 (VSSC). ms/CCS conversion and FAIMS rejection are Phase 2 [M-2].

    // Centroid IM frames that are still profile, keeping per-peak IM at a tolerance no coarser
    // than the trace IM tolerance [C-2]. Split MS1 vs MS2 into separate working maps.
    PeakPickerIM picker;
    Param pp = picker.getParameters();
    pp.setValue("pickIMCluster:im_tolerance_cluster", delta_im);
    pp.setValue("pickIMCluster:ppm_tolerance_cluster", mass_ppm);
    picker.setParameters(pp);

    // Peak-pick (centroid) all IM frames in PARALLEL. PeakPickerIM is firstprivate so each thread
    // has its own copy -> the mutable ccs_warning_shown_ flag is not shared (a shared const picker
    // would be a real data race). pickIMCluster mutates only its own spectrum; distinct indices are
    // independent. [par-Crit-1]
    writeLogInfo_("Loaded " + String(exp.size()) + " frames." + clk_() + rss_()); // raw residency peak? [mem]
    writeLogInfo_("Peak-picking " + String(exp.size()) + " frames (OpenMP)...");
    // Param is a deep-copy value type (mutable ParamNode root_, defaulted copy) so firstprivate
    // gives each thread a fully independent picker -> no refcount/COW race. pickIMCluster can throw
    // (e.g. unexpected IM format); an exception escaping an OMP region terminates the process, so
    // capture it and rethrow serially. [par-Crit-1, code-review]
    std::exception_ptr pick_err;
    #pragma omp parallel for firstprivate(picker) schedule(dynamic, 16)
    for (long i = 0; i < (long)exp.size(); ++i)
    {
      try
      {
        MSSpectrum& s = exp[i];
        if (!s.containsIMData()) continue;
        if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) picker.pickIMCluster(s);
        ensureIMArrayName(s); // so MassTraceDetection reads per-peak IM (IM-aware tracing) [C-1]
      }
      catch (...)
      {
        #pragma omp critical
        if (!pick_err) pick_err = std::current_exception();
      }
    }
    if (pick_err) std::rethrow_exception(pick_err);

    // Serial dispatch: MOVE (not copy) frames into the MS1 map / per-window MS2 maps to avoid
    // transient memory doubling; container inserts aren't thread-safe so this stays serial.
    if (!streamed)
    for (MSSpectrum& s : exp)
    {
      if (!s.containsIMData()) continue;
      if (s.getMSLevel() == 1)
      {
        ms1_map.addSpectrum(std::move(s));
      }
      else if (s.getMSLevel() == 2 && !s.getPrecursors().empty())
      {
        const OpenMS::Precursor& pr = s.getPrecursors()[0];
        double lo = pr.getMZ() - pr.getIsolationWindowLowerOffset();
        double hi = pr.getMZ() + pr.getIsolationWindowUpperOffset();
        ms2_by_window[winKey(lo, hi)].push_back(compactify(s, cstat));
      }
    }
    exp.clear(true); // frames moved out; release the container
    {
      size_t cb = 0, cp = 0;
      for (const auto& kv : ms2_by_window) for (const auto& f : kv.second) { cb += f.bytes(); cp += f.mzq.size(); }
      const size_t lost = cstat.no_im_array + cstat.size_mismatch + cstat.bad_mz + cstat.bad_im;
      if (lost)
      {
        writeLogWarn_("[compact] DROPPED " + String(lost) + " peaks: no_im_array=" + String(cstat.no_im_array)
                      + " size_mismatch=" + String(cstat.size_mismatch) + " bad_mz=" + String(cstat.bad_mz)
                      + " bad_im=" + String(cstat.bad_im) + " (kept " + String(cstat.kept) + ") -- if this is"
                      + " nonzero the compact store is NOT lossless and results differ from the old path.");
      }
      writeLogInfo_("Split into MS1 + " + String(ms2_by_window.size()) + " windows (picked); compact store "
                    + String(cb / (1024ULL * 1024ULL)) + " MB for " + String(cp) + " peaks (~"
                    + String(cp ? (double)cb / cp : 0.0) + " B/peak)." + clk_() + rss_()); // [mem]
    }

    if (ms1_map.empty())
    {
      writeLogError_("Error: no MS1 ion-mobility frames found; cannot seed precursors.");
      return INCOMPATIBLE_INPUT_DATA;
    }

    //-------------------------------------------------------------
    // MS1 traces + precursor/charge inference (shared, once) [H-8,H-10]
    //-------------------------------------------------------------
    const int agg_ms1 = getIntOption_("trace:frame_aggregation_ms1_n");
    ms1_map.sortSpectra();
    if (agg_ms1 > 1) { aggregateFrames_(ms1_map, agg_ms1, mass_ppm, delta_im); writeLogInfo_("MS1 cross-frame aggregation N=" + String(agg_ms1) + rss_()); }
    const double ms1_snr = getDoubleOption_("trace:ms1_chrom_peak_snr");
    const double ms2_snr = getDoubleOption_("trace:ms2_chrom_peak_snr");
    const double ms2_minlen = getDoubleOption_("trace:ms2_min_length_sec");
    const double ms2_msr = getDoubleOption_("trace:ms2_min_sample_rate");
    const double ms2_split = getDoubleOption_("trace:ms2_split_valleys");
    const double max_trace_len = getDoubleOption_("trace:max_trace_length_sec");
    const int agg_n = getIntOption_("trace:frame_aggregation_n");
    String ms1_span;
    vector<Trace> ms1_traces = detectTraces_(ms1_map, delta_im, getDoubleOption_("trace:noise_threshold_int"),
                                             ms1_snr, getDoubleOption_("trace:min_length_sec"), getDoubleOption_("trace:ms1_min_sample_rate"), max_trace_len, getDoubleOption_("trace:ms1_split_valleys"), &ms1_span);
    writeLogInfo_("MS1 " + ms1_span);   // [merged-trace] measure, do not assume
    const double iso_im_tol = getDoubleOption_("charge:iso_im_tolerance");
    const double ms2_noise = getDoubleOption_("trace:ms2_noise_threshold_int");
    writeLogInfo_("MS2 trace gate: apex > " + String(ms2_snr * ms2_noise) + " (snr " + String(ms2_snr)
                  + " x noise " + String(ms2_noise) + "), membership > " + String(ms2_noise)
                  + ", min_length " + String(ms2_minlen) + " s");
    const bool env_scoring = (getStringOption_("charge:scoring") == "envelope");
    const double ambig_margin = getDoubleOption_("charge:ambiguity_margin");
    vector<Precursor_> precursors = inferPrecursors_(ms1_traces, max_charge, delta_rt, iso_im_tol,
                                                     mass_ppm, env_scoring, ambig_margin);
    if (getFlag_("assembly:dedup_precursors"))
    {
      const size_t before = precursors.size();
      const size_t removed = dedupPrecursors_(precursors, ms1_traces, mass_ppm, delta_rt, delta_im);
      writeLogInfo_("Precursor dedup: removed " + String(removed) + " of " + String(before)
                    + " hypotheses (" + String((int)(1000.0 * removed / std::max<size_t>(before, 1)) / 10.0) + "%)");
    }
    // Emulate diaTracer, which assigns a charge to every precursor: give isotope-unsupported
    // precursors the default charge (0 = keep the strict isotope-only behavior and drop them).
    const int default_charge = getIntOption_("assembly:default_charge");
    // rp_max = per-fragment precursor cap. 0 = share-all (stream, memory-light). competitive flag = rp_max 1.
    int rp_max = getIntOption_("assembly:rp_max");
    const double apportion = getDoubleOption_("assembly:apportion");
    if (getFlag_("assembly:competitive") && rp_max <= 0) rp_max = 1;
    Size n_default = 0;
    const bool open_safe = getFlag_("assembly:open_search_safe");
    // [way-4] A GUESSED charge and an isotope-SUPPORTED charge must not be emitted identically.
    // Closed search hides the difference (engine enumerates); open search does not.
    for (auto& pc : precursors)
      if (pc.charge == 0) { pc.guessed = true; pc.charge = open_safe ? 0 : default_charge; ++n_default; }
    // [ms1-funnel] Measured recall against a DIA-NN truth set is 74.3% vs diaTracer's 84.4%, but a
    // single number cannot say WHERE a precursor is lost: no MS1 trace at all, a trace that never
    // became a hypothesis, or a hypothesis that never produced a spectrum. Dump both upstream
    // stages so the loss can be attributed to one of them instead of guessed at.
    const String ms1_dump = getStringOption_("diag:dump_ms1_tsv");
    if (!ms1_dump.empty())
    {
      ofstream ft((ms1_dump + ".traces.tsv").c_str());
      ft << "mz\trt\tim\tintensity\n";
      for (const auto& t : ms1_traces)
        ft << t.mz << '\t' << t.rt << '\t' << t.im << '\t' << t.intensity << '\n';
      ft.close();
      ofstream fp((ms1_dump + ".precursors.tsv").c_str());
      fp << "mono_mz\tcharge\trt\tim\tguessed\n";
      for (const auto& pc : precursors)
        fp << pc.mono_mz << '\t' << pc.charge << '\t' << pc.rt << '\t' << pc.im << '\t'
           << (pc.guessed ? 1 : 0) << '\n';
      fp.close();
      writeLogInfo_("MS1 funnel dump: " + String(ms1_traces.size()) + " traces, "
                    + String(precursors.size()) + " precursors -> " + ms1_dump + ".{traces,precursors}.tsv");
    }
    // sort precursors by mono m/z so each window can binary-search its members instead of
    // rescanning all ~1M precursors per window (total order for determinism). [code-review]
    sort(precursors.begin(), precursors.end(), [](const Precursor_& a, const Precursor_& b) {
      if (a.mono_mz != b.mono_mz) return a.mono_mz < b.mono_mz;
      if (a.rt != b.rt) return a.rt < b.rt;
      if (a.im != b.im) return a.im < b.im;
      return a.trace_idx < b.trace_idx;
    });
    vector<double> prec_mz(precursors.size());
    for (size_t i = 0; i < precursors.size(); ++i) prec_mz[i] = precursors[i].mono_mz;
    // MEMORY: ms1_map is dead the moment its traces exist (nothing below reads it). Release it
    // rather than letting it ride to end of scope through the whole parallel phase. [mem]
    ms1_map.clear(true);
    // MEMORY: scoring reads ms1_traces[pc.trace_idx].xic, so ms1_traces must STAY -- but only the
    // traces an actual precursor points at are ever dereferenced. Free the xic of the rest
    // (2.5M traces vs 1.2M precursors => roughly half are dead weight held through scoring). [mem]
    {
      vector<bool> needed(ms1_traces.size(), false);
      for (const auto& pc : precursors)
        if (pc.trace_idx < ms1_traces.size()) needed[pc.trace_idx] = true;
      Size freed = 0;
      for (size_t i = 0; i < ms1_traces.size(); ++i)
        if (!needed[i] && !ms1_traces[i].xic.empty())
        {
          vector<pair<double, double>>().swap(ms1_traces[i].xic);
          ++freed;
        }
      writeLogInfo_("Released XICs of " + String(freed) + " unreferenced MS1 traces." + rss_());
    }
    writeLogInfo_(clk_() + " Detected " + String(ms1_traces.size()) + " MS1 traces -> " + String(precursors.size())
                  + " precursor hypotheses (" + String(n_default) + " assigned default charge)." + rss_());
    std::exception_ptr worker_err;

    //-------------------------------------------------------------
    // PARALLEL over windows (independent; with ~24 windows this fully uses the cores, and each
    // window's trace detection AND gating run on its own core). Each window frees its raw frames
    // right after trace extraction. Per-window buckets are merged in fixed window order then
    // canonically sorted, so output is thread-count-independent. Exceptions (e.g. MassTraceDetection
    // on a degenerate window) are captured and rethrown serially. [C-3,H-5,C-5,par-Crit-2/3/4]
    //-------------------------------------------------------------
    vector<pair<pair<int, int>, vector<CompactFrame>*>> window_list;
    window_list.reserve(ms2_by_window.size());
    for (auto& kv : ms2_by_window) window_list.emplace_back(kv.first, &kv.second);
    vector<vector<MSSpectrum>> win_out(window_list.size());
    // Fan-out histogram (only filled when rp_max>0): buckets = #precursors each fragment gates,
    // BEFORE capping. Answers "does this rp_max actually bite?" Buckets: 1,2,3,4,5,6-10,11-25,26-50,51+.
    vector<uint64_t> fanout_hist(9, 0);
    Size win_done = 0;
    // MEMORY vs SPEED: every concurrently-processed window holds its own frames + traces + grid, so
    // peak RSS scales with the number of windows in flight, not with the window count. Capping
    // concurrency trades wall time for RAM. 0 = unlimited (previous behaviour). [mem]
    const int max_conc = getIntOption_("perf:max_concurrent_windows");
    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    // [bands] The window count is FIXED BY THE ACQUISITION METHOD (24 tiles here), so
    // min(windows, threads) pinned us to 24 cores no matter the machine -- measured: 24
    // concurrent windows on a 224-core node ran at load ~20. The ceiling is now the thread
    // count, with each window split into m/z bands to supply the extra work units.
    int n_bands = getIntOption_("perf:trace_bands");
    if (n_bands <= 0)  // auto: enough bands that windows x bands covers every thread
      n_bands = (int)std::ceil((double)n_threads / std::max<size_t>(window_list.size(), 1));
    n_bands = std::max(1, n_bands);
    int n_conc = std::min<int>((int)window_list.size(), n_threads);
    if (max_conc > 0) n_conc = std::min(n_conc, max_conc);
    if (n_conc < 1) n_conc = 1;
#ifdef _OPENMP
    // outer loop over windows + inner loop over bands = two live levels
    if (n_bands > 1) omp_set_max_active_levels(2);
#endif
    writeLogInfo_("Parallelism: " + String(window_list.size()) + " windows x " + String(n_bands)
                  + " m/z bands = " + String((int)window_list.size() * n_bands) + " trace units for "
                  + String(n_threads) + " threads");
    // [dyn-mem] Thread count is the UPPER bound; free RAM is the real one. A static cap either
    // wastes cores on a big machine or OOMs on a shared one, and this node is shared. Admit each
    // window only when its projected footprint fits in currently-free RAM, re-read from
    // /proc/meminfo at each admission so the run yields to other jobs instead of fighting them.
    // Always admit when nothing is in flight, so the loop cannot deadlock on a single huge window.
    const double mem_frac = getDoubleOption_("perf:mem_fraction");
    std::atomic<size_t> inflight{0};
    // Projected footprint per window, from its own compact size. Multiplier measured on S23:
    // ~5.8 GB peak against ~0.56 GB compact (materialised PeakMap at 20 B/peak, MassTrace at
    // 168 B header/trace, valley-split batch, FragGrid). Conservative by design.
    auto project = [&](const vector<CompactFrame>& fr) {
      size_t b = 0; for (const auto& f : fr) b += f.bytes();
      return (size_t)(b * 10.0) + (size_t)256 * 1024 * 1024;
    };
    writeLogInfo_("Processing " + String(window_list.size()) + " isolation windows (OpenMP over windows, <="
                  + String(n_conc) + " concurrent, admission bounded by " + String((int)(mem_frac * 100))
                  + "% of free RAM)..." + rss_());

    #pragma omp parallel for schedule(dynamic) num_threads(n_conc)
    for (long wi = 0; wi < (long)window_list.size(); ++wi)
    {
      try
      {
        const double win_lo = window_list[wi].first.first / 100.0;
        const double win_hi = window_list[wi].first.second / 100.0;
        // [dyn-mem] wait for room before materialising anything
        const size_t need = project(*window_list[wi].second);
        for (;;)
        {
          bool go = false;
          #pragma omp critical(admit)
          {
            const size_t budget = (size_t)(availableBytes_() * mem_frac);
            // "at least one": an empty pipeline always admits, however large the window
            if (inflight.load() == 0 || inflight.load() + need <= budget) { inflight += need; go = true; }
          }
          if (go) break;
          std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
        struct Release { std::atomic<size_t>& c; size_t n; ~Release() { c -= n; } } rel{inflight, need};
        // [compact] materialise ONLY this window, and free its compact frames as we go
        PeakMap wmap = materializeWindow(*window_list[wi].second);
        for (MSSpectrum& s : wmap) s.setMSLevel(1); // MassTraceDetection traces MS1-level only
        wmap.sortSpectra();
        // [cross-frame] sum adjacent cycles for THIS window before detecting fragment traces
        if (agg_n > 1) aggregateFrames_(wmap, agg_n, mass_ppm, delta_im);
        String ms2_span;
        vector<Trace> frag_traces = detectTraces_(wmap, delta_im, ms2_noise, ms2_snr, ms2_minlen, ms2_msr, max_trace_len, ms2_split, &ms2_span, n_bands);
        // (compact frames already released inside materializeWindow) [par-Crit-3]
        #pragma omp critical
        {
          if (win_done == 0) writeLogInfo_("MS2 " + ms2_span); // [merged-trace]
          writeLogInfo_("  window " + String(++win_done) + "/" + String(window_list.size()));
        }

        frag_traces.erase(remove_if(frag_traces.begin(), frag_traces.end(), [](const Trace& t) {
          return !std::isfinite(t.im) || !std::isfinite(t.mz) || !std::isfinite(t.rt);
        }), frag_traces.end());
        if (frag_traces.empty()) continue;

        sort(frag_traces.begin(), frag_traces.end(), [](const Trace& a, const Trace& b) {
          if (a.im != b.im) return a.im < b.im;
          if (a.mz != b.mz) return a.mz < b.mz;
          if (a.rt != b.rt) return a.rt < b.rt;
          return a.intensity < b.intensity;
        });
        vector<double> frag_im(frag_traces.size());
        for (size_t i = 0; i < frag_traces.size(); ++i) frag_im[i] = frag_traces[i].im;

        // Precompute the common correlation grid once per window; pdense is reused across the
        // window's precursors (this window runs on one thread). [perf]
        FragGrid fg = buildFragGrid(frag_traces);
        // MEMORY: the XICs are now resampled into fg; scoring below only reads the SCALAR fields
        // (mz/rt/im/intensity) of each Trace. Releasing the per-trace xic vectors here is the single
        // biggest per-window win: they would otherwise stay live through the whole scoring loop, and
        // up to N windows do this concurrently (the threshold-dependent bulk of peak RSS). NOTE the
        // traces themselves must NOT be freed -- mz/intensity are still needed. [mem]
        for (auto& t : frag_traces) { vector<pair<double, double>>().swap(t.xic); }
        vector<float> pdense(fg.rt.size(), 0.0f);

        vector<MSSpectrum>& bucket = win_out[wi];
        // only this window's precursors (by mono m/z); NaN mono_mz sorts to the end, excluded.
        // (charge 0 is NOT dropped: emitted with charge UNSET so the engine searches its range.)
        size_t plo = lower_bound(prec_mz.begin(), prec_mz.end(), win_lo) - prec_mz.begin();
        size_t phi = upper_bound(prec_mz.begin(), prec_mz.end(), win_hi) - prec_mz.begin();
        if (apportion > 0.0 && rp_max <= 0)
        {
          // [route-4] INTENSITY APPORTIONMENT — the only assembly lever that is not a count knob.
          //
          // Route 3 measured that for ~80% (chance-corrected) of the peptides diaTracer finds and we
          // miss, we DO have a precursor at the right (RT, m/z, IM): we emit a spectrum there, but the
          // right fragments never arrive with enough weight. Under share-all a fragment is copied at
          // FULL intensity into every gated precursor (mean fan-out 6.45), so a chimeric contribution
          // is indistinguishable from a private one and the true owner is diluted in rank.
          //
          // Every count-based remedy is falsified (rank-pruning, competitive, min_corr, min_corr_pts —
          // all monotonic losses). So: keep the fragment in ALL precursors (count unchanged) but split
          // its INTENSITY by relative correlation, w_i = corr_i^p / sum_j corr_j^p. The owner keeps
          // most of the signal, chimeric copies are down-weighted, and because max_fragments ranks on
          // score*intensity the ranking changes without any fragment being removed.
          // p = apportion (1 = linear, higher = sharper). p->inf degenerates to competitive, which is
          // falsified, so keep p modest.
          const size_t NF = frag_traces.size();
          vector<float> wsum(NF, 0.0f);
          for (size_t pi = plo; pi < phi; ++pi)
          {
            const Precursor_& pc = precursors[pi];
            if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
            scoreCandidates_(pc, frag_traces, frag_im, ms1_traces, fg, delta_im, delta_rt, min_corr,
                             min_corr_pts, pdense, [&](size_t fi, double c) {
              wsum[fi] += (float)std::pow(std::max(c, 0.0), apportion);
            });
          }
          for (size_t pi = plo; pi < phi; ++pi)
          {
            const Precursor_& pc = precursors[pi];
            if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
            vector<pair<double, double>> frags;
            vector<double> frag_scores;
            scoreCandidates_(pc, frag_traces, frag_im, ms1_traces, fg, delta_im, delta_rt, min_corr,
                             min_corr_pts, pdense, [&](size_t fi, double c) {
              const double w = wsum[fi] > 0.0f ? std::pow(std::max(c, 0.0), apportion) / wsum[fi] : 1.0;
              const double inten = frag_traces[fi].intensity * w;
              frags.emplace_back(frag_traces[fi].mz, inten);
              frag_scores.push_back(c * inten);
            });
            MSSpectrum ms2;
            assembleFromList_(pc, win_lo, win_hi, frags, frag_scores, min_frags, max_frags, ms2);
            if (!ms2.empty()) bucket.push_back(std::move(ms2));
          }
        }
        else if (rp_max <= 0)
        {
          // Share-all (default, current best): each fragment goes to EVERY gated precursor.
          // Streamed per precursor -> memory-light.
          for (size_t pi = plo; pi < phi; ++pi)
          {
            const Precursor_& pc = precursors[pi];
            if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
            MSSpectrum ms2;
            assembleOne_(pc, win_lo, win_hi, frag_traces, frag_im, ms1_traces, fg,
                         delta_im, delta_rt, min_corr, min_corr_pts, min_frags, max_frags, pdense, ms2);
            if (!ms2.empty()) bucket.push_back(std::move(ms2));
          }
        }
        else
        {
          // Soft RP rank-pruning (rp_max>=1; 1 == competitive): keep each fragment in only its
          // top-rp_max best-correlating precursors. Pass A: per fragment, a bounded top-rp_max heap of
          // (corr,pi) + total gate count for the fan-out histogram. Pass B: assemble per winning precursor.
          const size_t NF = frag_traces.size();
          vector<vector<pair<double, long>>> topk(NF); // ascending by corr; back() = weakest kept
          vector<uint32_t> nfan(NF, 0);                 // total precursors that gated this fragment
          for (size_t pi = plo; pi < phi; ++pi)
          {
            const Precursor_& pc = precursors[pi];
            if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
            scoreCandidates_(pc, frag_traces, frag_im, ms1_traces, fg, delta_im, delta_rt, min_corr,
                             min_corr_pts, pdense, [&](size_t fi, double c) {
              ++nfan[fi];
              auto& tk = topk[fi];
              // insertion sort into a size-<=rp_max list kept ascending; deterministic tie-break on pi.
              if ((int)tk.size() < rp_max) {
                tk.emplace_back(c, (long)pi);
                for (size_t j = tk.size() - 1; j > 0 && (tk[j].first < tk[j-1].first ||
                     (tk[j].first == tk[j-1].first && tk[j].second < tk[j-1].second)); --j)
                  std::swap(tk[j], tk[j-1]);
              } else if (c > tk.front().first || (c == tk.front().first && (long)pi < tk.front().second)) {
                tk.front() = {c, (long)pi};
                for (size_t j = 1; j < tk.size() && (tk[j-1].first > tk[j].first ||
                     (tk[j-1].first == tk[j].first && tk[j-1].second > tk[j].second)); ++j)
                  std::swap(tk[j-1], tk[j]);
              }
            });
          }
          // fan-out histogram (before capping)
          vector<uint64_t> h(9, 0);
          for (size_t fi = 0; fi < NF; ++fi) {
            uint32_t n = nfan[fi];
            if (n == 0) continue;
            size_t b = n <= 5 ? n - 1 : n <= 10 ? 5 : n <= 25 ? 6 : n <= 50 ? 7 : 8;
            ++h[b];
          }
          // invert to per-precursor fragment lists (ordered by pi -> deterministic)
          map<long, vector<size_t>> won;
          for (size_t fi = 0; fi < NF; ++fi) for (auto& cp : topk[fi]) won[cp.second].push_back(fi);
          for (auto& kv2 : won)
          {
            const Precursor_& pc = precursors[kv2.first];
            vector<pair<double, double>> frags;
            vector<double> frag_scores;
            for (size_t fi : kv2.second)
            {
              // this fragment's corr for THIS precursor (from its kept top-k list)
              double c = -1.0;
              for (auto& cp : topk[fi]) if (cp.second == kv2.first) { c = cp.first; break; }
              frags.emplace_back(frag_traces[fi].mz, frag_traces[fi].intensity);
              frag_scores.push_back(c * frag_traces[fi].intensity);
            }
            MSSpectrum ms2;
            assembleFromList_(pc, win_lo, win_hi, frags, frag_scores, min_frags, max_frags, ms2);
            if (!ms2.empty()) bucket.push_back(std::move(ms2));
          }
          #pragma omp critical
          for (size_t b = 0; b < 9; ++b) fanout_hist[b] += h[b];
        }
      }
      catch (...)
      {
        #pragma omp critical
        if (!worker_err) worker_err = std::current_exception();
      }
    }
    if (worker_err) std::rethrow_exception(worker_err);

    if (rp_max > 0)
    {
      uint64_t tot = 0; for (auto v : fanout_hist) tot += v;
      if (tot > 0)
      {
        const char* lbl[9] = {"1","2","3","4","5","6-10","11-25","26-50","51+"};
        String msg = "Fan-out histogram (precursors per fragment, before rp_max=" + String(rp_max) + " cap): ";
        for (size_t b = 0; b < 9; ++b)
          msg += String(lbl[b]) + "=" + String(fanout_hist[b]) + "(" + String((int)(1000.0 * fanout_hist[b] / tot) / 10.0) + "%) ";
        writeLogInfo_(msg);
      }
    }

    writeLogInfo_("All windows done." + clk_() + rss_()); // parallel-phase peak? [mem]

    vector<MSSpectrum> all_out;
    {
      size_t n_tot = 0;
      for (const auto& bucket : win_out) n_tot += bucket.size();
      all_out.reserve(n_tot); // avoid repeated reallocation of ~1M spectrum objects
    }
    for (auto& bucket : win_out) { for (auto& s : bucket) all_out.push_back(std::move(s)); bucket.clear(); bucket.shrink_to_fit(); }

    // Canonical total-order sort (RT, precursor m/z, charge): thread-count-independent and free of
    // RT-tie ambiguity (RT-only sortSpectra is NOT canonical). [par-Crit-2]
    sort(all_out.begin(), all_out.end(), [](const MSSpectrum& a, const MSSpectrum& b) {
      if (a.getRT() != b.getRT()) return a.getRT() < b.getRT();
      const double amz = a.getPrecursors().empty() ? 0.0 : a.getPrecursors()[0].getMZ();
      const double bmz = b.getPrecursors().empty() ? 0.0 : b.getPrecursors()[0].getMZ();
      if (amz != bmz) return amz < bmz;
      const int ac = a.getPrecursors().empty() ? 0 : a.getPrecursors()[0].getCharge();
      const int bc = b.getPrecursors().empty() ? 0 : b.getPrecursors()[0].getCharge();
      if (ac != bc) return ac < bc;
      // fragment-sequence tiebreaks so equal-precursor spectra have a canonical order [code-review]
      if (a.size() != b.size()) return a.size() < b.size();
      for (Size i = 0; i < a.size(); ++i)
      {
        if (a[i].getMZ() != b[i].getMZ()) return a[i].getMZ() < b[i].getMZ();
        if (a[i].getIntensity() != b[i].getIntensity()) return a[i].getIntensity() < b[i].getIntensity();
      }
      return false;
    });

    //-------------------------------------------------------------
    // [route-1] FEATURE-LEVEL CONSOLIDATION
    //
    // We emit 4.6 PSMs per unique peptide; diaTracer emits 1.46. We are not producing WORSE spectra
    // (purity 8.3% vs their 4.0%) - we are producing the SAME peptide many times while missing others.
    // Emission is per precursor HYPOTHESIS (per trace), not per chromatographic FEATURE, so one real
    // feature can be emitted from several RT slices / IM sub-ranges / charge hypotheses.
    //
    // Consolidate: spectra whose precursors agree in m/z (ppm), RT and IM are ONE feature -> keep the
    // single best (most fragments, then most intense) and drop the rest. This changes emission
    // multiplicity WITHOUT touching fragment content, so it is orthogonal to every falsified
    // fragment-count knob. It should also relieve FDR competition among near-duplicate spectra.
    //
    // NOTE this is a DELETION, and deletion has repeatedly cost peptides here - hence default OFF and
    // a hard kill criterion (see TRACING.md): keep only if unique peptides rise >=3%.
    //-------------------------------------------------------------
    const double cons_rt = getDoubleOption_("consolidate:delta_rt");
    if (cons_rt > 0.0)
    {
      const double cons_ppm = getDoubleOption_("consolidate:mass_ppm");
      const double cons_im = getDoubleOption_("consolidate:delta_im");
      const bool cons_z = getFlag_("consolidate:same_charge_only");
      vector<MSSpectrum> kept;
      kept.reserve(all_out.size());
      vector<char> dropped(all_out.size(), 0);
      // all_out is already sorted by (RT, mz, charge); a feature is a run within cons_rt of each other
      for (size_t i = 0; i < all_out.size(); ++i)
      {
        if (dropped[i]) continue;
        size_t best = i;
        const auto& pi = all_out[i].getPrecursors();
        if (pi.empty()) { kept.push_back(std::move(all_out[i])); continue; }
        const double mz_i = pi[0].getMZ();
        const int z_i = pi[0].getCharge();
        const double im_i = all_out[i].getDriftTime();
        // scan forward while RT is still within the window (sorted by RT)
        for (size_t j = i + 1; j < all_out.size(); ++j)
        {
          if (all_out[j].getRT() - all_out[i].getRT() > cons_rt) break;
          if (dropped[j]) continue;
          const auto& pj = all_out[j].getPrecursors();
          if (pj.empty()) continue;
          if (fabs(pj[0].getMZ() - mz_i) > mz_i * cons_ppm * 1e-6) continue;
          if (cons_z && pj[0].getCharge() != z_i) continue;
          const double im_j = all_out[j].getDriftTime();
          if (im_i > 0 && im_j > 0 && fabs(im_j - im_i) > cons_im) continue;
          // same feature: keep whichever carries more fragments (tie -> higher TIC)
          if (all_out[j].size() > all_out[best].size() ||
              (all_out[j].size() == all_out[best].size() &&
               all_out[j].calculateTIC() > all_out[best].calculateTIC()))
          {
            dropped[best] = 1; best = j;
          }
          else dropped[j] = 1;
        }
        kept.push_back(std::move(all_out[best]));
        dropped[best] = 1;
      }
      writeLogInfo_("Feature consolidation: " + String(all_out.size()) + " -> " + String(kept.size())
                    + " spectra (" + String((int)(1000.0 * kept.size() / std::max<size_t>(all_out.size(), 1)) / 10.0)
                    + "% kept)" + rss_());
      all_out.swap(kept);
    }

    PeakMap out_exp;
    for (auto& s : all_out) out_exp.addSpectrum(std::move(s));
    const Size n_out = out_exp.size();
    addDataProcessing_(out_exp, getProcessingInfo_(DataProcessing::DATA_PROCESSING));
    FileHandler().storeExperiment(out, out_exp, {FileTypes::MZML}, log_type_);
    writeLogInfo_("Wrote " + String(n_out) + " pseudo-MS2 spectra to " + out);
    return EXECUTION_OK;
  }
};

/// @endcond

int main(int argc, const char** argv)
{
  TOPPDIAPseudoSpectra tool;
  return tool.main(argc, argv);
}
