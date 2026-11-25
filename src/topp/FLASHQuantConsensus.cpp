// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <future>
#include <map>
#include <numeric>
#include <tuple>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
  @page TOPP_FLASHQuantConsensus TOPP_FLASHQuantConsensus

  @brief Quantifies the MS1 features of several runs together: FLASHDeconv features matched across the runs (retention times aligned),
  and every matched feature measured the same way in every run from its raw signal, with match between runs.

  Input: the MS1 mzML files and the FLASHDeconv feature tsv files (-out) of the runs. The retention times of each run are aligned to the
  first run by the features most runs share; features of the same mass (10 ppm, or one isotope apart) and aligned apex (1 min) form a
  consensus feature. A consensus feature detected in at least -min_run_fraction of the runs is measured in every run, in its window (the
  aligned apex mapped into the run, +- half its median width): it is there if a raw mass trace at one of its charges (the median range
  of its detections) and core isotopes lies within 10 ppm with its apex in the window, and its value is then the raw signal of all those
  charges and isotopes: per MS1 scan of the window the raw peaks within 5 ppm of them (shifted by the run's mass error), integrated over
  retention time (within 10 ppm if none is within 5 ppm; the traces' areas if none is within 10 ppm). The same charges and isotopes in
  every run, whether a weak one has a trace there or not: trace areas lost the weak ones at low amounts and widened ratios. In a run
  where it was not detected (match between runs),
  its mass is looked for in every MS1 scan of the window with FLASHDeconv's peak group scoring (PeakGroup: raw peaks at its charges +- 2,
  isotope offset search, charge fit, isotope cosine, charge SNR, Qscore; on all peaks of the spectra as read, the noise floored at
  FLASHDeconv's detection limit of the scan: the weakest peak FLASHDeconv keeps), kept if its mass stays within the tolerance of the
  target (or one isotope off), and the scans combined
  as FLASHDeconv combines the scans of a feature (Qscore2D: one minus the lowest product of (1 - Qscore) over consecutive scans, so the
  whole envelope must recur scan after scan; a scan without peaks is passed over and does not break the run). The score of the
  transfer is that Qscore2D if it has a raw abundance and the mass found (intensity-weighted over the scans) recurs close to the
  feature's own mass (measured the same way where it was detected; the feature table's mass if no detected run has evidence): found in 3
  or more scans, within ten times the median mass error of the detected features (at most FLASHDeconv's tolerance); in 2 scans, within
  three times; in 1 scan, only for a mass of 8 kDa or more with a Qscore of at least 0.5, within three times (FLASHDeconv's tolerance
  if no feature has evidence in two runs where it was detected); else 0. Chance matches are nearly all of one scan of a small mass,
  while true transfers recur scan after scan, some of them a little further off than detected features. The null: a decoy mass
  below the envelope by its width + 5.25 isotopes (a shift within the envelope, as in bottom-up, falls on the species' own isotopes; a
  whole number of isotopes below it, on its other forms: oxidized, methylated, ammonia or water lost, all near whole daltons; half an
  isotope off still collects signal; a quarter collects the least), measured and
  scored the same way. The FDR of a score threshold is the number of null scores at or above it over the number of target scores at or
  above it; values whose q-value exceeds -mbr_fdr are NaN. A detected run without raw abundance (no trace in its window) is NaN too.
  Ion mode: that of the features (negative charges: negative mode, m/z = mass / charge - proton mass); the feature files of one call must
  share it. A run with fewer than 3 MS1 spectra with peaks has no mass traces: its values are NaN (with a warning); MS1 spectra without
  peaks (e.g. every other scan) are left out of the mass traces.
  Identifications (-in_id, optional): FLASHTnT PrSMs of the runs. They mark the transfers whose values are as reliable as detected
  ones (status mbr_id); they add few values. A transfer with a raw abundance whose run identifies its mass in its window (within three
  times the median mass error of the identifications of the detected features; FLASHDeconv's tolerance if none is identified) is
  identified. The identified transfers have an FDR of their own: the decoy masses the runs identify the same way (of all the transfers)
  over the identified transfers. Within -mbr_fdr, they are kept with the status mbr_id, whatever FLASHDeconv's evidence, and leave the
  targets of the other transfers; else they are judged like the others. The other transfers are judged by FLASHDeconv's evidence
  against all the decoys, as without identifications but without the identified transfers among their targets: their q-values are
  never lower than without identifications (at a strict -mbr_fdr, a few fewer pass), so the identifications never let another
  transfer through. What they give (PXD042169, -mbr_fdr 0.01): MS2 is acquired and identified on strong signals, which FLASHDeconv's
  evidence finds too (97-99% of the identified transfers are kept without identifications), so the kept values grow by 0.1-0.2%
  (SpikeIn 23420 -> 23438, ProteomeMix 17218 -> 17244). An mbr_id value is a little more exact than an mbr one (SpikeIn, E. coli,
  median |log2 error| 0.10 vs 0.18; off by more than 2-fold 0.1% vs 2.6%, detected values 1.3%); mbr_id is 5% of the kept transfers.
  It costs a FLASHTnT search of every run. -mbr_min_score keeps only the mbr transfers of at least that Qscore2D.
  -out_id (with -in_id) writes which PrSMs identify each consensus feature, in detected runs and transfers alike, so that the values
  can be read per proteoform or protein; it changes no value.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFLASHQuantConsensus : public TOPPBase
{
public:
  TOPPFLASHQuantConsensus() : TOPPBase("FLASHQuantConsensus", "Quantification of FLASHDeconv features across runs, with match between runs.")
  {
  }

protected:
  /// the file name without its path and its last suffix (as QFileInfo::completeBaseName)
  static std::string completeBaseName_(const std::string& path)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    const std::string name = path.substr(path.find_last_of("\\/") + 1);
#else
    const std::string name = path.substr(path.find_last_of('/') + 1);
#endif
    const auto dot = name.rfind('.');
    return dot == std::string::npos ? name : name.substr(0, dot);
  }

  struct Feature
  {
    double mass, apex, start, end, aligned; // minutes
    int min_charge, max_charge;
  };

  struct Consensus
  {
    double mass, aligned;
    std::map<Size, Size> runs; // run -> feature
  };

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", {},
                           "The runs as given to FLASHDeconv: mzML files (centroided), or Thermo .raw files (recognised by their extension; read "
                           "with Thermo's peak picking, as ThermoRawFileParser converts them by default; needs OpenMS built with WITH_THERMO_RAW and "
                           "the .NET 8 runtime), or Bruker timsTOF .d directories or zipped ones (.d.zip; recognised by their names; read with "
                           "OpenMS's timsTOF reader, summed over the ion mobility scans and centroided, as FLASHDeconv reads them; needs OpenMS "
                           "built with WITH_OPENTIMS, the default).", true);
    setValidFormats_("in", ListUtils::create<std::string>("mzML,raw,d"));
    registerInputFileList_("in_feature", "<files>", {}, "FLASHDeconv feature tsv files (-out) of the runs, in the order of -in.", true);
    setValidFormats_("in_feature", ListUtils::create<std::string>("tsv"));
    registerOutputFile_("out", "<file>", "", "Consensus feature tsv: per feature its mass, aligned apex and charges, and per run (named by the base name of its "
                        "-in file; a name two inputs share gets the run's position appended, _1, _2, ...) its abundance, "
                        "status (detected; mbr: a match between runs; mbr_id: one the run identifies, -in_id; or nan) and, for a match between "
                        "runs, its score (FLASHDeconv's Qscore2D of the mass, see above) and q-value (for mbr_id, that of the identified "
                        "transfers).", true);
    setValidFormats_("out", ListUtils::create<std::string>("tsv"));
    registerDoubleOption_("min_run_fraction", "<fraction>", 0.25, "Consensus features detected in at least this fraction of the runs are quantified in every run "
                          "(where not detected, by match between runs at -mbr_fdr). Lower values quantify more features, but their many weak transfers "
                          "raise the score threshold of all transfers (PXD042169: 0.25 doubled the quantified features over 0.5 at the same accuracy; "
                          "features in >= 2 runs lost transfers of the strong features).", false);
    setMinFloat_("min_run_fraction", 0.0);
    setMaxFloat_("min_run_fraction", 1.0);
    registerOutputFile_("out_mbr_detail", "<file>", "", "Match between runs, per consensus feature and run (detected or transfer): the value, FLASHDeconv's "
                        "evidence for the mass (Qscore2D, best scan Qscore, scans found, mass error), the same for its decoy where not detected, "
                        "the feature's own mass error where detected (the reference of the mass errors) and its median abundance where detected; "
                        "fields not computed (no evidence with -mbr_fdr < 0, no decoy for a detected run) are empty. With -in_id, also whether the "
                        "run identifies the mass (Identified) and the decoy mass (DecoyIdentified) of a transfer.", false, true);
    setValidFormats_("out_mbr_detail", ListUtils::create<std::string>("tsv"));
    registerDoubleOption_("mbr_fdr", "<fdr>", 0.01, "Match between runs: values are kept up to this FDR (target-decoy, see above); 1 keeps every "
                          "value, a negative value switches match between runs off (NaN where not detected).", false);
    registerInputFileList_("in_id", "<files>", {}, "Optional: FLASHTnT PrSM tsv files (-out_prsm) of the runs, in the order of -in. They mark the "
                           "transfers whose values are as reliable as detected ones (status mbr_id) and add few values (see above; PXD042169: "
                           "+0.1-0.2% values, as 97-99% of the identified transfers pass on FLASHDeconv's evidence anyway; values off by more "
                           "than 2-fold: mbr_id 0.1%, mbr 2.6%, detected 1.3%). A transfer whose run identifies its mass (a PrSM in its window: target protein, not a decoy precursor, proteoform q-value <= -id_fdr; its "
                           "precursor mass equal to the transfer's, or one isotope off, within three times the median error of the identifications "
                           "of the features where detected; FLASHDeconv's tolerance if none is identified) and with a raw abundance is kept "
                           "with the status mbr_id, whatever FLASHDeconv's evidence, if the identified transfers pass -mbr_fdr as a class: the "
                           "decoy masses identified the same way over the identified transfers; else they are judged like the others. The "
                           "other transfers are judged by FLASHDeconv's evidence without the kept identified ones among their targets: their "
                           "q-values are never lower than without identifications. Needs a FLASHTnT search of every run, with q-values (a "
                           "database with decoys). Without it, nothing changes.", false);
    setValidFormats_("in_id", ListUtils::create<std::string>("tsv"));
    registerOutputFile_("out_id", "<file>", "", "With -in_id: the identifications of the consensus features. One row per PrSM of -in_id "
                        "(target protein, not a decoy precursor, proteoform q-value <= -id_fdr) assigned to a consensus feature: the PrSM lies "
                        "in the feature's window in its run and its precursor mass is the feature's mass, or one isotope off, within "
                        "FLASHDeconv's tolerance (the nearest feature if several qualify). Columns: ConsensusIndex, Run (as in -out), the "
                        "PrSM's PrSMIndex, Scan, RetentionTime, ProteinAccession, ProteinDescription, ProteoformIndex, Proforma, "
                        "ProteoformMass, PrecursorMass, ProteoformLevelQvalue, and the MassErrorPPM of the match (nearest isotope). "
                        "Detected runs and transfers alike; a feature without a row is not identified. Without it, nothing changes.", false);
    setValidFormats_("out_id", ListUtils::create<std::string>("tsv"));
    registerDoubleOption_("id_fdr", "<q-value>", 0.01, "With -in_id: the PrSMs used, by their proteoform-level q-value.", false);
    registerDoubleOption_("mbr_min_score", "<score>", 0.0, "Match between runs: transfers other than mbr_id (status mbr) are kept only with a "
                          "score (Qscore2D) at least this, besides -mbr_fdr; mbr_id ones always pass. The FDR says that the mass is "
                          "there, not how exact the value is, and a higher score does not make it more exact: PXD042169 SpikeIn, E. coli "
                          "values off by more than 2-fold: mbr of Qscore2D >= 0.95 2.9%, < 0.95 2.3% (mbr_id 0.1%, detected 1.3%); "
                          "0.5 / 0.8 / 0.95 kept 97 / 82 / 59% of the transfers.", false);
    setMinFloat_("mbr_min_score", 0.0);
    setMaxFloat_("mbr_min_score", 1.0);
    setMinFloat_("id_fdr", 0.0);
    setMaxFloat_("id_fdr", 1.0);
  }

  /// the identified masses of a FLASHTnT PrSM tsv, (retention time in s, mass), by retention time: target PrSMs (no DECOY accession, no
  /// decoy precursor) whose proteoform-level q-value is at most @p id_fdr (-1, a q-value not computed, is not used; counted in
  /// @p without_q); the mass is the precursor mass FLASHDeconv measured (the proteoform mass when there is none)
  std::vector<std::pair<double, double>> readIdentifications_(const std::string& file, double id_fdr, Size& without_q) const
  {
    without_q = 0;
    std::ifstream in(file);
    std::string line;
    std::getline(in, line);
    std::vector<std::string> header;
    StringUtils::split(line, '\t', header);
    std::map<std::string, Size> col;
    for (Size i = 0; i < header.size(); ++i) col[StringUtils::trim(header[i])] = i;
    for (const std::string& c : {"RetentionTime", "ProteinAccession", "PrecursorMass", "ProteoformMass", "ProteoformLevelQvalue"})
    {
      if (col.find(c) == col.end()) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file, std::string("no column ") + c);
    }
    const auto decoy_weight = col.find("DecoyPrecursorWeight");
    std::vector<std::pair<double, double>> ids;
    while (std::getline(in, line))
    {
      std::vector<std::string> v;
      StringUtils::split(line, '\t', v);
      if (v.size() < header.size() || StringUtils::hasSubstring(v[col["ProteinAccession"]], "DECOY")) continue;
      // a decoy precursor PrSM keeps its target protein's accession (FLASHTnT -tnt:decoy_precursor with -tnt:keep_decoy)
      if (decoy_weight != col.end() && !StringUtils::trim(v[decoy_weight->second]).empty() && StringUtils::toDouble(v[decoy_weight->second]) != 0) continue;
      const double q = StringUtils::toDouble(StringUtils::trim(v[col["ProteoformLevelQvalue"]]));
      if (q < 0) { ++without_q; continue; }
      if (q > id_fdr) continue;
      double mass = StringUtils::toDouble(v[col["PrecursorMass"]]);
      if (!(mass > 0)) mass = StringUtils::toDouble(v[col["ProteoformMass"]]);
      if (mass > 0) ids.emplace_back(StringUtils::toDouble(v[col["RetentionTime"]]), mass);
    }
    std::sort(ids.begin(), ids.end());
    return ids;
  }

  /// a PrSM of a FLASHTnT PrSM tsv, as -out_id reports it
  struct PrSM
  {
    double rt, mass; // seconds; the precursor mass FLASHDeconv measured (the proteoform mass when there is none)
    std::vector<std::string> fields; // PrSMIndex, Scan, RetentionTime, ProteinAccession, ProteinDescription, ProteoformIndex, Proforma, ProteoformMass, PrecursorMass, ProteoformLevelQvalue
  };

  /// the PrSMs of a FLASHTnT PrSM tsv that -out_id assigns: those readIdentifications_ uses (target protein, not a decoy precursor,
  /// proteoform-level q-value in [0, @p id_fdr]), with the columns -out_id reports (empty where the file has no such column)
  std::vector<PrSM> readPrSMs_(const std::string& file, double id_fdr) const
  {
    static const std::vector<std::string> reported = {"PrSMIndex", "Scan", "RetentionTime", "ProteinAccession", "ProteinDescription", "ProteoformIndex",
                                                 "Proforma", "ProteoformMass", "PrecursorMass", "ProteoformLevelQvalue"};
    std::ifstream in(file);
    std::string line;
    std::getline(in, line);
    std::vector<std::string> header;
    StringUtils::split(line, '\t', header);
    std::map<std::string, Size> col;
    for (Size i = 0; i < header.size(); ++i) col[StringUtils::trim(header[i])] = i;
    const auto decoy_weight = col.find("DecoyPrecursorWeight");
    std::vector<PrSM> prsms;
    while (std::getline(in, line))
    {
      std::vector<std::string> v;
      StringUtils::split(line, '\t', v);
      if (v.size() < header.size() || StringUtils::hasSubstring(v[col["ProteinAccession"]], "DECOY")) continue;
      if (decoy_weight != col.end() && !StringUtils::trim(v[decoy_weight->second]).empty() && StringUtils::toDouble(v[decoy_weight->second]) != 0) continue;
      const double q = StringUtils::toDouble(StringUtils::trim(v[col["ProteoformLevelQvalue"]]));
      if (q < 0 || q > id_fdr) continue;
      double mass = StringUtils::toDouble(v[col["PrecursorMass"]]);
      if (!(mass > 0)) mass = StringUtils::toDouble(v[col["ProteoformMass"]]);
      if (!(mass > 0)) continue;
      PrSM p{StringUtils::toDouble(v[col["RetentionTime"]]), mass, {}};
      for (const auto& c : reported)
      {
        const auto it = col.find(c);
        p.fields.push_back(it == col.end() ? std::string("") : StringUtils::trimmed(v[it->second]));
      }
      prsms.push_back(p);
    }
    return prsms;
  }

  /// the FLASHDeconv MS1 features of a feature tsv, and in @p source the mzML it was made from (its FileName column; empty without rows)
  std::vector<Feature> readFeatures_(const std::string& file, std::string& source) const
  {
    source.clear();
    std::ifstream in(file);
    std::string line;
    std::getline(in, line);
    std::vector<std::string> header;
    StringUtils::split(line, '\t', header);
    std::map<std::string, Size> col;
    for (Size i = 0; i < header.size(); ++i) col[header[i]] = i;
    for (const std::string& c : {"MSLevel", "MonoisotopicMass", "StartRetentionTime", "EndRetentionTime", "ApexRetentionTime", "MinCharge", "MaxCharge"})
    {
      if (col.find(c) == col.end()) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file, std::string("no column ") + c);
    }
    const auto decoy_col = col.find("IsDecoy"); // decoy features (FLASHDeconv with FDR reporting) are not features
    std::vector<Feature> features;
    while (std::getline(in, line))
    {
      std::vector<std::string> v;
      StringUtils::split(line, '\t', v);
      if (source.empty() && col.find("FileName") != col.end() && v.size() > col["FileName"]) source = v[col["FileName"]];
      if (v.size() < header.size() || v[col["MSLevel"]] != "1" || (decoy_col != col.end() && v[decoy_col->second] != "0")) continue;
      features.push_back({StringUtils::toDouble(v[col["MonoisotopicMass"]]), StringUtils::toDouble(v[col["ApexRetentionTime"]]), StringUtils::toDouble(v[col["StartRetentionTime"]]),
                          StringUtils::toDouble(v[col["EndRetentionTime"]]), StringUtils::toDouble(v[col["ApexRetentionTime"]]), StringUtils::toInt32(v[col["MinCharge"]]),
                          StringUtils::toInt32(v[col["MaxCharge"]])});
    }
    return features;
  }

  static bool sameMass_(double m1, double m2)
  {
    const double iso = Constants::ISOTOPE_MASSDIFF_55K_U, tol = 10e-6 * m1;
    return std::abs(m1 - m2) <= tol || std::abs(m1 - m2 - iso) <= tol || std::abs(m1 - m2 + iso) <= tol;
  }

  /// features of the same mass and retention time (by @p rt) across runs, greedily in mass order
  static std::vector<Consensus> match_(const std::vector<std::vector<Feature>>& runs, bool aligned, double rt_tol)
  {
    std::vector<std::tuple<double, Size, Size>> all;
    for (Size r = 0; r < runs.size(); ++r)
    {
      for (Size i = 0; i < runs[r].size(); ++i) all.emplace_back(runs[r][i].mass, r, i);
    }
    std::sort(all.begin(), all.end());
    std::vector<Consensus> cons;
    for (const auto& [m, r, i] : all)
    {
      const double rt = aligned ? runs[r][i].aligned : runs[r][i].apex;
      bool placed = false;
      for (Size k = cons.size() > 300 ? cons.size() - 300 : 0; k < cons.size(); ++k)
      {
        auto& c = cons[k];
        if (c.runs.count(r) == 0 && sameMass_(c.mass, m) && std::abs(c.aligned - rt) <= rt_tol)
        {
          c.runs[r] = i;
          placed = true;
          break;
        }
      }
      if (!placed) cons.push_back({m, rt, {{r, i}}});
    }
    return cons;
  }

  /// FLASHDeconv's evidence for a mass in a retention time window
  struct FDEvidence
  {
    double qscore2d = 0, max_qscore = 0, ppm = 0; // ppm: intensity-weighted signed error of the scans' masses to the mass (nearest isotope)
    int scans = 0;                                 // scans where the mass was found
  };

  /// FLASHDeconv's evidence for @p target: in every MS1 scan of its window, the peak group of its mass as FLASHDeconv forms and scores
  /// one (raw peaks recruited at its charges +- 2, isotope offset search, charge fit, isotope cosine >= @p min_cos, charge SNR >= @p min_snr,
  /// Qscore; on all peaks of the spectra as read, the noise floored at @p detection_limits per scan), found if its mass stays within @p tol
  /// of the target's (or one isotope off); the scans combined as FLASHDeconv combines the scans of a feature (Qscore2D: one minus the
  /// lowest product of (1 - Qscore) over consecutive scans)
  static FDEvidence fdEvidence_(const PeakMap& map, const std::vector<double>& rts, const std::vector<float>& detection_limits,
                                const FLASHQuantAlgorithm::QuantTarget& target, const FLASHHelperClasses::PrecalculatedAveragine& avg, double tol,
                                double min_cos, double min_snr)
  {
    FDEvidence e;
    const double iso = Constants::ISOTOPE_MASSDIFF_55K_U;
    int zmin, zmax;
    const bool positive = target.chargeRange(zmin, zmax);
    const int z1 = std::max(1, zmin - 2), z2 = zmax + 2;
    double chain = 1, lowest = 1, weights = 0, ppm_sum = 0;
    bool in_chain = false;
    std::vector<FLASHHelperClasses::LogMzPeak> noisy;
    for (Size s = std::lower_bound(rts.begin(), rts.end(), target.rt_start) - rts.begin(); s < map.size() && rts[s] <= target.rt_end; ++s)
    {
      const MSSpectrum& spec = map[s];
      // no peaks: the scan is passed over and does not break the run of consecutive scans, as FLASHDeconv's feature tracing passes over
      // the scans without masses
      if (spec.empty()) continue;
      PeakGroup pg(z1, z2, positive);
      pg.setIsotopeDaDistance(iso);
      pg.setDetectionLimit(detection_limits[s]);
      double next = target.mass, prev_cos = 0;
      int offset = 0;
      bool determined = false;
      for (int k = 0; k < 10; ++k) // as SpectralDeconvolution: recruit and score until the isotope offset settles
      {
        pg.recruitAllPeaksInSpectrum(spec, tol, avg, next, noisy);
        if (pg.empty()) break;
        auto [a, b] = pg.getAbsChargeRange();
        offset = pg.updateQscore(noisy, avg, min_cos, tol, (a + b) < 20, {}, false);
        if (prev_cos > pg.getIsotopeCosine() || offset == 0 || k == 9)
        {
          if (pg.getChargeSNR(pg.getRepAbsCharge()) < min_snr) break;
          if (offset != 0) pg.recruitAllPeaksInSpectrum(spec, tol, avg, pg.getMonoMass() + offset * iso, noisy);
          std::tie(a, b) = pg.getAbsChargeRange();
          pg.updateQscore(noisy, avg, min_cos, tol, (a + b) < 20, {}, true);
          determined = true;
          break;
        }
        prev_cos = pg.getIsotopeCosine();
        next = pg.getMonoMass() + offset * iso;
      }
      double error = 1e300;
      for (int i = -1; i <= 1 && determined; ++i)
      {
        const double d = (pg.getMonoMass() - target.mass - i * iso) / target.mass * 1e6;
        if (std::abs(d) < std::abs(error)) error = d;
      }
      const double q = (determined && std::abs(error) <= tol * 1e6) ? pg.getQscore() : 0.0;
      if (q <= 0)
      {
        in_chain = false;
        continue;
      }
      chain = in_chain ? chain * (1 - q) : 1 - q;
      in_chain = true;
      lowest = std::min(lowest, chain);
      e.scans++;
      e.max_qscore = std::max(e.max_qscore, q);
      weights += pg.getIntensity();
      ppm_sum += pg.getIntensity() * error;
    }
    e.qscore2d = 1 - lowest;
    e.ppm = weights > 0 ? ppm_sum / weights : 0;
    return e;
  }

  /// FLASHDeconv's detection limit of a centroided spectrum: the weakest peak its low-peak filter keeps (above the intensity of the
  /// 30,000th peak and above 1/1000 of the highest); it floors the noise power of a peak group. The peaks below it are not removed here:
  /// with them filtered as FLASHDeconv does, PXD042169 kept fewer transfers (SpikeIn 12,252 vs 12,674 with all peaks).
  static float detectionLimit_(const MSSpectrum& spec)
  {
    if (spec.empty()) return 0;
    std::vector<float> ints;
    ints.reserve(spec.size());
    float max_int = 0;
    for (const auto& p : spec)
    {
      ints.push_back(p.getIntensity());
      max_int = std::max(max_int, p.getIntensity());
    }
    double threshold = max_int / 1000.0;
    const Size count = 30000; // FLASHDeconv's peak count for centroided spectra
    if (ints.size() > count)
    {
      std::nth_element(ints.begin(), ints.begin() + (SignedSize)count, ints.end(), std::greater<float>());
      threshold = std::max(threshold, (double)ints[count]);
    }
    float limit = 0;
    for (const auto& p : spec)
    {
      if (p.getIntensity() > threshold && (limit == 0 || p.getIntensity() < limit)) limit = p.getIntensity();
    }
    return limit;
  }

  static double median_(std::vector<double> v)
  {
    if (v.empty()) return 0;
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
  }

  /// the median of @p v, the mean of the two middle values for an even count (median_ takes the upper one)
  static double centreMedian_(std::vector<double> v)
  {
    if (v.empty()) return 0;
    const Size h = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + h, v.end());
    if (v.size() % 2 == 1) return v[h];
    return 0.5 * (v[h] + *std::max_element(v.begin(), v.begin() + h));
  }

  /// the retention time offset of a run to the first one at retention time @p rt: the median offset of the anchor features within 2 min
  static double offset_(const std::vector<std::pair<double, double>>& anchors, double rt)
  {
    const auto lo = std::lower_bound(anchors.begin(), anchors.end(), std::make_pair(rt - 2.0, -1e300));
    const auto hi = std::upper_bound(anchors.begin(), anchors.end(), std::make_pair(rt + 2.0, 1e300));
    std::vector<double> v;
    for (auto it = lo; it != hi; ++it) v.push_back(it->second);
    if (v.empty())
    {
      for (const auto& a : anchors) v.push_back(a.second);
    }
    return median_(v);
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList in = getStringList_("in"), in_feature = getStringList_("in_feature");
    const std::string out = getStringOption_("out");
    const double min_fraction = getDoubleOption_("min_run_fraction"), mbr_fdr = getDoubleOption_("mbr_fdr"), mbr_min_score = getDoubleOption_("mbr_min_score");
    if (in.size() != in_feature.size() || in.empty())
    {
      OPENMS_LOG_ERROR << "-in and -in_feature need the same number of files" << endl;
      return ILLEGAL_PARAMETERS;
    }
    const StringList in_id = getStringList_("in_id");
    if (!in_id.empty() && in_id.size() != in.size())
    {
      OPENMS_LOG_ERROR << "-in_id needs one file per run (as -in), or none" << endl;
      return ILLEGAL_PARAMETERS;
    }
    const std::string out_id = getStringOption_("out_id");
    if (!out_id.empty() && in_id.empty())
    {
      OPENMS_LOG_ERROR << "-out_id needs the identifications of the runs (-in_id)" << endl;
      return ILLEGAL_PARAMETERS;
    }
    const Size run_count = in.size();
    std::vector<std::vector<Feature>> runs(run_count);
    bool positive_mode = false, negative_mode = false;
    for (Size r = 0; r < run_count; ++r)
    {
      std::string source;
      runs[r] = readFeatures_(in_feature[r], source);
      // the feature table names the mzML FLASHDeconv read: a different name hints at files given in another order than -in
      // (a path written on Windows has backslashes, which QFileInfo on other systems does not split at)
      const auto base_name = [](std::string path) { StringUtils::substitute(path, '\\', '/'); return completeBaseName_(path); };
      if (!source.empty() && base_name(source) != base_name(in[r]))
      {
        OPENMS_LOG_WARN << "Warning: " << in_feature[r] << " was made from " << source << ", not from " << in[r] << " (the files of -in and -in_feature pair up by their order)" << endl;
      }
      for (const auto& f : runs[r]) (f.max_charge > 0 ? positive_mode : negative_mode) = true;
    }
    if (positive_mode && negative_mode)
    {
      OPENMS_LOG_ERROR << "The feature files hold features of both ion modes (positive and negative charges); quantify the runs of each mode separately" << endl;
      return ILLEGAL_PARAMETERS;
    }

    // retention time alignment to the first run: anchors are the features matched (3 min, unaligned) in >= 80% of the runs and the first
    std::vector<std::vector<std::pair<double, double>>> anchors(run_count); // per run: (its apex, offset to the first run), by apex
    for (const auto& c : match_(runs, false, 3.0))
    {
      if (c.runs.size() < 0.8 * run_count || c.runs.count(0) == 0) continue;
      const double ref = runs[0][c.runs.at(0)].apex;
      for (const auto& [r, i] : c.runs) anchors[r].emplace_back(runs[r][i].apex, runs[r][i].apex - ref);
    }
    for (Size r = 0; r < run_count; ++r)
    {
      std::sort(anchors[r].begin(), anchors[r].end());
      for (auto& f : runs[r]) f.aligned = anchors[r].empty() ? f.apex : f.apex - offset_(anchors[r], f.apex);
    }
    std::vector<Consensus> cons = match_(runs, true, 1.0);
    std::vector<Consensus> kept;
    for (auto& c : cons)
    {
      if (!c.runs.empty() && (double)c.runs.size() / run_count >= min_fraction - 1e-9) kept.push_back(std::move(c));
    }
    OPENMS_LOG_INFO << "Consensus features: " << cons.size() << ", detected in at least " << min_fraction * 100 << "% of the " << run_count
                    << " runs: " << kept.size() << endl;

    // per consensus feature: median mass, aligned apex, half width and charge range of its detections
    struct Summary { double mass, aligned, half; int zmin, zmax; };
    std::vector<Summary> summary;
    double max_mass = 0;
    for (const auto& c : kept)
    {
      std::vector<double> m, a, h, z1, z2;
      for (const auto& [r, i] : c.runs)
      {
        const auto& f = runs[r][i];
        m.push_back(f.mass);
        a.push_back(f.aligned);
        h.push_back((f.end - f.start) / 2);
        // the median of the absolute charges (FLASHDeconv writes -largest .. -smallest in negative mode), signed again below
        z1.push_back(negative_mode ? -f.max_charge : f.min_charge);
        z2.push_back(negative_mode ? -f.min_charge : f.max_charge);
      }
      const int zlo = (int)median_(z1), zhi = (int)median_(z2);
      summary.push_back({median_(m), median_(a), std::max(0.5, median_(h)), negative_mode ? -zhi : zlo, negative_mode ? -zlo : zhi});
      max_mass = std::max(max_mass, summary.back().mass);
    }
    SpectralDeconvolution sd;
    Param sd_param = sd.getDefaults();
    sd_param.setValue("max_mass", std::max(1000.0, max_mass * 1.1 + 100));
    sd.setParameters(sd_param);
    sd.calculateAveragine(false);
    const auto& avg = sd.getAveragine();
    // FLASHDeconv's MS1 tolerance, minimum isotope cosine and minimum charge SNR (its defaults)
    const double fd_tol = sd_param.getValue("tol").toDoubleVector()[0] * 1e-6, fd_min_cos = sd_param.getValue("min_cos").toDoubleVector()[0],
                 fd_min_snr = sd_param.getValue("min_snr").toDoubleVector()[0];
    // the decoy mass of a feature: below its isotope envelope (core isotopes, >= 10% of the most abundant) by the envelope's width + 5.25
    // isotopes, where neither the species nor its other forms have signal. A whole number of isotopes (width + 5 Da before) falls on its
    // oxidized, methylated, ammonia- or water-loss forms (-14 to -18 Da), and half an isotope still collects signal: PXD042169
    // ProteomeMix, decoys with FLASHDeconv evidence at width + 5.0 / 5.25 / 5.4 / 5.5 / 5.6 / 5.75 / 6.0 isotopes: 16095 / 1403 / 1468 /
    // 5408 / 1367 / 1080 / 14925 (transfers kept at 1% FDR 8161 / 13739 / 13739 / 10426 / 13739 / 13739 / 8020).
    std::vector<double> decoy_shift(kept.size());
    for (Size k = 0; k < kept.size(); ++k)
    {
      const auto iso = avg.get(summary[k].mass);
      double iso_max = 0;
      for (Size i = 0; i < iso.size(); ++i) iso_max = std::max(iso_max, (double)iso[i].getIntensity());
      int first = -1, last = -1;
      for (Size i = 0; i < iso.size(); ++i)
      {
        if (iso[i].getIntensity() < 0.1 * iso_max) continue;
        if (first < 0) first = (int)i;
        last = (int)i;
      }
      decoy_shift[k] = (last - first + 1 + 5.25) * Constants::ISOTOPE_MASSDIFF_55K_U;
    }

    // every run: the features and, where not detected, their decoys; their raw abundance and FLASHDeconv's evidence for their masses
    std::vector<std::vector<FLASHQuantAlgorithm::MassAbundance>> value(kept.size(), std::vector<FLASHQuantAlgorithm::MassAbundance>(run_count));
    std::vector<std::vector<double>> decoy_value(kept.size(), std::vector<double>(run_count, 0));
    std::vector<std::vector<std::array<FDEvidence, 2>>> fd_evidence(kept.size(), std::vector<std::array<FDEvidence, 2>>(run_count)); // target, decoy
    // -in_id: per feature and run, the smallest |mass error| (ppm, nearest isotope) of an identification of its mass in its window, and
    // of its decoy mass (transfers only); and whether it is identified within the gate (set below)
    const double id_fdr = getDoubleOption_("id_fdr");
    std::vector<std::vector<std::array<double, 2>>> id_error(kept.size(), std::vector<std::array<double, 2>>(run_count, {-1.0, -1.0}));
    std::vector<std::vector<std::array<bool, 2>>> identified(kept.size(), std::vector<std::array<bool, 2>>(run_count, {false, false}));
    // the smallest |error| (ppm) of an identified mass to @p mass (or one isotope off) in [@p rt_start, @p rt_end]; -1 without
    const auto id_match = [&](const std::vector<std::pair<double, double>>& ids, double mass, double rt_start, double rt_end) {
      double best = -1;
      for (auto it = std::lower_bound(ids.begin(), ids.end(), std::make_pair(rt_start, -1e300)); it != ids.end() && it->first <= rt_end; ++it)
      {
        for (int i = -1; i <= 1; ++i)
        {
          const double e = std::abs(it->second - mass - i * Constants::ISOTOPE_MASSDIFF_55K_U) / mass * 1e6;
          if (e <= fd_tol * 1e6 && (best < 0 || e < best)) best = e;
        }
      }
      return best;
    };
    // the MS1 spectra of the next run are read while a run is measured (reading is mostly waiting for the disk: ~1.8 of ~4.3 s per run)
    const auto load_run = [&in](Size r) {
      PeakFileOptions options = MzMLFile().getOptions();
      options.setMSLevels({1});
      PeakMap map;
      FLASHHelperClasses::loadExperiment(in[r], map, options, ProgressLogger::NONE);
      map.sortSpectra(true);
      return map;
    };
    // -out_id: per run, the PrSMs assigned to a consensus feature: (consensus, |mass error| ppm, PrSM)
    std::vector<std::vector<std::tuple<Size, double, PrSM>>> assigned(run_count);
    std::future<PeakMap> next_map = std::async(std::launch::async, load_run, (Size)0);
    for (Size r = 0; r < run_count; ++r)
    {
      PeakMap map = next_map.get();
      if (r + 1 < run_count) next_map = std::async(std::launch::async, load_run, r + 1);
      std::vector<FLASHQuantAlgorithm::QuantTarget> targets;
      std::vector<std::pair<Size, bool>> owner; // (consensus, decoy)
      for (Size k = 0; k < kept.size(); ++k)
      {
        const auto& s = summary[k];
        const double center = (s.aligned + (anchors[r].empty() ? 0.0 : offset_(anchors[r], s.aligned))) * 60.0, half = s.half * 60.0;
        targets.push_back({s.mass, center - half, center + half, s.zmin, s.zmax});
        owner.emplace_back(k, false);
        if (kept[k].runs.count(r) == 0 && mbr_fdr >= 0)
        {
          targets.push_back({s.mass - decoy_shift[k], center - half, center + half, s.zmin, s.zmax});
          owner.emplace_back(k, true);
        }
      }
      FLASHQuantAlgorithm fq;
      std::vector<FLASHQuantAlgorithm::MassAbundance> abundances;
      fq.quantifyTargets(map, targets, avg, abundances);
      std::vector<FDEvidence> evidence(targets.size());
      if (mbr_fdr >= 0)
      {
        std::vector<double> rts(map.size());
        std::vector<float> limits(map.size());
        for (Size s = 0; s < map.size(); ++s)
        {
          rts[s] = map[s].getRT();
          limits[s] = detectionLimit_(map[s]);
        }
#pragma omp parallel for schedule(dynamic, 16)
        for (SignedSize t = 0; t < (SignedSize)targets.size(); ++t) evidence[t] = fdEvidence_(map, rts, limits, targets[t], avg, fd_tol, fd_min_cos, fd_min_snr);
      }
      std::vector<std::pair<double, double>> ids;
      if (!in_id.empty() && mbr_fdr >= 0)
      {
        Size without_q = 0;
        ids = readIdentifications_(in_id[r], id_fdr, without_q);
        if (without_q > 0)
        {
          OPENMS_LOG_WARN << "Warning: " << in_id[r] << ": " << without_q << " PrSMs without a q-value (FLASHTnT without decoys) are not used" << endl;
        }
      }
      for (Size t = 0; t < targets.size(); ++t)
      {
        const auto [k, decoy] = owner[t];
        if (decoy) decoy_value[k][r] = abundances[t].first;
        else value[k][r] = abundances[t];
        fd_evidence[k][r][decoy ? 1 : 0] = evidence[t];
        if (ids.empty()) continue;
        const auto& target = targets[t];
        id_error[k][r][decoy ? 1 : 0] = id_match(ids, target.mass, target.rt_start, target.rt_end);
      }
      if (!out_id.empty())
      {
        // the features' windows in this run, by mass: each PrSM goes to the feature of the smallest |error| (nearest isotope) whose
        // window holds it
        std::vector<std::pair<double, Size>> by_mass; // (mass, target)
        for (Size t = 0; t < targets.size(); ++t)
        {
          if (!owner[t].second) by_mass.emplace_back(targets[t].mass, t);
        }
        std::sort(by_mass.begin(), by_mass.end());
        for (const PrSM& p : readPrSMs_(in_id[r], id_fdr))
        {
          double best = -1;
          Size best_t = 0;
          for (int i = -1; i <= 1; ++i)
          {
            const double m = p.mass - i * Constants::ISOTOPE_MASSDIFF_55K_U; // the feature's mass if the PrSM is i isotopes off
            const double tol = m * fd_tol;
            for (auto it = std::lower_bound(by_mass.begin(), by_mass.end(), std::make_pair(m - tol, (Size)0));
                 it != by_mass.end() && it->first <= m + tol; ++it)
            {
              const auto& target = targets[it->second];
              if (p.rt < target.rt_start || p.rt > target.rt_end) continue;
              const double e = std::abs(m - target.mass) / target.mass * 1e6;
              if (best < 0 || e < best)
              {
                best = e;
                best_t = it->second;
              }
            }
          }
          if (best >= 0) assigned[r].emplace_back(owner[best_t].first, best, p);
        }
      }
      OPENMS_LOG_INFO << "Run " << (r + 1) << " of " << run_count << ": " << targets.size() << " targets measured" << endl;
    }

    // The score of a transfer: FLASHDeconv's Qscore2D of its mass in its window, if it has a raw abundance and FLASHDeconv's mass of it
    // recurs within a gate of the median mass error of the features where detected: in 3 or more scans within ten times that median
    // (at most FLASHDeconv's tolerance), in 2 scans within three times; in 1 scan only for a mass of 8 kDa or more with a Qscore of at least
    // 0.5, within three times; else 0 (kept only with -mbr_fdr 1, or as mbr_id). Its decoy is scored the same way (by the feature's mass).
    // Recurrence, not the mass error, separates the transfers from chance (PXD042169 SpikeIn + onlyEcoli, decoys vs transfers within
    // three times the median: 1 scan 66 vs 1330, 2 scans 8 vs 2109, 3 or more 4 vs 17978): with one gate of three times the median for
    // all, the chance matches were nearly all of 1 scan, and transfers found scan after scan but a little further off were lost (the
    // weaker transfers' median |error| is 0.33 ppm, the detected features' 0.13 ppm; FLASHDeconv itself traces masses within 10 ppm). A
    // single scan of a large species carries many peaks: no decoy of 8 kDa or more passed in 1 scan with a Qscore of at least 0.5.
    // PXD042169, by this rule (one gate of three times before): transfers kept SpikeIn 19601 -> 23420, SpikeIn + onlyEcoli 21417 ->
    // 25624, ProteomeMix 13739 -> 17218 with fewer nulls (73 / 77 / 78 -> 24 / 19 / 42); E. coli into the onlyEcoli runs 42.6 -> 52.5%,
    // ProtMix into them 2 -> 1 of 99.
    // FLASHDeconv's mass, from all the peaks of the envelope in a scan, is far more precise than any single trace. The error is taken
    // against the feature's own mass measured the same way where detected (the median signed error of those runs; for a detected run, of
    // the other ones): the mass of the FLASHDeconv feature table is another estimate, off by a fixed amount for a quarter of the features.
    // PXD042169: median |error| 0.12 ppm where detected, 2-3 ppm for the decoys.
    std::vector<double> reference(kept.size(), 0); // per feature: the median signed error (ppm) of its mass where detected
    double mass_gate = 0, wide_gate = 0; // three and ten times the median error (at most FLASHDeconv's tolerance)
    {
      std::vector<double> errors;
      for (Size k = 0; k < kept.size(); ++k)
      {
        std::vector<std::pair<Size, double>> own; // (run, signed error) where detected
        for (const auto& [r, i] : kept[k].runs)
        {
          if (fd_evidence[k][r][0].scans > 0) own.emplace_back(r, fd_evidence[k][r][0].ppm);
        }
        std::vector<double> all;
        for (const auto& [r, e] : own) all.push_back(e);
        reference[k] = centreMedian_(all);
        for (const auto& [r, e] : own)
        {
          std::vector<double> others;
          for (const auto& [r2, e2] : own) { if (r2 != r) others.push_back(e2); }
          if (!others.empty()) errors.push_back(std::abs(e - centreMedian_(others)));
        }
      }
      mass_gate = 3 * median_(errors);
      wide_gate = std::min(10 * median_(errors), fd_tol * 1e6);
      // no feature with FLASHDeconv's evidence in two detected runs: no spread to take the gate from; FLASHDeconv's tolerance instead
      // (the decoys pass the same gate, so the FDR still holds; a gate of 0 rejected every transfer)
      const bool any_transfer = std::any_of(kept.begin(), kept.end(), [run_count](const Consensus& c) { return c.runs.size() < run_count; });
      if (errors.empty() && mbr_fdr >= 0 && any_transfer)
      {
        mass_gate = wide_gate = fd_tol * 1e6;
        OPENMS_LOG_WARN << "Warning: no feature has FLASHDeconv's evidence in two runs where it was detected; the |mass error| gate of the match between runs "
                           "is FLASHDeconv's tolerance (" << mass_gate << " ppm)" << endl;
      }
      if (mbr_fdr >= 0)
      {
        OPENMS_LOG_INFO << "Match between runs: |mass error| gates from the detected features: " << mass_gate << " ppm (1-2 scans), " << wide_gate
                        << " ppm (3 or more scans)" << endl;
      }
    }
    // -in_id: the gate of the identifications' mass errors, three times the median error of those of the features where detected (true
    // identifications of the mass, nearly all: PXD042169 ProteomeMix median 0.34 ppm, 95% within 1.5 ppm; at FLASHDeconv's tolerance
    // of 10 ppm, chance matches of a mass were ~2.6% of the identified transfers); FLASHDeconv's tolerance without any
    if (!in_id.empty() && mbr_fdr >= 0)
    {
      std::vector<double> detected_errors;
      for (Size k = 0; k < kept.size(); ++k)
      {
        for (const auto& [r, i] : kept[k].runs)
        {
          if (id_error[k][r][0] >= 0) detected_errors.push_back(id_error[k][r][0]);
        }
      }
      const double id_gate = detected_errors.empty() ? fd_tol * 1e6 : 3 * median_(detected_errors);
      for (Size k = 0; k < kept.size(); ++k)
      {
        for (Size r = 0; r < run_count; ++r)
        {
          if (kept[k].runs.count(r) > 0) continue;
          for (int kind = 0; kind < 2; ++kind) identified[k][r][kind] = id_error[k][r][kind] >= 0 && id_error[k][r][kind] <= id_gate;
        }
      }
      OPENMS_LOG_INFO << "Match between runs: identifications (-in_id): |mass error| gate " << id_gate << " ppm (from " << detected_errors.size()
                      << " identified detections)" << endl;
    }
    // the score of a transfer: FLASHDeconv's evidence alone (see above)
    const auto score_of = [&](const FDEvidence& e, double v, double ref, double mass) {
      if (!(v > 0) || e.scans == 0) return 0.0;
      const double error = std::abs(e.ppm - ref);
      const bool passes = e.scans >= 3 ? error <= wide_gate : e.scans == 2 ? error <= mass_gate : (mass >= 8000 && e.qscore2d >= 0.5 && error <= mass_gate);
      return passes ? e.qscore2d : 0.0;
    };
    // -in_id: the identified transfers (with a raw abundance) are a class with an FDR of its own: the decoy masses the runs identify (of
    // all the transfers, with or without a raw abundance: the stricter count) over the identified transfers. Within -mbr_fdr, they are
    // kept (mbr_id) and leave the targets of the other transfers; else they stay among them. The other transfers are scored against
    // all the decoys by FLASHDeconv's evidence, as without identifications (one null per transfer, so a few more nulls than targets: the
    // identifications never let another transfer through). Pooled with the others, the identified transfers (2-3% of the transfers of
    // PXD042169, nearly no nulls) would loosen the threshold of the others.
    const auto is_identified = [&](Size k, Size r) { return value[k][r].first > 0 && identified[k][r][0]; };
    std::vector<std::vector<double>> score(kept.size(), std::vector<double>(run_count, 0));
    std::vector<double> target_scores, null_scores;
    Size id_targets = 0, id_nulls = 0;
    for (Size k = 0; k < kept.size() && mbr_fdr >= 0; ++k)
    {
      for (Size r = 0; r < run_count; ++r)
      {
        if (kept[k].runs.count(r) > 0) continue;
        score[k][r] = score_of(fd_evidence[k][r][0], value[k][r].first, reference[k], summary[k].mass);
        null_scores.push_back(score_of(fd_evidence[k][r][1], decoy_value[k][r], reference[k], summary[k].mass));
        id_targets += is_identified(k, r);
        id_nulls += identified[k][r][1];
      }
    }
    const double id_q = id_targets > 0 ? std::min(1.0, (double)id_nulls / (double)id_targets) : 1.0;
    const bool id_class = id_targets > 0 && id_q <= mbr_fdr;
    for (Size k = 0; k < kept.size() && mbr_fdr >= 0; ++k)
    {
      for (Size r = 0; r < run_count; ++r)
      {
        if (kept[k].runs.count(r) > 0 || (id_class && is_identified(k, r))) continue;
        target_scores.push_back(score[k][r]);
      }
    }
    if (!in_id.empty() && mbr_fdr >= 0)
    {
      OPENMS_LOG_INFO << "Match between runs: identified transfers " << id_targets << ", identified decoy masses " << id_nulls
                      << (id_targets == 0 ? std::string("")
                          : std::string(": FDR ") + StringUtils::toStr(id_q) + (id_class ? " (kept: mbr_id)" : " (above -mbr_fdr: judged by FLASHDeconv's evidence like the others)"))
                      << endl;
    }
    // q-values: the FDR of a threshold T is (null scores >= T) / (target scores >= T) (one null per transfer); the q-value of a
    // score is the lowest FDR of the thresholds at or below it
    std::sort(target_scores.begin(), target_scores.end());
    std::sort(null_scores.begin(), null_scores.end());
    std::map<double, double> qvalue; // distinct target score -> q-value
    {
      double lowest = std::numeric_limits<double>::max();
      for (Size i = 0; i < target_scores.size(); ++i)
      {
        if (i > 0 && target_scores[i] == target_scores[i - 1]) continue;
        const double t = target_scores[i];
        const double targets_above = (double)(target_scores.size() - i);
        const double nulls_above = (double)(null_scores.end() - std::lower_bound(null_scores.begin(), null_scores.end(), t));
        lowest = std::min(lowest, std::min(1.0, nulls_above / targets_above));
        qvalue[t] = lowest;
      }
    }
    const auto q_of = [&qvalue](double t) { const auto it = qvalue.find(t); return it == qvalue.end() ? 1.0 : it->second; };
    double threshold = std::numeric_limits<double>::infinity();
    for (const auto& [t, q] : qvalue)
    {
      if (t > 0 && q <= mbr_fdr) { threshold = t; break; }
    }
    if (mbr_fdr >= 0)
    {
      const Size null_above = null_scores.end() - std::lower_bound(null_scores.begin(), null_scores.end(), threshold);
      const Size target_above = target_scores.end() - std::lower_bound(target_scores.begin(), target_scores.end(), threshold);
      OPENMS_LOG_INFO << "Match between runs: " << target_scores.size() << " transfers" << (id_class ? " not identified" : "") << "; score threshold at FDR "
                      << mbr_fdr << ": " << threshold
                      << " (targets at or above it " << target_above << ", nulls " << null_above << ")" << endl;

    }

    const std::string detail_out = getStringOption_("out_mbr_detail");
    if (!detail_out.empty())
    {
      std::ofstream ds(detail_out);
      ds << "ConsensusIndex\tRun\tStatus\tValue\tQscore2D\tMaxQscore\tScans\tMassErrorPPM\tDecoyValue\tDecoyQscore2D\tDecoyMaxQscore\tDecoyScans\tDecoyMassErrorPPM"
            "\tReferenceErrorPPM\tDetectedMedian" << (in_id.empty() ? "" : "\tIdentified\tDecoyIdentified") << "\n";
      for (Size k = 0; k < kept.size(); ++k)
      {
        std::vector<double> det;
        for (const auto& [r, i] : kept[k].runs)
        {
          if (value[k][r].first > 0) det.push_back(value[k][r].first);
        }
        const double det_median = median_(det);
        for (Size r = 0; r < run_count; ++r)
        {
          const bool detected = kept[k].runs.count(r) > 0;
          ds << k << "\t" << r << "\t" << (detected ? "detected" : "transfer") << "\t" << std::to_string(value[k][r].first);
          for (int kind = 0; kind < 2; ++kind)
          {
            // empty where not computed: no evidence with -mbr_fdr < 0, no decoy for a detected run
            const bool computed = mbr_fdr >= 0 && (kind == 0 || !detected);
            const auto& e = fd_evidence[k][r][kind];
            if (kind == 1) ds << "\t" << (computed ? std::to_string(decoy_value[k][r]) : std::string(""));
            if (!computed)
            {
              ds << "\t\t\t\t";
              continue;
            }
            ds << "\t" << std::to_string(e.qscore2d) << "\t" << std::to_string(e.max_qscore) << "\t" << e.scans << "\t" << std::to_string(e.ppm);
          }
          ds << "\t" << (mbr_fdr >= 0 ? std::to_string(reference[k]) : std::string("")) << "\t" << std::to_string(det_median);
          // with -in_id, for a transfer: 1 if its run identifies the mass (its decoy mass) in its window
          if (!in_id.empty())
          {
            const bool id_computed = mbr_fdr >= 0 && !detected;
            ds << "\t" << (id_computed ? std::string(identified[k][r][0] ? "1" : "0") : std::string("")) << "\t" << (id_computed ? std::string(identified[k][r][1] ? "1" : "0") : std::string(""));
          }
          ds << "\n";
        }
      }
    }

    // output
    std::ofstream os(out);
    os << "ConsensusIndex\tMonoisotopicMass\tApexRetentionTime\tMinCharge\tMaxCharge\tDetectedRuns";
    // the run names (file base names); a name two inputs share gets the run's position appended (_1, _2, ...), again while it equals
    // another run's name (run, run, run_2: run_1, run_2_2, run_2), so that no column repeats
    std::vector<std::string> names;
    for (const auto& f : in) names.push_back(completeBaseName_(f));
    {
      std::vector<bool> renamed(run_count, false);
      for (bool repeated = true; repeated;)
      {
        repeated = false;
        std::map<std::string, Size> name_count;
        for (const auto& n : names) name_count[n]++;
        // first round: every run of a shared name; later: only the renamed ones (the names of the files stay)
        const bool first_round = std::none_of(renamed.begin(), renamed.end(), [](bool b) { return b; });
        for (Size r = 0; r < run_count; ++r)
        {
          if (name_count[names[r]] < 2) continue;
          if (!renamed[r] && !first_round) continue;
          names[r] += "_" + StringUtils::toStr(r + 1);
          renamed[r] = true;
          repeated = true;
        }
      }
    }
    if (!out_id.empty())
    {
      std::ofstream is(out_id);
      is << "ConsensusIndex\tRun\tPrSMIndex\tScan\tRetentionTime\tProteinAccession\tProteinDescription\tProteoformIndex\tProforma\tProteoformMass"
            "\tPrecursorMass\tProteoformLevelQvalue\tMassErrorPPM\n";
      std::vector<std::tuple<Size, Size, Size>> order; // (consensus, run, position): by feature, then run, in the PrSM file's order
      for (Size r = 0; r < run_count; ++r)
      {
        for (Size i = 0; i < assigned[r].size(); ++i) order.emplace_back(std::get<0>(assigned[r][i]), r, i);
      }
      std::sort(order.begin(), order.end());
      Size identified_features = 0;
      for (Size i = 0; i < order.size(); ++i)
      {
        const auto [k, r, pos] = order[i];
        if (i == 0 || std::get<0>(order[i - 1]) != k) ++identified_features;
        const double error = std::get<1>(assigned[r][pos]);
        const PrSM& p = std::get<2>(assigned[r][pos]);
        is << k << "\t" << names[r];
        for (const auto& f : p.fields) is << "\t" << f;
        is << "\t" << std::to_string(error) << "\n";
      }
      OPENMS_LOG_INFO << "Identifications (-out_id): " << order.size() << " PrSMs assigned to " << identified_features << " of " << kept.size()
                      << " consensus features" << endl;
    }
    for (const auto& n : names) os << "\tAbundance_" << n << "\tStatus_" << n << "\tMatchScore_" << n << "\tMatchQvalue_" << n;
    os << "\n";
    Size mbr_kept = 0, mbr_total = 0, mbr_id_kept = 0;
    for (Size k = 0; k < kept.size(); ++k)
    {
      const auto& s = summary[k];
      os << k << "\t" << std::to_string(s.mass) << "\t" << s.aligned << "\t" << s.zmin << "\t" << s.zmax << "\t" << kept[k].runs.size();
      for (Size r = 0; r < run_count; ++r)
      {
        const double v = value[k][r].first;
        if (kept[k].runs.count(r) > 0)
        {
          os << "\t" << (v > 0 ? std::to_string(v) : std::string("nan")) << "\tdetected\t\t";
          continue;
        }
        mbr_total++;
        const bool id = id_class && is_identified(k, r);
        const double q = (mbr_fdr >= 0 && v > 0) ? (id ? id_q : q_of(score[k][r])) : 1.0;
        const bool ok = mbr_fdr >= 0 && v > 0 && q <= mbr_fdr && (id || score[k][r] >= mbr_min_score);
        mbr_kept += ok;
        mbr_id_kept += ok && id;
        os << "\t" << (ok ? std::to_string(v) : std::string("nan")) << "\t" << (ok ? (id ? "mbr_id" : "mbr") : "nan") << "\t" << (mbr_fdr >= 0 ? std::to_string(score[k][r]) : std::string(""))
           << "\t" << (mbr_fdr >= 0 ? std::to_string(q) : std::string(""));
      }
      os << "\n";
    }
    OPENMS_LOG_INFO << "Match between runs: " << mbr_kept << " of " << mbr_total << " values kept (q-value <= " << mbr_fdr
                    << (mbr_min_score > 0 ? std::string(", score >= ") + StringUtils::toStr(mbr_min_score) + " unless mbr_id" : std::string("")) << ")"
                    << (in_id.empty() ? std::string("") : std::string(", identified in their run (mbr_id) ") + StringUtils::toStr(mbr_id_kept)) << endl;
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPFLASHQuantConsensus tool;
  return tool.main(argc, argv);
}

/// @endcond
