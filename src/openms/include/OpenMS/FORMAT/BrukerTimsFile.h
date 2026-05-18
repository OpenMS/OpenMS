// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <cstdint>
#include <string>

class TimsDataHandle;

namespace OpenMS
{
  class FullSwathFileConsumer;

  /**
   * @brief Reader for Bruker TimsTOF .d directories via opentims.
   *
   * Supports DDA-PASEF, DIA-PASEF, and raw frame-level 4D access.
   * Ion mobility data is stored in VSSC (1/K0) units using IM_PEAK format
   * for MS1 and DIA MS2, and scalar drift times for DDA MS2.
   *
   * In FRAME export mode, raw TOF indices and intensities are returned without
   * any signal processing. TOF-to-m/z and scan-to-IM conversions are applied,
   * but no smoothing or centroiding is performed.
   *
   * @section export_modes Export modes
   *
   * - **AUTO** (default): detects DDA vs DIA by querying for SWATH windows.
   *   DDA uses per-precursor MS2 spectra with scalar drift times; DIA uses
   *   per-frame MS2 spectra split by SWATH isolation window with per-peak IM.
   * - **SPECTRUM**: forces the DDA (per-precursor) path regardless of
   *   acquisition type.
   * - **FRAME**: returns raw 4D frames as IM_PEAK spectra (per-peak IM
   *   in FloatDataArray) for both MS1 and MS2, without precursor grouping or
   *   SWATH splitting. MS2 frame-level access may not be available for all
   *   acquisition types.
   */
  class OPENMS_DLLAPI BrukerTimsFile : public ProgressLogger
  {
  public:
    /// Processing and export configuration
    struct Config
    {
      double calibration_tolerance = 0.0;  ///< m/z recalibration tolerance in Da (0 = default 0.1 Da)
      bool calibrate = false;              ///< Enable m/z recalibration (off by default; may fail on some datasets)

      bool load_ms1 = true;                ///< Load MS1 spectra. Disable (false) for MS2-only workflows
                                           ///< (peptide database search) where MS1 surveys are not needed —
                                           ///< cuts memory and time substantially. Affects all export modes.

      /// Centroiding algorithms for IM-PASEF frames.
      ///   - Off:           no IM-axis centroiding (raw detector peaks pass through).
      ///   - Greedy2D:      2D greedy clustering in (m/z, IM) using
      ///                    ms1_centroid_mz_ppm + ms1_centroid_im_pct (the legacy path,
      ///                    adapted from Sage; Lazear 2023 doi:10.1021/acs.jproteome.3c00486).
      ///   - HillBased:     IM-axis hill detection (peaks linked across consecutive
      ///                    IM scans within a frame, with valley splitting); modeled on
      ///                    Biosaur2's hill detection but applied to the IM axis. Uses
      ///                    ms1_centroid_mz_ppm as the m/z tolerance and
      ///                    centroid_valley_factor / centroid_min_hill_length as the
      ///                    hill-specific knobs.
      ///   - PeakPickerIM:  mobilogram-based picking via OpenMS::PeakPickerIM
      ///                    (pickIMTraces). Builds an MSSpectrum from the frame's
      ///                    (m/z, intensity, IM) triples and runs the mobilogram
      ///                    extraction + per-trace centroiding pipeline. Parameters
      ///                    are read from peak_picker_im_params (subsection
      ///                    pickIMTraces:); the other centroiding-tolerance fields
      ///                    above are ignored. Also honors expose_hill_bounds —
      ///                    emits FWHM-derived bounding boxes under the same
      ///                    FloatDataArray names as HillBased. Not supported on
      ///                    DDA-MS2 (would warn and fall back to Off).
      enum class CentroidAlgo { OFF, GREEDY2D, HILL_BASED, PEAK_PICKER_IM };

      /// MS1 centroiding algorithm. Default Off. If left at Off but the legacy
      /// fields ms1_centroid_mz_ppm and ms1_centroid_im_pct are both > 0, the
      /// loader treats that as an implicit Greedy2D selection (back-compat).
      CentroidAlgo ms1_centroid_algo = CentroidAlgo::OFF;

      /// MS2 centroiding algorithm for both DIA-PASEF and DDA-PASEF.
      ///
      /// DIA-MS2: Off / Greedy2D / HillBased / PeakPickerIM. Greedy2D maps to the
      /// legacy 2D Gaussian + local-maxima path (equivalent to dia_ms2_centroid=true).
      /// HillBased and PeakPickerIM both run inside per-frame and aggregated
      /// (dia_ms2_n_neighbors > 0) DIA paths. Takes precedence over the legacy
      /// dia_ms2_centroid boolean.
      ///
      /// DDA-MS2: Off / HillBased. HillBased replaces the default TOF-domain
      /// processing pipeline (sort/group/smooth/local-max in TOF space) with
      /// IM-axis hill detection on per-precursor (m/z, intensity, IM) triples.
      /// Output schema is unchanged — DDA-MS2 still emits a single scalar
      /// drift_time per spectrum (no per-peak IM array). 'Greedy2D' on DDA-MS2
      /// falls through to the default TOF pipeline (DDA has no Greedy2D-style
      /// box clustering equivalent). 'PeakPickerIM' is not supported on DDA-MS2
      /// (its mobilogram pipeline expects per-peak IM arrays); the resolver
      /// warns once and falls back to Off.
      CentroidAlgo ms2_centroid_algo = CentroidAlgo::OFF;

      float ms1_centroid_mz_ppm = 10.0f; ///< MS1 centroiding m/z tolerance in ppm (used by Greedy2D and HillBased). Default 10 ppm is tuned for HillBased TIMS-PASEF MS1: detector-centroided peaks drift up to ~10 ppm in m/z between IM scans of the same ion, so 5 ppm under-links real multi-scan hills (only ~4% of MS1 hills end up multi-scan at 5 ppm vs ~9% at 10 ppm on bundled DIA fixture). Greedy2D still requires ms1_centroid_im_pct > 0 in addition. Set to 0 to disable MS1 centroiding outright.
      float ms1_centroid_im_pct = 0.0f;  ///< Greedy2D-only: MS1 ion-mobility tolerance in percent (0 = disabled, suggested: 3.0). Ignored by HillBased.
      int   ms1_centroid_max_peaks = 100000; ///< Greedy2D-only: upper bound on centroided peaks per MS1 spectrum. Caps retention at the top-intensity peaks; low-intensity tail is dropped when the limit is hit (a warning is emitted at that point). Sized for aggregated MS1 where frame-stacked grids commonly exceed 200k peaks; per-frame MS1 rarely hits even 10k.

      float  ms2_centroid_mz_ppm        = 20.0f; ///< HillBased MS2 m/z tolerance in ppm. 20 ppm is the DIA-PASEF-tuned default. Set to 0 to refuse to run HillBased MS2 centroiding (the algo helper falls back to Off and warns).
      double centroid_valley_factor     = 1.3;  ///< HillBased: hill valley factor (hvf). Both (left_max/valley) and (right_max/valley) must exceed this to split a hill at a valley. Smaller = more aggressive splitting. Shared between MS1 and MS2.
      Size   ms1_centroid_min_hill_length = 1;  ///< HillBased MS1: minimum number of IM scans a hill must span. Default 1 keeps single-IM-scan ions (common on detector-centroided TIMS-PASEF MS1: ~75% of peaks have no same-m/z partner in the previous IM scan within 100 ppm).
      Size   ms2_centroid_min_hill_length = 2;  ///< HillBased MS2: minimum number of IM scans a hill must span. Default 2 is the DIA-PASEF-tuned value (rejects single-scan singletons; ~33% of unfiltered hill output survives, bringing volume close to the legacy Gaussian-smooth + local-maxima path). DDA-PASEF users should override to 1: precursors span a narrow IM range so most fragments are seen in only one IM scan, and min=2 drops ~93% of peaks there.
      Size   centroid_max_scan_gap        = 0;  ///< HillBased (MS1 + MS2): maximum number of consecutive empty IM scans a hill may bridge while linking. 0 = strict consecutive-scan linking (Biosaur2 default). 1 = a single empty scan at the hill's m/z is tolerated. Useful on detector-centroided TIMS-PASEF where an ion may occasionally fail to register in one IM scan. Note: hill length counts only the scans where the ion was actually observed, not the bridged gap.
      bool   expose_hill_bounds           = false; ///< HillBased + PeakPickerIM (MS1 + DIA-MS2): when true, attach four extra FloatDataArrays per centroided spectrum ("im lower bound", "im upper bound", "m/z lower bound", "m/z upper bound") giving each centroid's bounding box. HillBased emits the source-peak extent (min/max IM and m/z of all peaks that contributed to the hill); PeakPickerIM emits the FWHM-derived support window (centroid_im ± IM_FWHM/2 and the spline m/z FWHM endpoints). Same array names so visualization tools (e.g. tools/scripts/plot_pasef_frames.py --show-hill-bounds) work uniformly. Bloats the mzML by roughly +25% on centroided spectra. Has no effect on DDA-MS2 hill (scalar drift_time schema).

      /// PeakPickerIM (MS1 + DIA-MS2): parameter overrides forwarded to
      /// PeakPickerIM. Treated as a sparse overlay on PeakPickerIM's own
      /// defaults: empty (the default) means "use all defaults"; any keys
      /// present here are merged in via Param::update with type/restriction
      /// validation. Only the `pickIMTraces:*` subsection is consumed
      /// (pickIMCluster / pickIMElutionProfiles are not currently wired into
      /// the loader). Set fields directly via
      /// peak_picker_im_params.setValue("pickIMTraces:gauss_ppm_tolerance", 5.0)
      /// etc. to override defaults.
      Param  peak_picker_im_params;

      /// MS1 + DIA-MS2 isotopic-partner prefilter (off by default).
      /// Applied after aggregation (or after raw extraction when aggregation
      /// is disabled), before the centroider dispatch. Drops peaks that lack
      /// at least one isotopic partner at m/z ± C13C12_MASSDIFF_U / q (q in
      /// {1, 2, 3, 4, 5}) within ±isotopic_prefilter_tol_da AND |Δscan_id|
      /// <= 1. Cleans up isolated detector-noise singletons; preserves both
      /// the monoisotopic peak and the isotopologue (mutual evidence).
      /// Pure existence check — no intensity-ratio or averagine model.
      /// Not applied to DDA-MS2 (which has no per-peak IM array).
      bool   isotopic_prefilter        = false;
      double isotopic_prefilter_tol_da = 0.05;   ///< Da tolerance for isotopic-partner matching (broad by design — survives per-scan calibration jitter).

      int dia_ms2_n_neighbors = 0;  ///< DIA MS2 frame aggregation: number of adjacent frames on each side (0 = disabled / raw per-frame export, 1 = 3-frame sum, 2 = 5-frame sum). Kept at 0 by default so the raw DIA-MS2 export path is the out-of-the-box behavior — this knob switches the entire DIA-MS2 export pipeline (sum + denoise) regardless of centroiding algo, so flipping it changes output for every caller. Set to 2 alongside ms2_centroid_algo=hillbased + min_hill_length=2 + mz_ppm=20 for the DIA-PASEF hill recipe; set to 1/2 with algo=off for raw RT-summed export.
      int dia_ms2_min_support = 1;  ///< DIA MS2 denoising: minimum occupied neighbors in 3x3 (mz x IM) grid to retain a point (center excluded)
      bool dia_ms2_centroid = false; ///< DIA MS2 2D peak picking: apply Gaussian smoothing + local maxima detection to produce IM_CENTROIDED spectra

      int    ms1_n_neighbors         = 0;    ///< MS1 frame aggregation: adjacent MS1 frames on each side
                                             ///< (0 = disabled, 1 = 3-frame sum, 2 = 5-frame sum). Applies
                                             ///< to both DIA and DDA; ignored in FRAME export mode.
      int    ms1_min_support         = 0;    ///< MS1 denoising after aggregation: min occupied 3x3 neighbors
                                             ///< (0 = disabled). Only effective when ms1_n_neighbors > 0.
                                             ///< Differs from dia_ms2_min_support (default 1) to preserve
                                             ///< zero-change defaults for the MS1 path.
      double ms1_max_rt_distance_sec = 0.0;  ///< Advanced: cap the RT distance (seconds) between a neighbor
                                             ///< MS1 frame and the center frame during aggregation
                                             ///< (0.0 = no cap). The center frame is always included
                                             ///< regardless of this cap. Recommended for DDA where MS1
                                             ///< frame cadence is irregular.

      enum ExportMode { AUTO, SPECTRUM, FRAME };
      ExportMode export_mode = AUTO;       ///< AUTO detects DDA vs DIA; SPECTRUM forces per-precursor; FRAME returns raw 4D frames

      /// Strategy for converting TIMS scan indices to 1/K0 values.
      /// AUTO (default): tries Bruker SDK → rational (TimsCalibration table) → linear.
      enum class TimsCalibrationStrategy { AUTO, BRUKER_SDK, RATIONAL, LINEAR };
      TimsCalibrationStrategy tims_calibration_strategy = TimsCalibrationStrategy::AUTO;

      /// Pressure compensation strategy (only effective with BRUKER_SDK calibration).
      /// Corrects for ambient gas pressure drift during acquisition.
      /// Ignored when using RATIONAL or LINEAR strategies; a warning is logged.
      enum class PressureCompensation { NONE, GLOBAL, PER_FRAME };
      PressureCompensation pressure_compensation = PressureCompensation::NONE;

      /// Path to Bruker SDK library (timsdata.dll / libtimsdata.so).
      /// Empty string (default): discover from OPENMS_BRUKER_SDK_PATH env var.
      std::string bruker_sdk_path;
    };

    /// Metadata for constructing a FullSwathFileConsumer for DIA streaming.
    struct DIAStreamingMetadata
    {
      std::vector<OpenSwath::SwathMap> boundaries;  ///< one SwathMap per DIAWindow (MS2 only)
      int nr_ms1_spectra = 0;                       ///< number of MS1 frames
      std::vector<int> nr_ms2_spectra;              ///< per-window spectrum counts (parallel to boundaries)
    };

    /// Read DIA SWATH boundaries and spectrum counts from a .d directory (SQL only, no peak data).
    /// Also populates exp_settings with source file metadata.
    DIAStreamingMetadata readDIAMetadata(const String& path, ExperimentalSettings& exp_settings);
    /// @overload with explicit configuration
    DIAStreamingMetadata readDIAMetadata(const String& path, ExperimentalSettings& exp_settings,
                                         const Config& config);

    /// Stream DIA spectra to a consumer one-at-a-time without accumulating.
    /// MS2 spectra are always raw (no aggregation/denoising/centroiding).
    /// MS1 centroiding respects the Config settings.
    void loadDIAStreaming(const String& path, FullSwathFileConsumer& consumer);
    /// @overload with explicit configuration
    void loadDIAStreaming(const String& path, FullSwathFileConsumer& consumer,
                          const Config& config);

    /// Load entire .d directory into MSExperiment
    void load(const String& path, MSExperiment& exp);
    /// @overload with explicit configuration
    void load(const String& path, MSExperiment& exp, const Config& config);

    /// Feed spectra from a .d directory to a consumer.
    /// @note Currently loads the full dataset into memory before streaming to the consumer.
    ///       A future optimization should iterate frame-by-frame for true constant-memory operation.
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer);
    /// @overload with explicit configuration
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config);

  private:
    /// Load DDA-PASEF data: MS1 frames (IM_PEAK) + MS2 spectra (scalar IM)
    void loadDDA_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);

    /// Load DIA-PASEF data: MS1 frames (IM_PEAK) + MS2 frames split by SWATH window
    void loadDIA_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);

    /// Load raw frames without precursor grouping or SWATH splitting
    void loadFrames_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);

    /// Detect DDA vs DIA by checking for SWATH windows
    bool isDIA_(const String& tdf_path) const;

    /// Populate SourceFile metadata from the .d path (no peak data read)
    void loadExperimentalSettings_(const String& path, ExperimentalSettings& settings);
  };

} // namespace OpenMS

#endif // WITH_OPENTIMS
