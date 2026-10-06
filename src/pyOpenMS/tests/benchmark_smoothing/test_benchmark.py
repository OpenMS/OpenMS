"""
Unit and regression tests for the chromatogram smoothing benchmark suite.
Verifies determinism, parameter validation, sub-scan interpolation, and peak picking.
"""

import unittest
import numpy as np
import pyopenms

try:
    from .synthetic_data import (
        SyntheticChromatogramGenerator,
        generate_gaussian_peak,
        generate_emg_peak,
        add_noise,
    )
    from .metrics import (
        evaluate_signal_fidelity,
        evaluate_noise_reduction,
        evaluate_peak_picking,
        evaluate_valley_to_peak_ratio,
        parabolic_apex_interpolation,
        interpolated_fwhm,
    )
    from .benchmark_engine import (
        BenchmarkEngine,
        apply_modified_sinc,
        apply_savitzky_golay,
        run_peak_picker_hires,
    )
    from .real_data import (
        DatasetInfo,
        PASS00779_INFO,
        MTBLS404_INFO,
        LOCAL_FIXTURES_INFO,
        RealChromatogramSelector,
        SelectionManifest,
        RawTraceProfile,
        estimate_mad_noise,
        profile_raw_chromatogram,
        load_chromatograms_from_file,
    )
    from .real_metrics import (
        estimate_flank_noise,
        compute_valley_to_peak_ratio_real,
        evaluate_single_real_trace,
        summarize_real_dataset_results,
        evaluate_qc_replicates,
        RealSmootherMetrics,
        RealTraceEvaluation,
        RealDatasetBenchmarkResult,
        QCReplicateConsistencyResult,
    )
except ImportError:
    from synthetic_data import (
        SyntheticChromatogramGenerator,
        generate_gaussian_peak,
        generate_emg_peak,
        add_noise,
    )
    from metrics import (
        evaluate_signal_fidelity,
        evaluate_noise_reduction,
        evaluate_peak_picking,
        evaluate_valley_to_peak_ratio,
        parabolic_apex_interpolation,
        interpolated_fwhm,
    )
    from benchmark_engine import (
        BenchmarkEngine,
        apply_modified_sinc,
        apply_savitzky_golay,
        run_peak_picker_hires,
    )
    from real_data import (
        DatasetInfo,
        PASS00779_INFO,
        MTBLS404_INFO,
        LOCAL_FIXTURES_INFO,
        RealChromatogramSelector,
        SelectionManifest,
        RawTraceProfile,
        estimate_mad_noise,
        profile_raw_chromatogram,
        load_chromatograms_from_file,
    )
    from real_metrics import (
        estimate_flank_noise,
        compute_valley_to_peak_ratio_real,
        evaluate_single_real_trace,
        summarize_real_dataset_results,
        evaluate_qc_replicates,
        RealSmootherMetrics,
        RealTraceEvaluation,
        RealDatasetBenchmarkResult,
        QCReplicateConsistencyResult,
    )



class TestSyntheticData(unittest.TestCase):

    def test_determinism(self):
        """Verify that identical random seeds produce byte-identical chromatograms."""
        gen1 = SyntheticChromatogramGenerator(seed=123)
        gen2 = SyntheticChromatogramGenerator(seed=123)

        ds1 = gen1.create_narrow_peak()
        ds2 = gen2.create_narrow_peak()

        np.testing.assert_array_equal(ds1.rt, ds2.rt)
        np.testing.assert_array_equal(ds1.intensity_true, ds2.intensity_true)
        np.testing.assert_array_equal(ds1.intensity_noisy, ds2.intensity_noisy)

    def test_sub_scan_phase_offset(self):
        """Verify phase offset shifts peak apex correctly."""
        gen = SyntheticChromatogramGenerator(seed=42)
        p0 = gen.create_narrow_peak(phase_offset=0.0)
        p5 = gen.create_narrow_peak(phase_offset=0.5)

        self.assertAlmostEqual(p5.peaks[0].true_apex_rt - p0.peaks[0].true_apex_rt, 0.5, delta=1e-5)

    def test_tuning_datasets_isolation(self):
        """Verify tuning datasets use separate seeds and configurations from test datasets."""
        gen = SyntheticChromatogramGenerator(seed=42)
        tuning = gen.generate_tuning_datasets()
        test = gen.generate_test_datasets()

        # Verify no overlapping references and distinct data
        self.assertNotEqual(len(tuning["tuning_narrow"].rt), len(test["narrow"].rt))
        tuning_rts = [p.true_apex_rt for p in tuning["tuning_composite"].peaks]
        test_rts = [p.true_apex_rt for p in test["composite"].peaks]
        self.assertNotEqual(tuning_rts, test_rts)

    def test_analytical_area_gaussian(self):
        """Verify analytical Gaussian area matches numerical integration."""
        rt = np.linspace(0, 100, 1000)
        profile, meta = generate_gaussian_peak(rt, center_rt=50.0, height=1000.0, fwhm=10.0)

        trapz_fn = getattr(np, "trapezoid", getattr(np, "trapz", None))
        num_area = float(trapz_fn(profile, rt))

        rel_diff = abs(num_area - meta.true_area) / meta.true_area
        self.assertLess(rel_diff, 0.001)

    def test_emg_tailing(self):
        """Verify EMG peak has valid apex and positive tailing."""
        rt = np.linspace(0, 100, 1000)
        profile, meta = generate_emg_peak(rt, center_rt=30.0, height=1000.0, fwhm=10.0, tau=5.0)

        self.assertEqual(meta.peak_type, "emg_tailing")
        self.assertGreater(meta.true_area, 0.0)
        self.assertGreaterEqual(meta.true_apex_rt, 30.0)


class TestMetrics(unittest.TestCase):

    def test_parabolic_apex_interpolation(self):
        """Verify parabolic apex interpolation finds true sub-scan apex."""
        # Apex placed at 5.25 (between grid points 5.0 and 6.0)
        rt = np.arange(0, 10, 1.0)
        true_apex = 5.25
        y = 1000.0 * np.exp(-0.5 * ((rt - true_apex) / 1.5) ** 2)

        est_rt, est_int = parabolic_apex_interpolation(rt, y)
        self.assertAlmostEqual(est_rt, true_apex, delta=0.08)
        self.assertAlmostEqual(est_int, 1000.0, delta=20.0)

    def test_flank_interpolated_fwhm(self):
        """Verify flank-interpolated FWHM achieves sub-percent accuracy."""
        rt = np.arange(0, 20, 0.5)
        sigma = 2.0
        true_fwhm = 2.354820045 * sigma  # ~ 4.7096
        y = 1000.0 * np.exp(-0.5 * ((rt - 10.0) / sigma) ** 2)

        fwhm = interpolated_fwhm(rt, y, apex_int=1000.0, baseline=0.0)
        self.assertAlmostEqual(fwhm, true_fwhm, delta=0.05)

    def test_perfect_fidelity(self):
        """Verify that identical true and smoothed signals yield 0 RMSE and r=1."""
        rt = np.linspace(0, 50, 100)
        y = np.exp(-0.5 * ((rt - 25.0) / 3.0) ** 2)

        rmse, mae, max_err, r, peaks = evaluate_signal_fidelity(rt, y, y, [])
        self.assertAlmostEqual(rmse, 0.0, places=6)
        self.assertAlmostEqual(mae, 0.0, places=6)
        self.assertAlmostEqual(r, 1.0, places=6)

    def test_peak_picking_evaluation(self):
        """Verify TP/FP/FN assignment in peak picking evaluation."""
        from synthetic_data import PeakMetadata
        true_peaks = [
            PeakMetadata("gauss", true_apex_rt=20.0, true_apex_intensity=1000.0, true_fwhm=4.0, true_area=4000.0, rt_start=14.0, rt_end=26.0),
            PeakMetadata("gauss", true_apex_rt=50.0, true_apex_intensity=800.0, true_fwhm=5.0, true_area=4000.0, rt_start=42.0, rt_end=58.0),
        ]
        detected = [
            (20.2, 980.0, 4.1),   # TP
            (80.0, 50.0, 2.0),    # FP
        ]
        pp_metrics = evaluate_peak_picking(detected, true_peaks, rt_tolerance_factor=0.5, sn_threshold=2.0)
        self.assertEqual(pp_metrics.true_positives, 1)
        self.assertEqual(pp_metrics.false_positives, 1)
        self.assertEqual(pp_metrics.false_negatives, 1)
        self.assertEqual(pp_metrics.precision, 0.5)
        self.assertEqual(pp_metrics.recall, 0.5)


class TestBenchmarkEngine(unittest.TestCase):

    def setUp(self):
        self.engine = BenchmarkEngine(seed=42, quick_mode=True)

    def test_apply_smoothers_on_chromatogram(self):
        """Verify both smoothers process OpenMS MSChromatogram without error."""
        chrom = pyopenms.MSChromatogram()
        rt = [float(x) for x in range(30)]
        intensity = [float(100.0 + x * (30 - x)) for x in range(30)]
        chrom.set_peaks((rt, intensity))

        # Modified Sinc
        rt_ms, int_ms, t_ms = apply_modified_sinc(chrom, degree=6, m=7, is_ms1=False)
        self.assertEqual(len(rt_ms), 30)
        self.assertGreater(t_ms, 0.0)

        # Savitzky-Golay
        rt_sg, int_sg, t_sg = apply_savitzky_golay(chrom, frame_length=11, polynomial_order=4)
        self.assertEqual(len(rt_sg), 30)
        self.assertGreater(t_sg, 0.0)

    def test_independent_tuning_workflow(self):
        """Verify tuning selects parameters and held-out evaluation executes."""
        best_ms, best_sg = self.engine.tune_parameters()
        self.assertIn("degree", best_ms)
        self.assertIn("frame_length", best_sg)

        held_out = self.engine.run_held_out_evaluation(best_ms, best_sg)
        self.assertIn("narrow", held_out)
        self.assertIn("composite", held_out)

    def test_peak_picker_hires_direct_sn2(self):
        """Verify PeakPickerHiRes picks chromatogram apex at realistic S/N=2.0 without double smoothing."""
        rt = np.arange(0, 30, 0.5)
        intensity = 50.0 + 1000.0 * np.exp(-0.5 * ((rt - 15.0) / 2.0) ** 2)
        detected = run_peak_picker_hires(rt, intensity, signal_to_noise=2.0)

        self.assertEqual(len(detected), 1)
        self.assertAlmostEqual(detected[0][0], 15.0, delta=0.5)

    def test_matched_bandwidth_consistency(self):
        """Verify savitzkyGolayBandwidth and bandwidthToM frequency matching."""
        # SG order 4, frame length 11 (m=5)
        bw = pyopenms.ModifiedSincSmoother.savitzkyGolayBandwidth(4, 5)
        self.assertGreater(bw, 0.0)
        self.assertLess(bw, 0.5)

        m_ms = pyopenms.ModifiedSincSmoother.bandwidthToM(False, 6, bw)
        self.assertGreaterEqual(m_ms, 5)

        matched = self.engine.run_matched_bandwidth_comparison()
        self.assertGreater(len(matched), 0)
        self.assertIn("equivalent_bandwidth", matched[0])

    def test_overlapping_doublet_total_area(self):
        """Verify that total doublet area matches sum of individual ground truth peak areas."""
        gen = SyntheticChromatogramGenerator(seed=42)
        syn = gen.create_overlapping_peaks()
        self.assertTrue(syn.is_overlapping_doublet)
        self.assertAlmostEqual(syn.doublet_total_area, syn.peaks[0].true_area + syn.peaks[1].true_area, places=4)

        # Check that evaluate_configuration evaluates total doublet area correctly
        res = self.engine.evaluate_configuration(syn, "ModifiedSincSmoother", {"degree": 6, "m": 7, "is_ms1": False})
        self.assertIsNotNone(res.doublet_total_area_error_pct)
        self.assertLess(abs(res.doublet_total_area_error_pct), 2.0)


class TestRealDataInfrastructure(unittest.TestCase):

    def setUp(self):
        self.engine = BenchmarkEngine(seed=42, quick_mode=True)

    def test_estimate_mad_noise(self):
        """Verify robust MAD noise estimation on Gaussian noise with linear baseline slope."""
        rng = np.random.default_rng(123)
        true_sigma = 15.0
        n_points = 500
        noise = rng.normal(0, true_sigma, n_points)
        linear_slope = np.linspace(100.0, 500.0, n_points)
        y = linear_slope + noise

        # MAD should estimate noise close to true_sigma despite steep slope
        sigma_est = estimate_mad_noise(y)
        self.assertAlmostEqual(sigma_est, true_sigma, delta=3.0)

    def test_estimate_flank_noise(self):
        """Verify flank noise isolates baseline and ignores central intense peak."""
        rt = np.linspace(0, 100, 200)
        rng = np.random.default_rng(456)
        noise = rng.normal(0, 5.0, 200)
        # Intense peak in center
        peak = 10000.0 * np.exp(-0.5 * ((rt - 50.0) / 3.0) ** 2)
        y = peak + noise

        flank_sigma = estimate_flank_noise(y, flank_fraction=0.20)
        self.assertAlmostEqual(flank_sigma, 5.0, delta=1.5)

    def test_raw_chromatogram_profiling(self):
        """Verify RawTraceProfile computes raw properties and assigns valid strata."""
        rt = np.linspace(0, 60, 120)
        intensity = 50.0 + 50000.0 * np.exp(-0.5 * ((rt - 30.0) / 4.0) ** 2)
        c = pyopenms.MSChromatogram()
        c.set_peaks((rt.tolist(), intensity.tolist()))
        c.setNativeID("test_trace_1")

        prof = profile_raw_chromatogram(c, trace_id="trace_1", source_file="test.mzML")
        self.assertTrue(prof.is_valid)
        self.assertEqual(prof.num_points, 120)
        self.assertAlmostEqual(prof.apex_rt, 30.0, delta=1.0)
        self.assertAlmostEqual(prof.apex_intensity, 50050.0, delta=100.0)
        self.assertEqual(prof.intensity_stratum, "medium")  # 50,000 is between 1e4 and 1e5
        self.assertIn(prof.width_stratum, ["medium", "broad"])

    def test_deterministic_selection_and_seed(self):
        """Verify RealChromatogramSelector produces identical selections for identical seeds."""
        chroms = []
        rt = np.linspace(0, 60, 60)
        for i in range(20):
            c = pyopenms.MSChromatogram()
            height = 1000.0 * (i + 1)
            sigma = 1.0 + (i % 5)
            y = 50.0 + height * np.exp(-0.5 * ((rt - 30.0) / sigma) ** 2)
            c.set_peaks((rt.tolist(), y.tolist()))
            c.setNativeID(f"native_id_{i}")
            chroms.append((c, f"trace_{i:02d}", "file_a.mzML"))

        selector_a = RealChromatogramSelector(target_sample_size=6, seed=42)
        selector_b = RealChromatogramSelector(target_sample_size=6, seed=42)
        selector_c = RealChromatogramSelector(target_sample_size=6, seed=999)

        _, manifest_a, _ = selector_a.select(chroms, dataset_accession="TEST")
        _, manifest_b, _ = selector_b.select(chroms, dataset_accession="TEST")
        _, manifest_c, _ = selector_c.select(chroms, dataset_accession="TEST")

        self.assertEqual(manifest_a.selected_trace_ids, manifest_b.selected_trace_ids)
        self.assertEqual(manifest_a.selection_seed, manifest_b.selection_seed)
        self.assertNotEqual(manifest_a.selected_trace_ids, manifest_c.selected_trace_ids)

    def test_manifest_serialization(self):
        """Verify SelectionManifest exports to dictionary and JSON cleanly."""
        manifest = SelectionManifest(
            dataset_accession="PASS00779",
            source_files=["run1.mzML", "run2.mzML"],
            selection_seed=42,
            target_sample_size=50,
            total_candidates_inspected=120,
            total_valid_candidates=100,
            total_selected=50,
            rejection_counts={"insufficient_points": 10, "low_snr": 10},
            stratum_counts_available={"low_narrow": 20, "high_broad": 80},
            stratum_counts_selected={"low_narrow": 15, "high_broad": 35},
            selected_trace_ids=["trace_0", "trace_1"],
            selection_criteria={"min_points": 20, "min_snr": 3.0},
        )
        d = manifest.to_dict()
        self.assertEqual(d["dataset_accession"], "PASS00779")
        self.assertEqual(len(d["selected_trace_ids"]), 2)
        self.assertEqual(d["total_selected"], 50)

    def test_real_trace_metrics_computation(self):
        """Verify evaluate_single_real_trace calculates empirical metrics without ground truth."""
        rt = np.linspace(0, 50, 100)
        y = 100.0 + 20000.0 * np.exp(-0.5 * ((rt - 25.0) / 3.0) ** 2)
        # Add slight synthetic noise to test noise reduction
        rng = np.random.default_rng(77)
        y_noisy = y + rng.normal(0, 30.0, 100)

        c = pyopenms.MSChromatogram()
        c.set_peaks((rt.tolist(), y_noisy.tolist()))
        c.setNativeID("trace_eval_1")

        ms = pyopenms.ModifiedSincSmoother()
        sg = pyopenms.SavitzkyGolayFilter()

        res = evaluate_single_real_trace(
            chrom=c,
            trace_id="eval_1",
            source_file="test.mzML",
            intensity_stratum="medium",
            width_stratum="medium",
            ms_params={"degree": 6, "m": 12, "is_ms1": False},
            sg_params={"frame_length": 15, "polynomial_order": 4},
            ms_smoother_instance=ms,
            sg_smoother_instance=sg,
        )

        self.assertAlmostEqual(res.raw_apex_rt, 25.0, delta=1.0)
        self.assertLess(abs(res.modified_sinc.apex_rt_shift), 0.5)
        self.assertLess(abs(res.savitzky_golay.apex_rt_shift), 0.5)
        self.assertGreater(res.modified_sinc.noise_reduction_pct, 10.0)
        self.assertGreater(res.savitzky_golay.noise_reduction_pct, 10.0)
        self.assertGreater(res.modified_sinc.peak_area_ratio, 0.9)
        self.assertGreater(res.savitzky_golay.peak_area_ratio, 0.9)

    def test_qc_replicate_consistency_calculation(self):
        """Verify evaluate_qc_replicates computes CV% across replicate traces."""
        rt = np.linspace(0, 50, 100)
        reps = []
        rng = np.random.default_rng(99)
        for i in range(4):
            # Jitter apex RT by +/- 0.05 s and intensity by +/- 2%
            rt_shift = float(rng.normal(0, 0.05))
            int_scale = float(1.0 + rng.normal(0, 0.02))
            y = 50.0 + 10000.0 * int_scale * np.exp(-0.5 * ((rt - (25.0 + rt_shift)) / 2.5) ** 2)
            c = pyopenms.MSChromatogram()
            c.set_peaks((rt.tolist(), y.tolist()))
            c.setNativeID(f"qc_rep_{i}")
            reps.append(c)

        qc_res = evaluate_qc_replicates(
            replicate_chroms=reps,
            feature_id="test_qc_metabolite",
        )
        self.assertEqual(qc_res.num_replicates, 4)
        self.assertGreater(qc_res.raw_apex_rt_cv_pct, 0.0)
        self.assertLess(qc_res.raw_apex_rt_cv_pct, 2.0)
        self.assertGreater(qc_res.raw_apex_int_cv_pct, 0.0)
        self.assertLess(qc_res.raw_apex_int_cv_pct, 5.0)

    def test_graceful_missing_data_handling(self):
        """Verify load_chromatograms_from_file and benchmark handle missing paths gracefully."""
        res_empty = load_chromatograms_from_file("path/to/nonexistent_file.mzML")
        self.assertEqual(res_empty, [])

        res, manifest, overlays = self.engine.run_real_data_benchmark(
            dataset="pass00779",
            data_dir="path/to/nonexistent_pass_dir",
        )
        self.assertIsNone(res)
        self.assertIsNone(manifest)
        self.assertEqual(overlays, [])

    def test_real_data_engine_local(self):
        """Verify engine executes real-data benchmark offline using in-tree fixtures."""
        res, manifest, overlays = self.engine.run_real_data_benchmark(
            dataset="local",
            target_sample_size=10,
        )
        self.assertIsNotNone(res)
        self.assertIsNotNone(manifest)
        self.assertGreater(res.total_traces_evaluated, 0)
        self.assertEqual(manifest.dataset_accession, "LOCAL_FIXTURES")
        self.assertGreater(len(overlays), 0)


if __name__ == "__main__":
    unittest.main()

