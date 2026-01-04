"""
Test suite and examples for Targeted Quantitation algorithms.

This module demonstrates the usage of the SmartPeaks-derived algorithms
for MRM/SRM and SWATH/DIA quantitation workflows.

Classes covered:
- MRMFeatureSelector, MRMFeatureSelectorScore, MRMFeatureSelectorQMIP
- MRMBatchFeatureSelector
- MRMFeatureFilter, MRMFeatureQC
- AbsoluteQuantitation, AbsoluteQuantitationMethod
- PeakIntegrator
- EmgGradientDescent
- TargetedSpectraExtractor
"""

import unittest
import os

import pyopenms as oms


class TestMRMFeatureQC(unittest.TestCase):
    """Test MRMFeatureQC data structure."""

    def test_component_qcs(self):
        """Test ComponentQCs creation and attributes."""
        qc = oms.MRMFeatureQC()

        # Create component QC
        comp_qc = oms.MRMFQC_ComponentQCs()
        comp_qc.component_name = "glucose_quantifier"
        comp_qc.retention_time_l = 5.0
        comp_qc.retention_time_u = 7.0
        comp_qc.intensity_l = 1000.0
        comp_qc.intensity_u = 1e9
        comp_qc.overall_quality_l = 0.5
        comp_qc.overall_quality_u = 1.0

        qc.component_qcs.append(comp_qc)

        self.assertEqual(len(qc.component_qcs), 1)
        self.assertEqual(qc.component_qcs[0].component_name, "glucose_quantifier")
        self.assertAlmostEqual(qc.component_qcs[0].retention_time_l, 5.0)

    def test_component_group_qcs(self):
        """Test ComponentGroupQCs creation with ion ratio."""
        qc = oms.MRMFeatureQC()

        group_qc = oms.MRMFQC_ComponentGroupQCs()
        group_qc.component_group_name = "glucose"
        group_qc.n_heavy_l = 1
        group_qc.n_light_l = 1
        group_qc.ion_ratio_pair_name_1 = "glucose_qualifier"
        group_qc.ion_ratio_pair_name_2 = "glucose_quantifier"
        group_qc.ion_ratio_l = 0.3
        group_qc.ion_ratio_u = 0.5

        qc.component_group_qcs.append(group_qc)

        self.assertEqual(len(qc.component_group_qcs), 1)
        self.assertEqual(qc.component_group_qcs[0].component_group_name, "glucose")

    def test_component_group_pair_qcs(self):
        """Test ComponentGroupPairQCs for isotopic pairs."""
        qc = oms.MRMFeatureQC()

        pair_qc = oms.MRMFQC_ComponentGroupPairQCs()
        pair_qc.component_group_name = "glucose"
        pair_qc.resolution_pair_name = "glucose_13C6"
        pair_qc.resolution_l = 1.5
        pair_qc.resolution_u = 10.0
        pair_qc.rt_diff_l = 0.0
        pair_qc.rt_diff_u = 0.5

        qc.component_group_pair_qcs.append(pair_qc)

        self.assertEqual(len(qc.component_group_pair_qcs), 1)


class TestMRMFeatureFilter(unittest.TestCase):
    """Test MRMFeatureFilter QC filtering."""

    def test_filter_creation(self):
        """Test MRMFeatureFilter instantiation."""
        filter_obj = oms.MRMFeatureFilter()
        params = filter_obj.getParameters()

        # Check that flag_or_filter parameter exists
        self.assertTrue(params.exists("flag_or_filter"))

    def test_calculate_ion_ratio(self):
        """Test ion ratio calculation between features."""
        filter_obj = oms.MRMFeatureFilter()

        # Create two features with peak intensities
        f1 = oms.Feature()
        f1.setMetaValue("peak_apex_int", 1000.0)

        f2 = oms.Feature()
        f2.setMetaValue("peak_apex_int", 500.0)

        ratio = filter_obj.calculateIonRatio(f1, f2, "peak_apex_int")
        self.assertAlmostEqual(ratio, 2.0, places=5)


class TestAbsoluteQuantitationMethod(unittest.TestCase):
    """Test AbsoluteQuantitationMethod configuration."""

    def test_method_configuration(self):
        """Test setting up a quantitation method."""
        method = oms.AbsoluteQuantitationMethod()

        method.setComponentName("glucose")
        method.setFeatureName("peak_apex_int")
        method.setISName("glucose_13C6")
        method.setConcentrationUnits("uM")
        method.setTransformationModel("linear")

        method.setLLOD(0.1)
        method.setULOD(1000.0)
        method.setLLOQ(1.0)
        method.setULOQ(500.0)

        self.assertEqual(method.getComponentName(), "glucose")
        self.assertEqual(method.getISName(), "glucose_13C6")
        self.assertAlmostEqual(method.getLLOQ(), 1.0)
        self.assertAlmostEqual(method.getULOQ(), 500.0)


class TestAbsoluteQuantitation(unittest.TestCase):
    """Test AbsoluteQuantitation workflow."""

    def test_quantitation_setup(self):
        """Test setting up quantitation methods."""
        aq = oms.AbsoluteQuantitation()

        # Create a method
        method = oms.AbsoluteQuantitationMethod()
        method.setComponentName("test_compound")
        method.setFeatureName("peak_apex_int")
        method.setTransformationModel("linear")

        # Set methods
        aq.setQuantMethods([method])

        # Retrieve methods
        methods = aq.getQuantMethods()
        self.assertEqual(len(methods), 1)
        self.assertEqual(methods[0].getComponentName(), "test_compound")

    def test_calculate_ratio(self):
        """Test ratio calculation between analyte and IS."""
        aq = oms.AbsoluteQuantitation()

        # Create analyte feature
        analyte = oms.Feature()
        analyte.setMetaValue("peak_apex_int", 2000.0)

        # Create IS feature
        internal_std = oms.Feature()
        internal_std.setMetaValue("peak_apex_int", 1000.0)

        ratio = aq.calculateRatio(analyte, internal_std, "peak_apex_int")
        self.assertAlmostEqual(ratio, 2.0, places=5)

    def test_calculate_bias(self):
        """Test bias calculation."""
        aq = oms.AbsoluteQuantitation()

        # Bias = (calculated - actual) / actual * 100
        bias = aq.calculateBias(100.0, 110.0)
        self.assertAlmostEqual(bias, 10.0, places=5)


class TestPeakIntegrator(unittest.TestCase):
    """Test PeakIntegrator for peak area calculation."""

    def test_integrator_creation(self):
        """Test PeakIntegrator instantiation and parameters."""
        pi = oms.PeakIntegrator()
        params = pi.getParameters()

        # Check default parameters
        self.assertTrue(params.exists("integration_type"))
        self.assertTrue(params.exists("baseline_type"))

    def test_integrate_chromatogram(self):
        """Test peak integration on a chromatogram."""
        pi = oms.PeakIntegrator()

        # Create a simple Gaussian-like chromatogram
        chrom = oms.MSChromatogram()
        import math

        # Generate Gaussian peak centered at RT=10
        for i in range(20):
            rt = 5.0 + i * 0.5
            intensity = 1000.0 * math.exp(-0.5 * ((rt - 10.0) / 1.0) ** 2)
            peak = oms.ChromatogramPeak()
            peak.setRT(rt)
            peak.setIntensity(intensity)
            chrom.push_back(peak)

        # Integrate peak
        result = pi.integratePeak(chrom, 7.0, 13.0)

        # Check results
        self.assertGreater(result.area, 0)
        self.assertGreater(result.height, 0)
        self.assertGreater(result.apex_pos, 0)

    def test_peak_shape_metrics(self):
        """Test peak shape metrics calculation."""
        pi = oms.PeakIntegrator()

        # Create a chromatogram
        chrom = oms.MSChromatogram()
        import math

        for i in range(40):
            rt = 5.0 + i * 0.25
            intensity = 1000.0 * math.exp(-0.5 * ((rt - 10.0) / 1.0) ** 2)
            peak = oms.ChromatogramPeak()
            peak.setRT(rt)
            peak.setIntensity(intensity)
            chrom.push_back(peak)

        # Get peak area first
        area_result = pi.integratePeak(chrom, 7.0, 13.0)

        # Calculate shape metrics
        metrics = pi.calculatePeakShapeMetrics(
            chrom, 7.0, 13.0, area_result.height, area_result.apex_pos
        )

        # Check metrics
        self.assertGreater(metrics.width_at_50, 0)
        self.assertGreater(metrics.total_width, 0)


class TestEmgGradientDescent(unittest.TestCase):
    """Test EmgGradientDescent for EMG peak fitting."""

    def test_emg_creation(self):
        """Test EmgGradientDescent instantiation."""
        emg = oms.EmgGradientDescent()
        params = oms.Param()
        emg.getDefaultParameters(params)

        # Check that parameters exist
        self.assertTrue(params.exists("max_gd_iter"))


class TestMRMFeatureSelector(unittest.TestCase):
    """Test MRMFeatureSelector LP-based selection."""

    def test_selector_parameters(self):
        """Test SelectorParameters configuration."""
        params = oms.SelectorParameters()

        # Check default values
        self.assertEqual(params.nn_threshold, 4)
        self.assertEqual(params.segment_window_length, 8)
        self.assertEqual(params.segment_step_length, 4)
        self.assertAlmostEqual(params.optimal_threshold, 0.5)

        # Modify parameters
        params.nn_threshold = 6
        params.optimal_threshold = 0.7
        params.select_transition_group = False

        self.assertEqual(params.nn_threshold, 6)
        self.assertAlmostEqual(params.optimal_threshold, 0.7)

    def test_score_selector(self):
        """Test MRMFeatureSelectorScore instantiation."""
        selector = oms.MRMFeatureSelectorScore()
        self.assertIsNotNone(selector)

    def test_qmip_selector(self):
        """Test MRMFeatureSelectorQMIP instantiation."""
        selector = oms.MRMFeatureSelectorQMIP()
        self.assertIsNotNone(selector)


class TestTargetedSpectraExtractor(unittest.TestCase):
    """Test TargetedSpectraExtractor for DDA processing."""

    def test_extractor_creation(self):
        """Test TargetedSpectraExtractor instantiation."""
        extractor = oms.TargetedSpectraExtractor()
        params = oms.Param()
        extractor.getDefaultParameters(params)

        # Check key parameters
        self.assertTrue(params.exists("rt_window"))
        self.assertTrue(params.exists("mz_tolerance"))

    def test_pick_spectrum(self):
        """Test spectrum peak picking."""
        extractor = oms.TargetedSpectraExtractor()

        # Create a simple spectrum
        spectrum = oms.MSSpectrum()
        import math

        for i in range(100):
            mz = 100.0 + i * 1.0
            # Add some peaks
            if i in [20, 50, 80]:
                intensity = 1000.0
            else:
                intensity = 10.0 + 5.0 * math.sin(i * 0.1)

            peak = oms.Peak1D()
            peak.setMZ(mz)
            peak.setIntensity(intensity)
            spectrum.push_back(peak)

        picked = oms.MSSpectrum()
        extractor.pickSpectrum(spectrum, picked)

        # Should have picked some peaks
        self.assertGreater(len(picked), 0)


# Example usage demonstration
def example_targeted_quantitation_workflow():
    """
    Complete example of targeted quantitation workflow.

    This demonstrates the typical workflow for MRM/SRM quantitation:
    1. Load chromatograms and transitions
    2. Pick peaks with MRMTransitionGroupPicker
    3. Integrate peaks with PeakIntegrator
    4. Apply QC filtering with MRMFeatureFilter
    5. Select optimal features with MRMFeatureSelector
    6. Quantify with AbsoluteQuantitation
    """
    print("=== Targeted Quantitation Workflow Example ===\n")

    # Step 1: Set up QC criteria
    print("Step 1: Setting up QC criteria...")
    qc = oms.MRMFeatureQC()

    # Component-level QC
    comp_qc = oms.MRMFQC_ComponentQCs()
    comp_qc.component_name = "glucose_quantifier"
    comp_qc.retention_time_l = 5.0
    comp_qc.retention_time_u = 7.0
    comp_qc.intensity_l = 1000.0
    comp_qc.intensity_u = 1e9
    qc.component_qcs.append(comp_qc)
    print(f"  Added QC for: {comp_qc.component_name}")

    # Step 2: Configure quantitation method
    print("\nStep 2: Configuring quantitation method...")
    method = oms.AbsoluteQuantitationMethod()
    method.setComponentName("glucose")
    method.setFeatureName("peak_apex_int")
    method.setISName("glucose_13C6")
    method.setConcentrationUnits("uM")
    method.setTransformationModel("linear")
    method.setLLOQ(1.0)
    method.setULOQ(500.0)
    print(f"  Component: {method.getComponentName()}")
    print(f"  Internal Standard: {method.getISName()}")
    print(f"  LLOQ-ULOQ: {method.getLLOQ()}-{method.getULOQ()} {method.getConcentrationUnits()}")

    # Step 3: Set up AbsoluteQuantitation
    print("\nStep 3: Setting up AbsoluteQuantitation...")
    aq = oms.AbsoluteQuantitation()
    aq.setQuantMethods([method])
    print(f"  Configured {len(aq.getQuantMethods())} quantitation method(s)")

    # Step 4: Configure feature selector
    print("\nStep 4: Configuring feature selector...")
    params = oms.SelectorParameters()
    params.nn_threshold = 4
    params.segment_window_length = 8
    params.optimal_threshold = 0.5
    print(f"  Nearest neighbors: {params.nn_threshold}")
    print(f"  Window length: {params.segment_window_length}")
    print(f"  Optimal threshold: {params.optimal_threshold}")

    # Step 5: Demonstrate PeakIntegrator
    print("\nStep 5: Demonstrating PeakIntegrator...")
    pi = oms.PeakIntegrator()
    pi_params = pi.getParameters()
    pi_params.setValue("integration_type", "trapezoid")
    pi_params.setValue("baseline_type", "base_to_base")
    pi.setParameters(pi_params)
    print(f"  Integration type: {pi_params.getValue('integration_type')}")
    print(f"  Baseline type: {pi_params.getValue('baseline_type')}")

    print("\n=== Workflow setup complete! ===")
    print("\nIn a real workflow, you would now:")
    print("  1. Load your mzML/featureXML data")
    print("  2. Run MRMTransitionGroupPicker for peak picking")
    print("  3. Apply MRMFeatureFilter with your QC criteria")
    print("  4. Use MRMFeatureSelector for optimal feature selection")
    print("  5. Call aq.quantifyComponents() on your unknowns")


if __name__ == "__main__":
    # Run the example
    example_targeted_quantitation_workflow()

    # Run unit tests
    print("\n\n=== Running Unit Tests ===\n")
    unittest.main(verbosity=2)
