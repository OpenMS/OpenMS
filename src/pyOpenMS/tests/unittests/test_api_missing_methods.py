"""
Tests for missing methods on existing classes (Task 3).

Verifies that methods present in the Cython 3.5 API
are available in the nanobind build.
"""
import pytest
import pyopenms


# === Tier 1: Large gaps (9+ missing methods) ===

class TestInstrumentMethods:
    """Instrument is missing 19 of its own methods (only has MetaInfoInterface)."""

    def test_name(self):
        inst = pyopenms.Instrument()
        inst.setName("test_instrument")
        assert inst.getName() == "test_instrument"

    def test_vendor(self):
        inst = pyopenms.Instrument()
        inst.setVendor("Thermo")
        assert inst.getVendor() == "Thermo"

    def test_model(self):
        inst = pyopenms.Instrument()
        inst.setModel("Orbitrap")
        assert inst.getModel() == "Orbitrap"

    def test_customizations(self):
        inst = pyopenms.Instrument()
        inst.setCustomizations("custom")
        assert inst.getCustomizations() == "custom"

    def test_ion_sources(self):
        inst = pyopenms.Instrument()
        sources = inst.getIonSources()
        assert isinstance(sources, list)
        inst.setIonSources(sources)

    def test_mass_analyzers(self):
        inst = pyopenms.Instrument()
        analyzers = inst.getMassAnalyzers()
        assert isinstance(analyzers, list)
        inst.setMassAnalyzers(analyzers)

    def test_ion_detectors(self):
        inst = pyopenms.Instrument()
        detectors = inst.getIonDetectors()
        assert isinstance(detectors, list)
        inst.setIonDetectors(detectors)

    def test_software(self):
        inst = pyopenms.Instrument()
        sw = inst.getSoftware()
        assert sw is not None
        inst.setSoftware(sw)

    def test_ion_optics(self):
        inst = pyopenms.Instrument()
        optics = inst.getIonOptics()
        assert optics is not None
        inst.setIonOptics(optics)


class TestNASequenceMethods:
    """NASequence missing 11 core methods."""

    def test_get(self):
        seq = pyopenms.NASequence.fromString("ACGU")
        r = seq.get(0)
        assert r is not None

    def test_get_prefix(self):
        seq = pyopenms.NASequence.fromString("ACGU")
        prefix = seq.getPrefix(2)
        assert prefix is not None

    def test_get_suffix(self):
        seq = pyopenms.NASequence.fromString("ACGU")
        suffix = seq.getSuffix(2)
        assert suffix is not None

    def test_get_subsequence(self):
        seq = pyopenms.NASequence.fromString("ACGU")
        sub = seq.getSubsequence(1, 2)
        assert sub is not None

    def test_set(self):
        seq = pyopenms.NASequence.fromString("ACGU")
        assert hasattr(seq, 'set')

    def test_five_prime_mod(self):
        seq = pyopenms.NASequence.fromString("ACGU")
        mod = seq.getFivePrimeMod()
        # May be None if no mod
        assert hasattr(seq, 'setFivePrimeMod')

    def test_three_prime_mod(self):
        seq = pyopenms.NASequence.fromString("ACGU")
        mod = seq.getThreePrimeMod()
        assert hasattr(seq, 'setThreePrimeMod')


class TestSpectrumMetaDataLookupMethods:
    """SpectrumMetaDataLookup missing 12 methods."""

    def test_exists(self):
        assert hasattr(pyopenms, 'SpectrumMetaDataLookup')

    def test_read_spectra(self):
        lookup = pyopenms.SpectrumMetaDataLookup()
        assert hasattr(lookup, 'readSpectra')

    def test_find_by_native_id(self):
        lookup = pyopenms.SpectrumMetaDataLookup()
        assert hasattr(lookup, 'findByNativeID')

    def test_find_by_rt(self):
        lookup = pyopenms.SpectrumMetaDataLookup()
        assert hasattr(lookup, 'findByRT')

    def test_find_by_index(self):
        lookup = pyopenms.SpectrumMetaDataLookup()
        assert hasattr(lookup, 'findByIndex')


# === Tier 2: Medium gaps ===

class TestMSExperimentMethods:
    """MSExperiment/PeakMap missing get2DPeakData*, getSpectrum, get_df_columns."""

    def test_get_spectrum_by_index(self):
        """getSpectrum(index) should return spectrum by index."""
        exp = pyopenms.MSExperiment()
        spec = pyopenms.MSSpectrum()
        exp.addSpectrum(spec)
        s = exp.getSpectrum(0)
        assert s is not None

    @pytest.mark.xfail(reason="get2DPeakData needs Python addon implementation")
    def test_get_2d_peak_data(self):
        exp = pyopenms.MSExperiment()
        assert hasattr(exp, 'get2DPeakData')

    @pytest.mark.xfail(reason="get2DPeakDataIM needs Python addon implementation")
    def test_get_2d_peak_data_im(self):
        exp = pyopenms.MSExperiment()
        assert hasattr(exp, 'get2DPeakDataIM')

    @pytest.mark.xfail(reason="get_df_columns deprecated, use df_columns instead")
    def test_get_df_columns(self):
        exp = pyopenms.MSExperiment()
        assert hasattr(exp, 'get_df_columns')


class TestMSSpectrumMethods:
    """MSSpectrum missing findHighestInWindow, select, intensityInRange, get_df_columns."""

    def test_find_highest_in_window(self):
        spec = pyopenms.MSSpectrum()
        assert hasattr(spec, 'findHighestInWindow')

    def test_select(self):
        spec = pyopenms.MSSpectrum()
        assert hasattr(spec, 'select')

    @pytest.mark.xfail(reason="intensityInRange does not exist in C++ MSSpectrum")
    def test_intensity_in_range(self):
        spec = pyopenms.MSSpectrum()
        assert hasattr(spec, 'intensityInRange')

    @pytest.mark.xfail(reason="get_df_columns deprecated, use df_columns instead")
    def test_get_df_columns(self):
        spec = pyopenms.MSSpectrum()
        assert hasattr(spec, 'get_df_columns')


class TestMSChromatogramMethods:
    """MSChromatogram missing clearRanges, get_df_columns, reserve, resize."""

    def test_reserve(self):
        chrom = pyopenms.MSChromatogram()
        assert hasattr(chrom, 'reserve')

    def test_resize(self):
        chrom = pyopenms.MSChromatogram()
        assert hasattr(chrom, 'resize')

    def test_clear_ranges(self):
        chrom = pyopenms.MSChromatogram()
        assert hasattr(chrom, 'clearRanges')

    @pytest.mark.xfail(reason="get_df_columns deprecated, use df_columns instead")
    def test_get_df_columns(self):
        chrom = pyopenms.MSChromatogram()
        assert hasattr(chrom, 'get_df_columns')


class TestSetCVTerms:
    """setCVTerms is missing on 14 CVTermList-derived classes."""

    @pytest.mark.parametrize("cls_name", [
        "CVTermList", "CVTermListInterface", "Compound", "Configuration",
        "Contact", "IncludeExcludeTarget", "Peptide", "Precursor",
        "Prediction", "Protein", "Publication",
        "ReactionMonitoringTransition", "RetentionTime", "TraMLProduct",
    ])
    def test_set_cv_terms(self, cls_name):
        cls = getattr(pyopenms, cls_name)
        obj = cls()
        assert hasattr(obj, 'setCVTerms'), \
            f"{cls_name} missing setCVTerms"


class TestTargetedExperimentMethods:
    """TargetedExperiment missing contact/instrument/publication management."""

    def test_add_contact(self):
        te = pyopenms.TargetedExperiment()
        assert hasattr(te, 'addContact')

    def test_get_contacts(self):
        te = pyopenms.TargetedExperiment()
        assert hasattr(te, 'getContacts')

    def test_set_contacts(self):
        te = pyopenms.TargetedExperiment()
        assert hasattr(te, 'setContacts')

    def test_add_instrument(self):
        te = pyopenms.TargetedExperiment()
        assert hasattr(te, 'addInstrument')

    def test_get_instruments(self):
        te = pyopenms.TargetedExperiment()
        assert hasattr(te, 'getInstruments')

    def test_add_publication(self):
        te = pyopenms.TargetedExperiment()
        assert hasattr(te, 'addPublication')


class TestOnDiscMSExperimentMethods:
    """OnDiscMSExperiment missing getSpectrumById, getChromatogramById, etc."""

    def test_get_spectrum_by_id(self):
        assert hasattr(pyopenms.OnDiscMSExperiment, 'getSpectrumById')

    def test_get_chromatogram_by_id(self):
        assert hasattr(pyopenms.OnDiscMSExperiment, 'getChromatogramById')

    def test_get_experimental_settings(self):
        assert hasattr(pyopenms.OnDiscMSExperiment, 'getExperimentalSettings')

    def test_get_meta_data(self):
        assert hasattr(pyopenms.OnDiscMSExperiment, 'getMetaData')


class TestConvexHull2DMethods:
    """ConvexHull2D missing numpy interop methods."""

    def test_add_points_npy(self):
        hull = pyopenms.ConvexHull2D()
        assert hasattr(hull, 'addPointsNPY')

    def test_get_hull_points_npy(self):
        hull = pyopenms.ConvexHull2D()
        assert hasattr(hull, 'getHullPointsNPY')

    def test_set_hull_points_npy(self):
        hull = pyopenms.ConvexHull2D()
        assert hasattr(hull, 'setHullPointsNPY')

    def test_encloses_xy(self):
        hull = pyopenms.ConvexHull2D()
        assert hasattr(hull, 'enclosesXY')

    def test_get_bounding_box_2d(self):
        hull = pyopenms.ConvexHull2D()
        assert hasattr(hull, 'getBoundingBox2D')


class TestDataValueMethods:
    """DataValue missing unit getter/setter methods.
    DataValue uses a nanobind type caster (transparent C++ <-> Python conversion),
    so there is no nb::class_<DataValue> to add methods to. Unit info is lost in conversion."""

    @pytest.mark.xfail(reason="DataValue uses type caster, unit info lost in conversion")
    def test_get_unit(self):
        dv = pyopenms.DataValue(42)
        assert hasattr(dv, 'getUnit')

    @pytest.mark.xfail(reason="DataValue uses type caster, unit info lost in conversion")
    def test_set_unit(self):
        dv = pyopenms.DataValue(42)
        assert hasattr(dv, 'setUnit')

    @pytest.mark.xfail(reason="DataValue uses type caster, unit info lost in conversion")
    def test_get_unit_type(self):
        dv = pyopenms.DataValue(42)
        assert hasattr(dv, 'getUnitType')

    @pytest.mark.xfail(reason="DataValue uses type caster, unit info lost in conversion")
    def test_set_unit_type(self):
        dv = pyopenms.DataValue(42)
        assert hasattr(dv, 'setUnitType')

    @pytest.mark.xfail(reason="DataValue uses type caster, unit info lost in conversion")
    def test_has_unit(self):
        dv = pyopenms.DataValue(42)
        assert hasattr(dv, 'hasUnit')


class TestConsensusMapMethods:
    """ConsensusMap missing clearRanges, empty, reserve."""

    def test_empty(self):
        cm = pyopenms.ConsensusMap()
        assert hasattr(cm, 'empty')
        assert cm.empty()

    def test_reserve(self):
        cm = pyopenms.ConsensusMap()
        assert hasattr(cm, 'reserve')

    def test_clear_ranges(self):
        cm = pyopenms.ConsensusMap()
        assert hasattr(cm, 'clearRanges')


class TestFeatureMapMethods:
    """FeatureMap missing clearRanges, get_df_columns."""

    def test_clear_ranges(self):
        fm = pyopenms.FeatureMap()
        assert hasattr(fm, 'clearRanges')

    @pytest.mark.xfail(reason="get_df_columns deprecated, use df_columns instead")
    def test_get_df_columns(self):
        fm = pyopenms.FeatureMap()
        assert hasattr(fm, 'get_df_columns')


class TestMetaInfoInterfaceMethods:
    """MetaInfoInterface missing getMetaValues, setMetaValues."""

    def test_get_meta_values(self):
        mii = pyopenms.MetaInfoInterface()
        assert hasattr(mii, 'getMetaValues')

    def test_set_meta_values(self):
        mii = pyopenms.MetaInfoInterface()
        assert hasattr(mii, 'setMetaValues')


class TestFASTAFileMethods:
    """FASTAFile missing readNext, writeNext."""

    def test_read_next(self):
        ff = pyopenms.FASTAFile()
        assert hasattr(ff, 'readNext')

    def test_write_next(self):
        ff = pyopenms.FASTAFile()
        assert hasattr(ff, 'writeNext')


class TestBase64Methods:
    """Base64 missing encode/decode methods."""

    def test_encode64(self):
        assert hasattr(pyopenms.Base64, 'encode64')

    def test_decode64(self):
        assert hasattr(pyopenms.Base64, 'decode64')

    def test_encode32(self):
        assert hasattr(pyopenms.Base64, 'encode32')

    def test_decode32(self):
        assert hasattr(pyopenms.Base64, 'decode32')


class TestDeisotoperMethods:
    """Deisotoper missing convenience methods."""

    def test_deisotope_and_single_charge_default(self):
        assert hasattr(pyopenms.Deisotoper, 'deisotopeAndSingleChargeDefault')

    def test_deisotope_with_averagine_model(self):
        assert hasattr(pyopenms.Deisotoper, 'deisotopeWithAveragineModel')


class TestResidueModificationMethods:
    """ResidueModification missing getDiffFormula."""

    def test_get_diff_formula(self):
        mod = pyopenms.ResidueModification()
        assert hasattr(mod, 'getDiffFormula')


class TestProteinIdentificationMethods:
    """ProteinIdentification missing computeCoverage, setPrimaryMSRunPath."""

    def test_compute_coverage(self):
        pi = pyopenms.ProteinIdentification()
        assert hasattr(pi, 'computeCoverage')

    def test_set_primary_ms_run_path(self):
        pi = pyopenms.ProteinIdentification()
        assert hasattr(pi, 'setPrimaryMSRunPath')


class TestPeptideHitMethods:
    """PeptideHit missing setAnalysisResults."""

    def test_set_analysis_results(self):
        ph = pyopenms.PeptideHit()
        assert hasattr(ph, 'setAnalysisResults')


class TestMobilogramMethods:
    """Mobilogram missing data array and range methods."""

    def test_find_nearest(self):
        mob = pyopenms.Mobilogram()
        assert hasattr(mob, 'findNearest')

    def test_get_integer_data_arrays(self):
        mob = pyopenms.Mobilogram()
        assert hasattr(mob, 'getIntegerDataArrays')

    def test_set_integer_data_arrays(self):
        mob = pyopenms.Mobilogram()
        assert hasattr(mob, 'setIntegerDataArrays')

    def test_get_string_data_arrays(self):
        mob = pyopenms.Mobilogram()
        assert hasattr(mob, 'getStringDataArrays')

    def test_set_string_data_arrays(self):
        mob = pyopenms.Mobilogram()
        assert hasattr(mob, 'setStringDataArrays')


class TestTextFileMethods:
    """TextFile missing addLine."""

    def test_add_line(self):
        tf = pyopenms.TextFile()
        assert hasattr(tf, 'addLine')


class TestMassTraceDetectionMethods:
    """MassTraceDetection missing run."""

    def test_run(self):
        mtd = pyopenms.MassTraceDetection()
        assert hasattr(mtd, 'run')


class TestFeatureFinderAlgorithmMetaboIdentMethods:
    """FeatureFinderAlgorithmMetaboIdent missing run."""

    def test_run(self):
        ff = pyopenms.FeatureFinderAlgorithmMetaboIdent()
        assert hasattr(ff, 'run')


# === Tier 3: Additional gaps filled ===

class TestLinearInterpolationMethods:
    """LinearInterpolation missing 17 methods."""

    def test_value(self):
        li = pyopenms.LinearInterpolation()
        assert hasattr(li, 'value')

    def test_add_value(self):
        li = pyopenms.LinearInterpolation()
        assert hasattr(li, 'addValue')

    def test_get_set_data(self):
        li = pyopenms.LinearInterpolation()
        assert hasattr(li, 'getData')
        assert hasattr(li, 'setData')

    def test_get_set_scale(self):
        li = pyopenms.LinearInterpolation()
        assert hasattr(li, 'getScale')
        assert hasattr(li, 'setScale')

    def test_get_set_offset(self):
        li = pyopenms.LinearInterpolation()
        assert hasattr(li, 'getOffset')
        assert hasattr(li, 'setOffset')

    def test_set_mapping(self):
        li = pyopenms.LinearInterpolation()
        assert hasattr(li, 'setMapping')

    def test_support_range(self):
        li = pyopenms.LinearInterpolation()
        assert hasattr(li, 'supportMin')
        assert hasattr(li, 'supportMax')

    def test_key_index_transform(self):
        li = pyopenms.LinearInterpolation()
        assert hasattr(li, 'key2index')
        assert hasattr(li, 'index2key')


class TestDistanceMatrixMethods:
    """DistanceMatrix missing 9 methods."""

    def test_core_operations(self):
        dm = pyopenms.DistanceMatrix(3, 0.0)
        dm.setValue(0, 1, 1.5)
        assert dm.getValue(0, 1) == pytest.approx(1.5)
        assert dm.dimensionsize() == 3

    def test_min_element(self):
        dm = pyopenms.DistanceMatrix(3, 0.0)
        dm.setValue(0, 1, 1.5)
        dm.setValue(0, 2, 0.5)
        dm.updateMinElement()
        coords = dm.getMinElementCoordinates()
        assert coords is not None

    def test_clear_resize(self):
        dm = pyopenms.DistanceMatrix(3, 0.0)
        assert hasattr(dm, 'clear')
        assert hasattr(dm, 'resize')
        assert hasattr(dm, 'reduce')


class TestCVMappingsMethods:
    """CVMappings missing 7 methods."""

    def test_cv_references(self):
        cvm = pyopenms.CVMappings()
        refs = cvm.getCVReferences()
        assert isinstance(refs, list)
        assert hasattr(cvm, 'setCVReferences')
        assert hasattr(cvm, 'addCVReference')
        assert hasattr(cvm, 'hasCVReference')

    def test_mapping_rules(self):
        cvm = pyopenms.CVMappings()
        rules = cvm.getMappingRules()
        assert isinstance(rules, list)
        assert hasattr(cvm, 'setMappingRules')
        assert hasattr(cvm, 'addMappingRule')


class TestCompomerMethods:
    """Compomer missing 4 methods."""

    def test_get_labels(self):
        c = pyopenms.Compomer()
        assert hasattr(c, 'getLabels')

    def test_is_conflicting(self):
        c = pyopenms.Compomer()
        assert hasattr(c, 'isConflicting')

    def test_is_single_adduct(self):
        c = pyopenms.Compomer()
        assert hasattr(c, 'isSingleAdduct')

    def test_remove_adduct(self):
        c = pyopenms.Compomer()
        assert hasattr(c, 'removeAdduct')


class TestCoarseIsotopePatternGeneratorFragmentMethods:
    """CoarseIsotopePatternGenerator missing 6 fragment methods."""

    def test_fragment_methods(self):
        gen = pyopenms.CoarseIsotopePatternGenerator()
        assert hasattr(gen, 'calcFragmentIsotopeDist')
        assert hasattr(gen, 'estimateForFragmentFromPeptideWeight')
        assert hasattr(gen, 'estimateForFragmentFromPeptideWeightAndS')
        assert hasattr(gen, 'estimateForFragmentFromDNAWeight')
        assert hasattr(gen, 'estimateForFragmentFromRNAWeight')
        assert hasattr(gen, 'estimateForFragmentFromWeightAndComp')


class TestExperimentalDesignMethods:
    """ExperimentalDesign missing 4 section methods."""

    def test_ms_file_section(self):
        ed = pyopenms.ExperimentalDesign()
        assert hasattr(ed, 'getMSFileSection')
        assert hasattr(ed, 'setMSFileSection')

    def test_sample_section(self):
        ed = pyopenms.ExperimentalDesign()
        assert hasattr(ed, 'getSampleSection')
        assert hasattr(ed, 'setSampleSection')


class TestControlledVocabularyMethods:
    """ControlledVocabulary missing 3 methods."""

    def test_get_term(self):
        cv = pyopenms.ControlledVocabulary()
        assert hasattr(cv, 'getTerm')

    def test_get_term_by_name(self):
        cv = pyopenms.ControlledVocabulary()
        assert hasattr(cv, 'getTermByName')

    def test_get_all_child_terms(self):
        cv = pyopenms.ControlledVocabulary()
        assert hasattr(cv, 'getAllChildTerms')


class TestFalseDiscoveryRateMethods:
    """FalseDiscoveryRate missing applyBasic, applyEstimated, rocN."""

    def test_apply_basic(self):
        fdr = pyopenms.FalseDiscoveryRate()
        assert hasattr(fdr, 'applyBasic')

    def test_apply_estimated(self):
        fdr = pyopenms.FalseDiscoveryRate()
        assert hasattr(fdr, 'applyEstimated')

    def test_roc_n(self):
        fdr = pyopenms.FalseDiscoveryRate()
        assert hasattr(fdr, 'rocN')


class TestAbsoluteQuantitationMethods:
    """AbsoluteQuantitation missing calculateBiasAndR, fitCalibration."""

    def test_fit_calibration(self):
        aq = pyopenms.AbsoluteQuantitation()
        assert hasattr(aq, 'fitCalibration')

    def test_calculate_bias_and_r(self):
        aq = pyopenms.AbsoluteQuantitation()
        assert hasattr(aq, 'calculateBiasAndR')


class TestEmpiricalFormulaMethods:
    """EmpiricalFormula missing __iadd__, getConditionalFragmentIsotopeDist."""

    def test_iadd(self):
        f1 = pyopenms.EmpiricalFormula("C6H12O6")
        f2 = pyopenms.EmpiricalFormula("H2O")
        f1 += f2
        assert f1.toString() != ""

    def test_get_conditional_fragment_isotope_dist(self):
        ef = pyopenms.EmpiricalFormula("C6H12O6")
        assert hasattr(ef, 'getConditionalFragmentIsotopeDist')


class TestAASequenceIadd:
    """AASequence missing __iadd__."""

    def test_iadd(self):
        seq = pyopenms.AASequence.fromString("PEP")
        seq2 = pyopenms.AASequence.fromString("TIDE")
        seq += seq2
        assert seq.toString() == "PEPTIDE"


class TestSpectraMergerMethods:
    """SpectraMerger missing average, mergeSpectraBlockWise, mergeSpectraPrecursors."""

    def test_methods(self):
        sm = pyopenms.SpectraMerger()
        assert hasattr(sm, 'average')
        assert hasattr(sm, 'mergeSpectraBlockWise')
        assert hasattr(sm, 'mergeSpectraPrecursors')


class TestMRMTransitionGroupPickerMethods:
    """MRMTransitionGroupPicker missing findLargestPeak, findWidestPeakIndices."""

    def test_find_methods(self):
        picker = pyopenms.MRMTransitionGroupPicker()
        assert hasattr(picker, 'findLargestPeak')
        assert hasattr(picker, 'findWidestPeakIndices')


class TestSignalToNoiseEstimatorMethods:
    """SignalToNoise estimators missing methods."""

    def test_mean_iterative(self):
        sne = pyopenms.SignalToNoiseEstimatorMeanIterative()
        assert hasattr(sne, 'getSignalToNoise')
        assert hasattr(sne, 'init')

    def test_median(self):
        sne = pyopenms.SignalToNoiseEstimatorMedian()
        assert hasattr(sne, 'getSparseWindowPercent')
        assert hasattr(sne, 'getHistogramRightmostPercent')


class TestPeptideIdentificationListMethods:
    """PeptideIdentificationList missing container methods."""

    def test_iter(self):
        pil = pyopenms.PeptideIdentificationList()
        pi = pyopenms.PeptideIdentification()
        pil.append(pi)
        count = sum(1 for _ in pil)
        assert count == 1

    def test_setitem(self):
        pil = pyopenms.PeptideIdentificationList()
        pi = pyopenms.PeptideIdentification()
        pil.append(pi)
        pil[0] = pi  # __setitem__

    def test_at_front_back(self):
        pil = pyopenms.PeptideIdentificationList()
        pi = pyopenms.PeptideIdentification()
        pil.append(pi)
        assert pil.at(0) is not None
        assert pil.front() is not None
        assert pil.back() is not None


class TestParamNodeMethods:
    """ParamNode missing findEntryRecursive, findParentOf, insert."""

    def test_methods(self):
        pn = pyopenms.ParamNode()
        assert hasattr(pn, 'findEntryRecursive')
        assert hasattr(pn, 'findParentOf')
        assert hasattr(pn, 'insert')


class TestIDFilterMethods:
    """IDFilter missing countHits."""

    def test_count_hits(self):
        assert hasattr(pyopenms.IDFilter, 'countHits')


class TestAccurateMassSearchResultMethods:
    """AccurateMassSearchResult missing setIndividualIntensities."""

    def test_set_individual_intensities(self):
        r = pyopenms.AccurateMassSearchResult()
        r.setIndividualIntensities([1.0, 2.0, 3.0])
        assert r.getIndividualIntensities() == [1.0, 2.0, 3.0]


class TestSpectrumAccessMethods:
    """SpectrumAccess classes missing getChromatogramNativeID, getSpectraByRT."""

    def test_openms(self):
        assert hasattr(pyopenms.SpectrumAccessOpenMS, 'getChromatogramNativeID')
        assert hasattr(pyopenms.SpectrumAccessOpenMS, 'getSpectraByRT')

    def test_in_memory(self):
        assert hasattr(pyopenms.SpectrumAccessOpenMSInMemory, 'getChromatogramNativeID')
        assert hasattr(pyopenms.SpectrumAccessOpenMSInMemory, 'getSpectraByRT')


class TestFloatIntegerDataArrayMethods:
    """Float/IntegerDataArray missing reserve, getDataProcessing."""

    def test_float_reserve(self):
        fda = pyopenms.FloatDataArray()
        fda.reserve(100)
        assert hasattr(fda, 'getDataProcessing')

    def test_integer_reserve(self):
        ida = pyopenms.IntegerDataArray()
        ida.reserve(100)
        assert hasattr(ida, 'getDataProcessing')
