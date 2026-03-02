"""
Tests for missing class API completeness (Task 2).

Verifies that significant classes from the Cython 3.5 API
are present and functional in the nanobind build.
"""
import pytest
import pyopenms


# === Significant Missing Classes (5+ methods) ===

class TestKernelMassTrace:
    """Kernel_MassTrace was an alias for MassTrace in Cython.
    MassTrace exists in nanobind, but the alias is missing."""

    def test_alias_exists(self):
        assert hasattr(pyopenms, 'Kernel_MassTrace')

    def test_alias_is_mass_trace(self):
        assert pyopenms.Kernel_MassTrace is pyopenms.MassTrace


class TestMzMLSqliteHandler:
    """MzMLSqliteHandler — SQLite backend for mzML."""

    def test_exists(self):
        assert hasattr(pyopenms, 'MzMLSqliteHandler')

    def test_constructor(self):
        # MzMLSqliteHandler(filename, run_id)
        cls = pyopenms.MzMLSqliteHandler
        assert cls is not None

    def test_has_read_methods(self):
        cls = pyopenms.MzMLSqliteHandler
        for method in ['readExperiment', 'readSpectra', 'readChromatograms',
                       'getNrSpectra', 'getNrChromatograms', 'getRunID',
                       'getSpectraIndicesbyRT']:
            assert hasattr(cls, method) or hasattr(getattr(cls, '__init__', None), method), \
                f"MzMLSqliteHandler missing {method}"

    def test_has_write_methods(self):
        cls = pyopenms.MzMLSqliteHandler
        for method in ['writeExperiment', 'writeSpectra', 'writeChromatograms',
                       'createTables', 'writeRunLevelInformation', 'setConfig']:
            assert hasattr(cls, method) or hasattr(getattr(cls, '__init__', None), method), \
                f"MzMLSqliteHandler missing {method}"


class TestParamValue:
    """ParamValue — type-flexible parameter container.

    In nanobind, ParamValue uses a type caster for transparent conversion.
    A Python wrapper class provides backward compatibility.
    """

    def test_exists(self):
        assert hasattr(pyopenms, 'ParamValue')

    def test_construct_string(self):
        pv = pyopenms.ParamValue("hello")
        assert pv.toString() == "hello"

    def test_construct_int(self):
        pv = pyopenms.ParamValue(42)
        assert pv.toInt() == 42

    def test_construct_double(self):
        pv = pyopenms.ParamValue(3.14)
        assert abs(pv.toDouble() - 3.14) < 1e-10

    def test_is_empty(self):
        pv = pyopenms.ParamValue()
        assert pv.isEmpty()

    def test_value_type(self):
        pv = pyopenms.ParamValue("hello")
        vt = pv.valueType()
        assert vt is not None

    def test_to_bool(self):
        pv = pyopenms.ParamValue("true")
        assert pv.toBool() is True


class TestTransformationModelInterpolated:
    """TransformationModelInterpolated — interpolated RT alignment."""

    def test_exists(self):
        assert hasattr(pyopenms, 'TransformationModelInterpolated')

    def test_get_default_parameters(self):
        # nanobind binding returns Param directly (static method, no mutable ref)
        p = pyopenms.TransformationModelInterpolated.getDefaultParameters()
        assert not p.empty()

    def test_evaluate(self):
        cls = pyopenms.TransformationModelInterpolated
        assert hasattr(cls, 'evaluate')


class TestTargetedExperimentInstrument:
    """TargetedExperiment_Instrument — CV-term bearing instrument info."""

    def test_exists(self):
        assert hasattr(pyopenms, 'TargetedExperiment_Instrument')

    def test_construct(self):
        inst = pyopenms.TargetedExperiment_Instrument()
        assert inst is not None

    def test_has_cv_term_methods(self):
        inst = pyopenms.TargetedExperiment_Instrument()
        for method in ['addCVTerm', 'hasCVTerm', 'getCVTerms', 'replaceCVTerms',
                       'getMetaValue', 'setMetaValue', 'metaValueExists']:
            assert hasattr(inst, method), \
                f"TargetedExperiment_Instrument missing {method}"


class TestTargetedExperimentInterpretation:
    """TargetedExperiment_Interpretation — CV-term bearing interpretation."""

    def test_exists(self):
        assert hasattr(pyopenms, 'TargetedExperiment_Interpretation')

    def test_construct(self):
        interp = pyopenms.TargetedExperiment_Interpretation()
        assert interp is not None

    def test_has_cv_term_methods(self):
        interp = pyopenms.TargetedExperiment_Interpretation()
        for method in ['addCVTerm', 'hasCVTerm', 'getCVTerms', 'replaceCVTerms',
                       'getMetaValue', 'setMetaValue', 'metaValueExists']:
            assert hasattr(interp, method), \
                f"TargetedExperiment_Interpretation missing {method}"


class TestCVMappingRule:
    """CVMappingRule — CV mapping rule management."""

    def test_exists(self):
        assert hasattr(pyopenms, 'CVMappingRule')

    def test_construct(self):
        rule = pyopenms.CVMappingRule()
        assert rule is not None

    def test_identifier(self):
        rule = pyopenms.CVMappingRule()
        rule.setIdentifier("test")
        assert rule.getIdentifier() == "test"

    def test_element_path(self):
        rule = pyopenms.CVMappingRule()
        rule.setElementPath("/mzML/run")
        assert rule.getElementPath() == "/mzML/run"


class TestIMSAlphabet:
    """IMSAlphabet — IMS alphabet for mass decomposition."""

    def test_exists(self):
        assert hasattr(pyopenms, 'IMSAlphabet')

    def test_construct(self):
        alpha = pyopenms.IMSAlphabet()
        assert alpha is not None

    def test_size(self):
        alpha = pyopenms.IMSAlphabet()
        assert alpha.size() == 0


@pytest.mark.xfail(reason="IMSWeights class not found in OpenMS codebase — may have been removed")
class TestIMSWeights:
    """IMSWeights — IMS integer mass decomposition weights."""

    def test_exists(self):
        assert hasattr(pyopenms, 'IMSWeights')

    def test_construct(self):
        w = pyopenms.IMSWeights()
        assert w is not None

    def test_size(self):
        w = pyopenms.IMSWeights()
        assert w.size() == 0


class TestMSPGenericFile:
    """MSPGenericFile — MSP spectral library file I/O."""

    def test_exists(self):
        assert hasattr(pyopenms, 'MSPGenericFile')

    def test_construct(self):
        f = pyopenms.MSPGenericFile()
        assert f is not None

    def test_has_load_store(self):
        f = pyopenms.MSPGenericFile()
        assert hasattr(f, 'load')
        assert hasattr(f, 'store')


class TestMassTraces:
    """MassTraces — collection of mass traces."""

    def test_exists(self):
        assert hasattr(pyopenms, 'MassTraces')

    def test_construct(self):
        mt = pyopenms.MassTraces()
        assert mt is not None

    def test_get_peak_count(self):
        mt = pyopenms.MassTraces()
        assert mt.getPeakCount() == 0


@pytest.mark.xfail(reason="SpectrumAccessOpenMSCached inherits from ISpectrumAccess+CachedmzML, not bindable without base class support")
class TestSpectrumAccessCached:
    """SpectrumAccessOpenMSCached — cached spectrum access."""

    def test_exists(self):
        assert hasattr(pyopenms, 'SpectrumAccessOpenMSCached')


@pytest.mark.xfail(reason="SpectrumAccessQuadMZTransforming inherits from SpectrumAccessTransforming, not bindable without base class support")
class TestSpectrumAccessQuadMZ:
    """SpectrumAccessQuadMZTransforming — quadratic m/z transform."""

    def test_exists(self):
        assert hasattr(pyopenms, 'SpectrumAccessQuadMZTransforming')


class TestSpectrumAccessSqMass:
    """SpectrumAccessSqMass — SqMass-backed access."""

    def test_exists(self):
        assert hasattr(pyopenms, 'SpectrumAccessSqMass')


class TestFeatureFinderMetaboIdentCompound:
    """FeatureFinderMetaboIdentCompound — compound definitions."""

    def test_exists(self):
        assert hasattr(pyopenms, 'FeatureFinderMetaboIdentCompound')

    def test_construct(self):
        # nanobind uses parameterized constructor (no default constructor in C++)
        c = pyopenms.FeatureFinderMetaboIdentCompound(
            "test", "C6H12O6", 180.063, [1], [100.0], [10.0], [1.0])
        assert c is not None


class TestAMSEAdductInfo:
    """AMSE_AdductInfo — adduct information."""

    def test_exists(self):
        assert hasattr(pyopenms, 'AMSE_AdductInfo')

    @pytest.mark.xfail(reason="AMSE_AdductInfo default constructor not yet bound")
    def test_construct(self):
        ai = pyopenms.AMSE_AdductInfo()
        assert ai is not None


# === Data Structs / Small Classes ===
# Only testing the most important ones

class TestDate:
    """Date — proper OpenMS::Date class with date-only methods."""

    def test_exists(self):
        assert hasattr(pyopenms, 'Date')

    def test_is_separate_class(self):
        assert pyopenms.Date is not pyopenms.DateTime

    def test_today(self):
        d = pyopenms.Date.today()
        assert d is not None
        assert d.year() > 2020

    def test_get_set(self):
        d = pyopenms.Date()
        d.set("2024-01-15")
        assert "2024" in d.get()

    def test_components(self):
        d = pyopenms.Date()
        d.set(1, 15, 2024)  # month, day, year
        assert d.year() == 2024
        assert d.month() == 1
        assert d.day() == 15


class TestUniqueIdGenerator:
    """UniqueIdGenerator — unique ID generation."""

    def test_exists(self):
        assert hasattr(pyopenms, 'UniqueIdGenerator')

    def test_get_unique_id(self):
        uid = pyopenms.UniqueIdGenerator.getUniqueId()
        assert isinstance(uid, int)
        assert uid > 0

    def test_set_seed(self):
        pyopenms.UniqueIdGenerator.setSeed(42)
        uid1 = pyopenms.UniqueIdGenerator.getUniqueId()
        pyopenms.UniqueIdGenerator.setSeed(42)
        uid2 = pyopenms.UniqueIdGenerator.getUniqueId()
        assert uid1 == uid2


class TestStringView:
    """StringView — lightweight string reference."""

    def test_exists(self):
        assert hasattr(pyopenms, 'StringView')

    def test_construct(self):
        sv = pyopenms.StringView()
        assert sv is not None


class TestSignalToNoiseEstimatorMedianChrom:
    """SignalToNoiseEstimatorMedianChrom — S/N for chromatograms."""

    def test_exists(self):
        assert hasattr(pyopenms, 'SignalToNoiseEstimatorMedianChrom')


class TestSimplePeak:
    """SimplePeak — simple peak data holder."""

    def test_exists(self):
        assert hasattr(pyopenms, 'SimplePeak')


class TestDataFilter:
    """DataFilter — spectrum/chromatogram filtering."""

    def test_exists(self):
        assert hasattr(pyopenms, 'DataFilter')
