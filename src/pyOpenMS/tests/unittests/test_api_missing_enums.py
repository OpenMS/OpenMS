"""
Tests for missing enum API completeness (Task 4).

Verifies that enums from the Cython 3.5 API are accessible
at module level in the nanobind build.
"""
import pytest
import pyopenms


# === 8 Completely Missing Enums ===

class TestDRangeIntersection:
    """DRangeIntersection — returned by DRange intersects()."""

    def test_exists(self):
        assert hasattr(pyopenms, 'DRangeIntersection')

    def test_values(self):
        assert hasattr(pyopenms.DRangeIntersection, 'Disjoint')
        assert hasattr(pyopenms.DRangeIntersection, 'Intersects')
        assert hasattr(pyopenms.DRangeIntersection, 'Inside')


class TestMTQuantMethod:
    """MT_QUANTMETHOD — MassTrace quantification method."""

    def test_exists(self):
        assert hasattr(pyopenms, 'MT_QUANTMETHOD')

    def test_values(self):
        assert hasattr(pyopenms.MT_QUANTMETHOD, 'MT_QUANT_AREA')
        assert hasattr(pyopenms.MT_QUANTMETHOD, 'MT_QUANT_MEDIAN')
        assert hasattr(pyopenms.MT_QUANTMETHOD, 'MT_QUANT_HEIGHT')


class TestMZTrafoModelModelType:
    """MZTrafoModel_MODELTYPE — m/z transformation model types."""

    def test_exists(self):
        assert hasattr(pyopenms, 'MZTrafoModel_MODELTYPE')

    def test_values(self):
        mt = pyopenms.MZTrafoModel_MODELTYPE
        assert hasattr(mt, 'LINEAR')
        assert hasattr(mt, 'LINEAR_WEIGHTED')
        assert hasattr(mt, 'QUADRATIC')
        assert hasattr(mt, 'QUADRATIC_WEIGHTED')


class TestQuotingMethod:
    """QuotingMethod — string quoting methods for CSV."""

    def test_exists(self):
        assert hasattr(pyopenms, 'QuotingMethod')

    def test_values(self):
        qm = pyopenms.QuotingMethod
        assert hasattr(qm, 'NONE')
        assert hasattr(qm, 'ESCAPE')
        assert hasattr(qm, 'DOUBLE')


class TestSIDE:
    """SIDE — used in adduct decomposition."""

    def test_exists(self):
        assert hasattr(pyopenms, 'SIDE')


class TestUnitType:
    """UnitType — DataValue unit ontology type."""

    def test_exists(self):
        assert hasattr(pyopenms, 'UnitType')

    def test_values(self):
        ut = pyopenms.UnitType
        assert hasattr(ut, 'UNIT_ONTOLOGY')
        assert hasattr(ut, 'MS_ONTOLOGY')
        assert hasattr(ut, 'OTHER')


class TestValueType:
    """ValueType — ParamValue/DataValue type tags."""

    def test_exists(self):
        assert hasattr(pyopenms, 'ValueType')

    def test_values(self):
        vt = pyopenms.ValueType
        assert hasattr(vt, 'STRING_VALUE')
        assert hasattr(vt, 'INT_VALUE')
        assert hasattr(vt, 'DOUBLE_VALUE')
        assert hasattr(vt, 'EMPTY_VALUE')


class TestXRefTypeCVTermControlledVocabulary:
    """XRefType_CVTerm_ControlledVocabulary — CV term cross-reference types."""

    def test_exists(self):
        assert hasattr(pyopenms, 'XRefType_CVTerm_ControlledVocabulary')


# === 8 Nested Enums needing module-level aliases ===
# These exist as nested enums but should also be accessible at module level
# for backward compatibility with Cython 3.5 API

class TestAnnotationStateModuleLevel:
    """AnnotationState should be accessible at module level."""

    def test_module_level(self):
        assert hasattr(pyopenms, 'AnnotationState')

    def test_nested_still_works(self):
        """Should also be accessible via BaseFeature."""
        assert hasattr(pyopenms.BaseFeature, 'AnnotationState')

    def test_values(self):
        # Access via whichever path works
        cls = getattr(pyopenms, 'AnnotationState',
                      getattr(pyopenms.BaseFeature, 'AnnotationState', None))
        assert cls is not None
        assert hasattr(cls, 'FEATURE_ID_NONE')


class TestBoundaryConditionModuleLevel:
    """BoundaryCondition should be accessible at module level."""

    def test_module_level(self):
        assert hasattr(pyopenms, 'BoundaryCondition')


class TestChecksumTypeModuleLevel:
    """ChecksumType should be accessible at module level."""

    def test_module_level(self):
        assert hasattr(pyopenms, 'ChecksumType')

    def test_values(self):
        cls = getattr(pyopenms, 'ChecksumType',
                      getattr(pyopenms.SourceFile, 'ChecksumType', None))
        assert cls is not None
        assert hasattr(cls, 'SHA1')


class TestDecoyTransitionTypeModuleLevel:
    """DecoyTransitionType should be accessible at module level."""

    def test_module_level(self):
        assert hasattr(pyopenms, 'DecoyTransitionType')


class TestDimensionDescriptionModuleLevel:
    """DimensionDescription should be accessible at module level."""

    def test_module_level(self):
        assert hasattr(pyopenms, 'DimensionDescription')


class TestMeasureModuleLevel:
    """Measure should be accessible at module level."""

    def test_module_level(self):
        assert hasattr(pyopenms, 'Measure')


class TestNormalizationMethodModuleLevel:
    """NormalizationMethod should be accessible at module level."""

    def test_module_level(self):
        assert hasattr(pyopenms, 'NormalizationMethod')


class TestOriginAnnotationFormatModuleLevel:
    """OriginAnnotationFormat should be accessible at module level."""

    def test_module_level(self):
        assert hasattr(pyopenms, 'OriginAnnotationFormat')


# === 5 Already at module level (verify they work like the Cython version) ===

class TestExistingEnumsBackwardCompat:
    """Enums that exist at module level — verify int conversion works."""

    @pytest.mark.parametrize("enum_name,value_name", [
        ("DriftTimeUnit", "MILLISECOND"),
        ("FileType", "MZML"),
        ("IMFormat", "NONE"),
        ("IMPeakType", "IM_PROFILE"),
        ("LogType", "CMD"),
    ])
    def test_enum_value_accessible(self, enum_name, value_name):
        enum_cls = getattr(pyopenms, enum_name)
        val = getattr(enum_cls, value_name)
        assert val is not None

    @pytest.mark.parametrize("enum_name,value_name", [
        ("DriftTimeUnit", "MILLISECOND"),
        ("FileType", "MZML"),
    ])
    def test_enum_int_conversion(self, enum_name, value_name):
        """Enums should be convertible to int (via int() or .value)."""
        enum_cls = getattr(pyopenms, enum_name)
        val = getattr(enum_cls, value_name)
        # nanobind enums use .value; Cython enums support int()
        try:
            int_val = int(val)
        except TypeError:
            int_val = val.value
        assert isinstance(int_val, int)
