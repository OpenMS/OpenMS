"""
Tests for the pyOpenMS generator components.

These tests verify that the generator correctly parses .pxd files
and emits nanobind C++ code.
"""

import pytest
import textwrap


class TestTypeRegistry:
    """Tests for TypeRegistry."""

    def test_primitive_types(self):
        """Test resolution of primitive types."""
        from generator.type_registry import TypeRegistry, TypeCategory

        registry = TypeRegistry()

        # Test integer types
        info = registry.resolve_type("int")
        assert info.cpp_type == "int"
        assert info.python_type == "int"
        assert info.category == TypeCategory.PRIMITIVE

        # Test floating types
        info = registry.resolve_type("double")
        assert info.cpp_type == "double"
        assert info.python_type == "float"

        # Test boolean
        info = registry.resolve_type("bool")
        assert info.cpp_type == "bool"
        assert info.python_type == "bool"

    def test_string_type(self):
        """Test resolution of OpenMS::String type."""
        from generator.type_registry import TypeRegistry, TypeCategory

        registry = TypeRegistry()

        info = registry.resolve_type("String")
        assert info.cpp_type == "OpenMS::String"
        assert info.python_type == "str"
        assert info.needs_custom_caster

    def test_container_types(self):
        """Test resolution of container types."""
        from generator.type_registry import TypeRegistry

        registry = TypeRegistry()

        # Vector
        info = registry.resolve_type("libcpp_vector[int]")
        assert "std::vector<int>" in info.cpp_type
        assert "List[int]" in info.python_type

        # Map
        info = registry.resolve_type("libcpp_map[String, int]")
        assert "std::map" in info.cpp_type

    def test_dposition_types(self):
        """Test resolution of DPosition types."""
        from generator.type_registry import TypeRegistry

        registry = TypeRegistry()

        info = registry.resolve_type("DPosition1")
        assert "DPosition<1>" in info.cpp_type
        assert info.needs_custom_caster

        info = registry.resolve_type("DPosition2")
        assert "DPosition<2>" in info.cpp_type

    def test_class_registration(self):
        """Test class registration and lookup."""
        from generator.type_registry import TypeRegistry, TypeCategory

        registry = TypeRegistry()
        registry.register_class("MSSpectrum", "OpenMS")

        assert registry.is_known_class("MSSpectrum")
        assert registry.is_known_class("OpenMS::MSSpectrum")

        info = registry.resolve_type("MSSpectrum")
        assert info.category == TypeCategory.CLASS
        assert info.cpp_type == "OpenMS::MSSpectrum"

    def test_type_qualifiers(self):
        """Test handling of const, reference, pointer qualifiers."""
        from generator.type_registry import TypeRegistry

        registry = TypeRegistry()

        # Const
        info = registry.resolve_type("const int")
        assert info.is_const

        # Reference
        info = registry.resolve_type("int&")
        assert info.is_reference

        # Pointer
        info = registry.resolve_type("int*")
        assert info.is_pointer


class TestPxdParser:
    """Tests for PxdParser."""

    @pytest.fixture
    def sample_pxd(self, tmp_path):
        """Create a sample .pxd file for testing."""
        pxd_content = textwrap.dedent('''
            from libcpp cimport *
            from Types cimport *

            cdef extern from "<OpenMS/KERNEL/Peak1D.h>" namespace "OpenMS":

                cdef cppclass Peak1D:
                    # wrap-hash:
                    #  std
                    # wrap-doc:
                    #  A 1-dimensional raw data point or peak.

                    Peak1D() except + nogil
                    Peak1D(Peak1D &) except + nogil

                    float getIntensity() nogil
                    double getMZ() nogil
                    void setMZ(double) nogil
                    void setIntensity(float) nogil
                    bool operator==(Peak1D) except + nogil
        ''')

        pxd_file = tmp_path / "Peak1D.pxd"
        pxd_file.write_text(pxd_content)
        return pxd_file

    def test_parse_class(self, sample_pxd):
        """Test parsing of class declaration."""
        from generator.pxd_parser import PxdParser
        from generator.type_registry import TypeRegistry

        registry = TypeRegistry()
        parser = PxdParser(registry)

        classes = parser.parse_file(sample_pxd)

        assert "Peak1D" in classes
        cls = classes["Peak1D"]
        assert cls.name == "Peak1D"
        assert cls.namespace == "OpenMS"
        assert cls.header_file == "<OpenMS/KERNEL/Peak1D.h>"

    def test_parse_constructors(self, sample_pxd):
        """Test parsing of constructors."""
        from generator.pxd_parser import PxdParser
        from generator.type_registry import TypeRegistry

        registry = TypeRegistry()
        parser = PxdParser(registry)

        classes = parser.parse_file(sample_pxd)
        cls = classes["Peak1D"]

        # Should have default and copy constructors
        assert len(cls.constructors) >= 1

    def test_parse_methods(self, sample_pxd):
        """Test parsing of methods."""
        from generator.pxd_parser import PxdParser
        from generator.type_registry import TypeRegistry

        registry = TypeRegistry()
        parser = PxdParser(registry)

        classes = parser.parse_file(sample_pxd)
        cls = classes["Peak1D"]

        method_names = [m.name for m in cls.methods]
        assert "getIntensity" in method_names
        assert "getMZ" in method_names
        assert "setMZ" in method_names

    def test_parse_wrap_directives(self, sample_pxd):
        """Test parsing of wrap- directives."""
        from generator.pxd_parser import PxdParser
        from generator.type_registry import TypeRegistry

        registry = TypeRegistry()
        parser = PxdParser(registry)

        classes = parser.parse_file(sample_pxd)
        cls = classes["Peak1D"]

        assert cls.wrap_hash is not None
        assert "1-dimensional" in cls.doc


class TestNanobindEmitter:
    """Tests for NanobindEmitter."""

    def test_emitter_creation(self):
        """Test that emitter can be created."""
        from generator.nanobind_emitter import NanobindEmitter

        emitter = NanobindEmitter()
        assert emitter.single_module is False

    def test_caster_owned_types_detection(self):
        """Test auto-detection of caster-owned types."""
        from generator.nanobind_emitter import get_caster_owned_types

        types = get_caster_owned_types()
        # Should detect types from type_casters/*.h
        assert "String" in types
        assert "DataValue" in types

    def test_container_detection(self):
        """Test AST-based container detection."""
        from generator.nanobind_emitter import NanobindEmitter
        from generator.cpp_parser import CppClass, CppMethod, MergedClass, MergedMethod

        emitter = NanobindEmitter()

        # Create a mock container class
        cpp_class = CppClass(
            name="TestContainer",
            qualified_name="OpenMS::TestContainer",
            namespace="OpenMS",
            header_file="test.h",
        )
        cpp_class.methods = [
            CppMethod(name="size", return_type="size_t"),
            CppMethod(name="begin", return_type="iterator"),
            CppMethod(name="end", return_type="iterator"),
        ]
        merged = MergedClass(
            cpp_class=cpp_class,
            methods=[MergedMethod(cpp_method=m) for m in cpp_class.methods],
            constructors=[],
        )

        assert emitter._has_size_method(merged)
        assert emitter._has_iterator_methods(merged)
