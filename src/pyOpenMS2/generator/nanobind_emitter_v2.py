"""
Nanobind C++ Code Emitter for pyOpenMS2 (v2 - using libclang info)

This module generates nanobind C++ binding code from merged C++/pxd information.
It uses accurate C++ type information from libclang for correct method signatures.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .cpp_parser import CppMethod, CppParameter, MergedClass, MergedMethod


# Methods to skip for specific classes (complex return types or signature mismatches)
SKIP_METHODS = {
    "MSSpectrum": {
        "getFloatDataArrays", "setFloatDataArrays",
        "getIntegerDataArrays", "setIntegerDataArrays",
        "getStringDataArrays", "setStringDataArrays",
        "getIMData",  # Returns pair
        "calculateTIC",  # Return type mismatch (double in .pxd vs float in C++)
        "findNearest",  # Return type mismatch (Size vs int in .pxd)
        "findHighestInWindow",  # Return type mismatch
        "getName",  # Returns const String& vs String in .pxd
        "getType",  # Complex return type
        "select",  # Complex parameter types
        "getSourceFile",  # const/non-const overload
    },
    "MSChromatogram": {
        "getFloatDataArrays", "setFloatDataArrays",
        "getIntegerDataArrays", "setIntegerDataArrays",
        "getStringDataArrays", "setStringDataArrays",
        "getSourceFile",  # const/non-const overload
    },
    "MSExperiment": {
        "getSpectrum", "getSpectra",  # const/non-const overloads
        "getChromatogram", "getChromatograms",  # const/non-const overloads
        "getSourceFiles",  # const/non-const overloads
        "setSourceFile",  # const/non-const overloads
    },
}

# Classes to skip due to incomplete type dependencies or other issues
# These will need manual handling or fixes in the C++ headers
SKIP_CLASSES = {
    # Type caster conflict - have casters for these so can't bind as class
    "String",               # Has type caster (OpenMS::String <-> str)
    "DataValue",            # Has type caster (OpenMS::DataValue <-> Python types)

    # Incomplete type issues
    "MassExplainer",        # References Compomer which is forward-declared
    "ILPDCWrapper",         # Complex ILP dependencies
    "SwathFileConsumer",    # Template complexity
    "MSDataWritingConsumer", # Template complexity
    "SpectrumAccessOpenMSInMemory",  # Complex template
    # Abstract classes that libclang may not detect
    "BaseGroupFinder",
    "BaseSuperimposer",
    "ConsensusIDAlgorithm",
    "ConsensusIDAlgorithmIdentity",
    # Nested classes or special cases
    "XMLHandler",
    "IndexedMzMLHandler",
    # Classes with static methods not marked in .pxd (need @staticmethod or wrap-static)
    "PercolatorFeatureSetHelper",
    "PercolatorInfile",
    "PercolatorOutfile",
    "TransformationXMLFile",
    "OpenSwathDataAccessHelper",
    # Classes with no default constructor
    "QTCluster",
    "SpectrumAccessOpenMS",
    # Classes with Cython template syntax issues
    "TraMLFile",
    "MZTrafoModel",
    # Private inheritance issues (inherit from std::vector but not accessible)
    "AcquisitionInfo",

    # Classes with complex constructors that reference unbound types
    "TransformationModelBSpline",
    "TransformationModelLinear",
    "TransformationModelLowess",
    "TransformationModelInterpolated",

    # Classes with copy constructors that cause overload ambiguity
    "CVMappings",
    "ConsensusMapNormalizerAlgorithmMedian",
    "ConsensusMapNormalizerAlgorithmQuantile",

    # Classes where the type is actually a template alias (DPosition2 -> DPosition<2>)
    "DPosition2",

    # Classes with overloaded methods that can't be resolved without libclang
    "Biosaur2Algorithm",  # setMSData has const ref and move versions
    "ConsensusMap",       # push_back has multiple overloads
    "Feature",            # getSubordinates has const/non-const overloads
    "FeatureMap",         # push_back has multiple overloads
    "FeatureFinderAlgorithmMetaboIdent",  # setMSData overloads
    "FeatureFinderIdentificationAlgorithm",  # constructor parameter issues
    "Instrument",         # getIonSources has const/non-const overloads
    "ExperimentalDesign", # getters with overloads
    "ExperimentalSettings",  # similar pattern
    "CVMappingRule",      # copy constructor issues

    # Classes with template parameters or complex constructors
    "OpenSwathWorkflowSonar",
    "OpenSwathWorkflow",
    "MRMTransitionGroup",
    "MRMTransitionGroupCP",
    "LightMRMTransitionGroup",

    # Classes with private copy constructors
    "CVMappingFile",

    # Classes with unresolved overloads
    "KroenikFile",
    "MascotGenericFile",
    "MzTabFile",
    "NLargest",
    "RNaseDigestion",
    "RankScaler",
    "PeakGroup",
    "Mobilogram",

    # Classes with type alias issues
    "OSW_ChromExtractParams",

    # TODO: Add wrap-static support in .pxd files and remove from SKIP_CLASSES
}

# Additional headers needed for specific classes
ADDITIONAL_INCLUDES = {
    "MSSpectrum": ["<OpenMS/KERNEL/Peak1D.h>"],
    "MSChromatogram": ["<OpenMS/KERNEL/ChromatogramPeak.h>"],
    "MSExperiment": ["<OpenMS/KERNEL/MSSpectrum.h>", "<OpenMS/KERNEL/MSChromatogram.h>"],
    "Feature": ["<OpenMS/KERNEL/Feature.h>"],
    "FeatureMap": ["<OpenMS/KERNEL/Feature.h>"],
    "ConsensusMap": ["<OpenMS/KERNEL/ConsensusFeature.h>"],
}

# Classes that need special __len__ support (container-like)
CONTAINER_CLASSES = {
    "MSSpectrum", "MSChromatogram", "MSExperiment", "FeatureMap",
    "ConsensusMap", "PeakMap", "Mobilogram",
    "FloatDataArray", "IntegerDataArray", "StringDataArray",
}

# Classes that need iterator support
ITERABLE_CLASSES = {
    "MSSpectrum", "MSChromatogram", "MSExperiment", "FeatureMap",
    "ConsensusMap", "PeakMap", "Mobilogram", "AASequence",
}

# Core classes to bind first - these are well-tested and have simple APIs
# All other classes require more work on type normalization
CORE_CLASSES = {
    # Basic peak types
    "Peak1D", "Peak2D", "ChromatogramPeak", "MobilityPeak1D",
    # Spectrum and chromatogram
    "MSSpectrum", "MSChromatogram",
    # Basic data types (DataValue has type caster, excluded)
    "DateTime",
    # Note: Normalizer excluded - filterSpectrum is a template method not compatible with Doxygen parsing
}

# Classes that inherit from std::vector (have size() method)
VECTOR_BASED_CLASSES = {
    "MSSpectrum", "MSChromatogram", "Mobilogram",
    "FloatDataArray", "IntegerDataArray", "StringDataArray",
}

# Methods to skip for vector-based classes (inherited from private std::vector base)
# We provide lambda implementations instead via CONTAINER_CLASSES/__len__
VECTOR_INHERITED_METHODS = {
    "size", "reserve", "resize", "push_back", "pop_back", "clear",
    "empty", "front", "back", "begin", "end", "cbegin", "cend",
    "at", "operator[]", "data", "capacity", "max_size", "shrink_to_fit",
    "insert", "erase", "emplace", "emplace_back", "assign", "swap",
}

# Special method implementations for critical classes
# These are C++ lambdas that implement Python methods
SPECIAL_METHODS = {
    "MSSpectrum": {
        "get_peaks": '''
        .def("get_peaks", [](const OpenMS::MSSpectrum& self) {
            // Return (mz_array, intensity_array) as vectors (converted to numpy by nanobind)
            const size_t n = self.size();
            std::vector<double> mz(n);
            std::vector<double> intensity(n);
            for (size_t i = 0; i < n; ++i) {
                mz[i] = self[i].getMZ();
                intensity[i] = self[i].getIntensity();
            }
            return nb::make_tuple(mz, intensity);
        }, "Returns a tuple of (mz_array, intensity_array)")''',
        "set_peaks": '''
        .def("set_peaks", [](OpenMS::MSSpectrum& self, std::vector<double> mz, std::vector<double> intensity) {
            // Accept mz_array and intensity_array separately
            const size_t n = mz.size();
            if (intensity.size() != n) {
                throw std::runtime_error("mz and intensity arrays must have same length");
            }
            self.resize(n);
            for (size_t i = 0; i < n; ++i) {
                self[i].setMZ(mz[i]);
                self[i].setIntensity(intensity[i]);
            }
        }, "mz"_a, "intensity"_a, "Set peaks from mz and intensity arrays")''',
    },
    "MSChromatogram": {
        "get_peaks": '''
        .def("get_peaks", [](const OpenMS::MSChromatogram& self) {
            const size_t n = self.size();
            std::vector<double> rt(n);
            std::vector<double> intensity(n);
            for (size_t i = 0; i < n; ++i) {
                rt[i] = self[i].getRT();
                intensity[i] = self[i].getIntensity();
            }
            return nb::make_tuple(rt, intensity);
        }, "Returns a tuple of (rt_array, intensity_array)")''',
    },
    "FloatDataArray": {
        "get_data": '''
        .def("get_data", [](const OpenMS::DataArrays::FloatDataArray& self) {
            std::vector<float> arr(self.begin(), self.end());
            return arr;
        }, "Returns a copy of the data as a list")''',
        "set_data": '''
        .def("set_data", [](OpenMS::DataArrays::FloatDataArray& self, std::vector<float> arr) {
            self.assign(arr.begin(), arr.end());
        }, "data"_a, "Set data from a list")''',
    },
    "IntegerDataArray": {
        "get_data": '''
        .def("get_data", [](const OpenMS::DataArrays::IntegerDataArray& self) {
            std::vector<OpenMS::Int> arr(self.begin(), self.end());
            return arr;
        }, "Returns a copy of the data as a list")''',
        "set_data": '''
        .def("set_data", [](OpenMS::DataArrays::IntegerDataArray& self, std::vector<OpenMS::Int> arr) {
            self.assign(arr.begin(), arr.end());
        }, "data"_a, "Set data from a list")''',
    },
    "MSExperiment": {
        "get_spectra": '''
        .def("get_spectra", [](OpenMS::MSExperiment& self) -> std::vector<OpenMS::MSSpectrum>& {
            return self.getSpectra();
        }, nb::rv_policy::reference_internal, "Returns reference to spectra vector")''',
        "get_chromatograms": '''
        .def("get_chromatograms", [](OpenMS::MSExperiment& self) -> std::vector<OpenMS::MSChromatogram>& {
            return self.getChromatograms();
        }, nb::rv_policy::reference_internal, "Returns reference to chromatograms vector")''',
    },
}


@dataclass
class ModuleContent:
    """Content for a single C++ module file."""
    includes: Set[str]
    forward_declarations: List[str]
    class_bindings: List[str]
    enum_bindings: List[str]


class NanobindEmitterV2:
    """
    Generates nanobind C++ binding code from merged C++/pxd declarations.
    Uses accurate C++ type information from libclang.
    """

    def __init__(self, num_modules: int = 8, core_only: bool = True):
        self.num_modules = num_modules
        self.core_only = core_only  # Only bind CORE_CLASSES for now

        # Standard includes for all modules
        self._standard_includes = {
            "<nanobind/nanobind.h>",
            "<nanobind/operators.h>",
            "<nanobind/stl/string.h>",
            "<nanobind/stl/vector.h>",
            "<nanobind/stl/map.h>",
            "<nanobind/stl/set.h>",
            "<nanobind/stl/shared_ptr.h>",
            "<nanobind/stl/optional.h>",
            "<nanobind/ndarray.h>",
            "<nanobind/make_iterator.h>",
            "<sstream>",       # For std::ostringstream in __repr__
            "<iomanip>",       # For std::setprecision in __repr__
            '"all_casters.h"',
        }

        # Operator name mappings
        self._operator_map = {
            "operator==": "__eq__",
            "operator!=": "__ne__",
            "operator<": "__lt__",
            "operator<=": "__le__",
            "operator>": "__gt__",
            "operator>=": "__ge__",
            "operator+": "__add__",
            "operator-": "__sub__",
            "operator*": "__mul__",
            "operator/": "__truediv__",
            "operator+=": "__iadd__",
            "operator-=": "__isub__",
            "operator*=": "__imul__",
            "operator/=": "__itruediv__",
            "operator[]": "__getitem__",
            "operator()": "__call__",
        }

    def emit(
        self,
        classes: Dict[str, MergedClass],
        output_dir: Path,
    ) -> None:
        """
        Generate nanobind C++ files.

        Parameters
        ----------
        classes : dict
            Merged class declarations with accurate C++ info.
        output_dir : Path
            Output directory for generated files.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Distribute classes across modules
        modules = self._distribute_classes(classes)

        # Generate module files
        for i, module_classes in enumerate(modules, 1):
            content = self._generate_module_content(module_classes, classes)
            output_file = output_dir / f"module_{i}.cpp"
            self._write_module_file(output_file, i, content)

        # Generate main module file
        self._write_main_module(output_dir / "main_module.cpp", len(modules))

    def _distribute_classes(
        self, classes: Dict[str, MergedClass]
    ) -> List[List[MergedClass]]:
        """Distribute classes across modules for parallel compilation."""
        sorted_classes = sorted(classes.values(), key=lambda c: c.name)

        modules: List[List[MergedClass]] = [[] for _ in range(self.num_modules)]

        for merged_class in sorted_classes:
            hash_val = int(hashlib.md5(merged_class.name.encode()).hexdigest(), 16)
            module_idx = hash_val % self.num_modules
            modules[module_idx].append(merged_class)

        return modules

    def _generate_module_content(
        self,
        module_classes: List[MergedClass],
        all_classes: Dict[str, MergedClass],
    ) -> ModuleContent:
        """Generate content for a module."""
        content = ModuleContent(
            includes=self._standard_includes.copy(),
            forward_declarations=[],
            class_bindings=[],
            enum_bindings=[],
        )

        for merged_class in module_classes:
            class_name = merged_class.name

            # Check if this class will actually be bound (avoid including unnecessary headers)
            # In core_only mode, only CORE_CLASSES are bound
            if self.core_only and class_name not in CORE_CLASSES:
                continue

            # Skip problematic classes
            if class_name in SKIP_CLASSES:
                continue

            # Skip nested classes
            qualified_name = merged_class.qualified_name
            if "::" in qualified_name.replace("OpenMS::", "", 1):
                continue

            # Skip abstract classes
            if merged_class.is_abstract:
                continue

            # Generate class binding first - only add headers if successful
            binding = self._generate_class_binding(merged_class)
            if binding:
                content.class_bindings.append(binding)

                # Add OpenMS header - extract relative path from full path
                header_file = merged_class.header_file
                if header_file:
                    # Extract OpenMS-relative path (e.g., OpenMS/KERNEL/Peak1D.h)
                    # Look for the last occurrence of "OpenMS/" to handle paths like
                    # /path/to/OpenMS/src/openms/include/OpenMS/KERNEL/Peak1D.h
                    if "OpenMS/" in header_file:
                        idx = header_file.rfind("/include/")
                        if idx != -1:
                            header = header_file[idx + len("/include/"):]
                        else:
                            idx = header_file.rfind("OpenMS/")
                            header = header_file[idx:]
                    else:
                        header = header_file
                    content.includes.add(f"<{header}>")

                # Add additional includes for this class if needed
                if class_name in ADDITIONAL_INCLUDES:
                    for inc in ADDITIONAL_INCLUDES[class_name]:
                        content.includes.add(inc)

        return content

    def _generate_class_binding(self, merged_class: MergedClass) -> Optional[str]:
        """Generate nanobind binding code for a class."""
        lines = []

        class_name = merged_class.name
        qualified_name = merged_class.qualified_name

        # In core_only mode, only bind CORE_CLASSES
        if self.core_only and class_name not in CORE_CLASSES:
            return None

        # Skip problematic classes
        if class_name in SKIP_CLASSES:
            return None

        # Skip nested classes
        if "::" in qualified_name.replace("OpenMS::", "", 1):
            return None

        # Skip abstract classes (can't instantiate)
        if merged_class.is_abstract:
            return None

        # Handle base classes - nanobind supports single inheritance
        # Note: We can only specify base classes that are also bound to nanobind.
        # Specifying an unbound base class causes a runtime error.
        # For now, we skip specifying base classes since RangeManagerContainer,
        # MetaInfoInterface etc. are not in CORE_CLASSES.
        base_spec = ""
        # TODO: When we bind more classes, we can enable base class specification
        # for bases that are actually bound. For now, skip to avoid runtime errors.

        # Start class definition
        lines.append(f'    nb::class_<{qualified_name}{base_spec}>(m, "{class_name}")')

        # Add docstring
        if merged_class.doc:
            doc = self._escape_string(merged_class.doc)
            lines[-1] = lines[-1][:-1]  # Remove trailing paren
            lines[-1] += f', "{doc}")'

        # Generate constructors
        for ctor in merged_class.constructors:
            ctor_code = self._generate_constructor(ctor)
            if ctor_code:
                lines.append(f"        {ctor_code}")

        # Generate methods
        for merged_method in merged_class.methods:
            if merged_method.wrap_ignore:
                continue
            method_code = self._generate_method(merged_method, qualified_name, merged_class)
            if method_code:
                lines.append(f"        {method_code}")

        # Add hash support
        if merged_class.wrap_hash:
            lines.append(f'        .def("__hash__", [](const {qualified_name}& self) {{ return std::hash<{qualified_name}>{{}}(self); }})')

        # Add iterator support (from wrap-iter directive or known iterable classes)
        if merged_class.wrap_iter or class_name in ITERABLE_CLASSES:
            lines.append(f'        .def("__iter__", []({qualified_name}& self) {{ return nb::make_iterator<nb::rv_policy::reference_internal>(nb::handle(), "{class_name}_iter", self.begin(), self.end()); }})')

        # Add __len__ for container classes
        if class_name in CONTAINER_CLASSES or class_name in VECTOR_BASED_CLASSES:
            lines.append(f'        .def("__len__", []({qualified_name}& self) {{ return self.size(); }})')

        # Add __getitem__ for vector-based classes
        if class_name in VECTOR_BASED_CLASSES:
            # Determine element type from class name
            elem_type = self._get_element_type(class_name)
            if elem_type:
                lines.append(f'        .def("__getitem__", []({qualified_name}& self, size_t i) -> {elem_type}& {{ ')
                lines.append(f'            if (i >= self.size()) throw nb::index_error();')
                lines.append(f'            return self[i];')
                lines.append(f'        }}, nb::rv_policy::reference_internal)')

        # Add special methods for this class (get_peaks, set_peaks, etc.)
        if class_name in SPECIAL_METHODS:
            for method_name, method_code in SPECIAL_METHODS[class_name].items():
                lines.append(method_code)

        # Add __repr__ for key classes
        repr_code = self._generate_repr(class_name, qualified_name)
        if repr_code:
            lines.append(repr_code)

        # Close class definition
        lines.append("        ;")

        return "\n".join(lines)

    def _get_element_type(self, class_name: str) -> Optional[str]:
        """Get the element type for a vector-based class."""
        element_types = {
            "MSSpectrum": "OpenMS::Peak1D",
            "MSChromatogram": "OpenMS::ChromatogramPeak",
            "Mobilogram": "OpenMS::MobilityPeak1D",
            "FloatDataArray": "float",
            "IntegerDataArray": "OpenMS::Int",
            "StringDataArray": "OpenMS::String",
        }
        return element_types.get(class_name)

    def _generate_repr(self, class_name: str, qualified_name: str) -> Optional[str]:
        """Generate __repr__ method for a class."""
        # Define repr formats for key classes
        repr_formats = {
            "Peak1D": f'''        .def("__repr__", []({qualified_name}& self) {{
            return "<Peak1D mz=" + std::to_string(self.getMZ()) + " intensity=" + std::to_string(self.getIntensity()) + ">";
        }})''',
            "Peak2D": f'''        .def("__repr__", []({qualified_name}& self) {{
            return "<Peak2D rt=" + std::to_string(self.getRT()) + " mz=" + std::to_string(self.getMZ()) + " intensity=" + std::to_string(self.getIntensity()) + ">";
        }})''',
            "ChromatogramPeak": f'''        .def("__repr__", []({qualified_name}& self) {{
            return "<ChromatogramPeak rt=" + std::to_string(self.getRT()) + " intensity=" + std::to_string(self.getIntensity()) + ">";
        }})''',
            "MSSpectrum": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<MSSpectrum ms_level=" << self.getMSLevel()
                << " rt=" << std::fixed << std::setprecision(2) << self.getRT()
                << " num_peaks=" << self.size() << ">";
            return oss.str();
        }})''',
            "MSChromatogram": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<MSChromatogram native_id='" << self.getNativeID()
                << "' num_peaks=" << self.size() << ">";
            return oss.str();
        }})''',
            "MSExperiment": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<MSExperiment num_spectra=" << self.getNrSpectra()
                << " num_chromatograms=" << self.getNrChromatograms() << ">";
            return oss.str();
        }})''',
            "Feature": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<Feature rt=" << std::fixed << std::setprecision(2) << self.getRT()
                << " mz=" << std::setprecision(4) << self.getMZ()
                << " intensity=" << std::setprecision(0) << self.getIntensity()
                << " charge=" << self.getCharge() << ">";
            return oss.str();
        }})''',
            "FeatureMap": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<FeatureMap num_features=" << self.size() << ">";
            return oss.str();
        }})''',
            "ConsensusFeature": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<ConsensusFeature rt=" << std::fixed << std::setprecision(2) << self.getRT()
                << " mz=" << std::setprecision(4) << self.getMZ()
                << " intensity=" << std::setprecision(0) << self.getIntensity()
                << " size=" << self.size() << ">";
            return oss.str();
        }})''',
            "ConsensusMap": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<ConsensusMap num_consensus_features=" << self.size() << ">";
            return oss.str();
        }})''',
        }
        return repr_formats.get(class_name)

    def _generate_constructor(self, ctor: CppMethod) -> Optional[str]:
        """Generate constructor binding."""
        # Skip copy and move constructors
        if len(ctor.parameters) == 1:
            param_type = ctor.parameters[0].type_str
            if "&&" in param_type or (ctor.name in param_type and "&" in param_type):
                return None

        if not ctor.parameters:
            return ".def(nb::init<>())"

        param_types = [self._normalize_type(p.type_str) for p in ctor.parameters]
        types_str = ", ".join(param_types)
        return f".def(nb::init<{types_str}>())"

    def _generate_method(
        self, merged_method: MergedMethod, qualified_name: str, merged_class: MergedClass
    ) -> Optional[str]:
        """Generate method binding."""
        method = merged_method.cpp_method
        method_name = merged_method.wrap_as or method.name
        class_name = merged_class.name

        # Skip methods inherited from private std::vector base
        if class_name in VECTOR_BASED_CLASSES and method.name in VECTOR_INHERITED_METHODS:
            return None

        # Skip specific methods with complex return types
        if class_name in SKIP_METHODS and method.name in SKIP_METHODS[class_name]:
            return None

        # Handle operators
        if method.name.startswith("operator"):
            py_op = self._operator_map.get(method.name)
            if py_op:
                return self._generate_operator_binding(method, qualified_name, py_op)
            return None

        # Skip static methods for now (need different binding approach)
        if method.is_static:
            return self._generate_static_method(method, qualified_name, method_name)

        # Regular methods
        return self._generate_regular_method(merged_method, qualified_name, method_name, merged_class)

    def _generate_regular_method(
        self,
        merged_method: MergedMethod,
        qualified_name: str,
        method_name: str,
        merged_class: MergedClass,
    ) -> Optional[str]:
        """Generate binding for a regular method using lambda wrapper.

        We use lambdas for all methods because:
        1. We can't reliably detect C++ overloads without libclang
        2. Lambdas let us explicitly call the method with correct signature
        3. This avoids all overload resolution issues
        """
        method = merged_method.cpp_method

        # Check if method is overloaded by parameter count/types
        overloads = [m for m in merged_class.methods if m.cpp_method.name == method.name]
        is_overloaded = len(overloads) > 1

        # For const/non-const overloads with same signature, prefer const version
        if is_overloaded:
            # Check if this is a const/non-const pair (same params, different const-ness)
            same_params_overloads = [
                m for m in overloads
                if len(m.cpp_method.parameters) == len(method.parameters)
            ]
            if len(same_params_overloads) == 2:
                const_versions = [m for m in same_params_overloads if m.cpp_method.is_const]
                non_const_versions = [m for m in same_params_overloads if not m.cpp_method.is_const]
                if const_versions and non_const_versions:
                    # This is a const/non-const pair - skip the non-const version
                    if not method.is_const:
                        return None

        # Build lambda parameters
        param_decls = []
        param_names_call = []
        param_names_arg = []

        # C++ reserved keywords that can't be used as parameter names
        cpp_keywords = {
            'int', 'char', 'float', 'double', 'void', 'bool', 'long', 'short',
            'unsigned', 'signed', 'const', 'static', 'class', 'struct', 'enum',
            'union', 'virtual', 'public', 'private', 'protected', 'new', 'delete',
            'return', 'if', 'else', 'for', 'while', 'do', 'switch', 'case',
            'break', 'continue', 'default', 'goto', 'try', 'catch', 'throw',
            'namespace', 'using', 'typedef', 'typename', 'template', 'auto',
            'register', 'extern', 'volatile', 'mutable', 'inline', 'explicit',
            'export', 'true', 'false', 'nullptr', 'this', 'operator', 'sizeof',
            'alignof', 'decltype', 'constexpr', 'noexcept', 'override', 'final',
        }

        for p in method.parameters:
            ptype = self._normalize_type(p.type_str)
            # Check if name is valid (not empty, not "arg*", not a C++ keyword)
            valid_name = (
                p.name
                and not p.name.startswith("arg")
                and p.name not in cpp_keywords
            )
            pname = p.name if valid_name else f"p{len(param_decls)}"
            param_decls.append(f"{ptype} {pname}")
            param_names_call.append(pname)
            if valid_name:
                param_names_arg.append(f'"{p.name}"_a')

        params_decl = ", ".join(param_decls)
        params_call = ", ".join(param_names_call)

        # Determine const qualifier for lambda capture
        if method.is_const:
            self_decl = f"const {qualified_name}& self"
        else:
            self_decl = f"{qualified_name}& self"

        # Build the lambda
        if params_decl:
            lambda_sig = f"[]({self_decl}, {params_decl})"
            call_expr = f"self.{method.name}({params_call})"
        else:
            lambda_sig = f"[]({self_decl})"
            call_expr = f"self.{method.name}()"

        # Build def call
        result = f'.def("{method_name}", {lambda_sig} {{ return {call_expr}; }}'
        if param_names_arg:
            result += f", {', '.join(param_names_arg)}"
        if merged_method.doc:
            doc = self._escape_string(merged_method.doc)
            result += f', "{doc}"'
        result += ")"

        return result

    def _generate_static_method(
        self, method: CppMethod, qualified_name: str, method_name: str
    ) -> Optional[str]:
        """Generate binding for a static method."""
        # Build parameter types
        param_types = [self._normalize_type(p.type_str) for p in method.parameters]
        params_str = ", ".join(param_types)

        method_ptr = f"&{qualified_name}::{method.name}"

        return f'.def_static("{method_name}", {method_ptr})'

    def _generate_operator_binding(
        self, method: CppMethod, qualified_name: str, py_op: str
    ) -> str:
        """Generate binding for an operator."""
        op_name = method.name.replace("operator", "")

        if op_name == "[]":
            return f'.def("__getitem__", []({qualified_name}& self, size_t i) {{ return self[i]; }})'
        elif op_name in ["==", "!=", "<", "<=", ">", ">="]:
            return f'.def(nb::self {op_name} nb::self)'
        elif op_name in ["+", "-", "*", "/"]:
            return f'.def(nb::self {op_name} nb::self)'
        elif op_name in ["+=", "-=", "*=", "/="]:
            return f'.def(nb::self {op_name} nb::self)'

        return ""

    def _normalize_type(self, type_str: str, preserve_reference: bool = False) -> str:
        """Normalize C++ type string for use in bindings.

        Parameters
        ----------
        type_str : str
            The type string to normalize
        preserve_reference : bool
            If True, preserve reference qualifiers for return types in static_cast
        """
        if not type_str:
            return "void"

        result = type_str

        # Convert Cython template syntax to C++: libcpp_vector[T] -> std::vector<T>
        import re
        result = re.sub(r'libcpp_vector\[([^\]]+)\]', r'std::vector<\1>', result)
        result = re.sub(r'libcpp_map\[([^,]+),\s*([^\]]+)\]', r'std::map<\1, \2>', result)
        result = re.sub(r'libcpp_set\[([^\]]+)\]', r'std::set<\1>', result)
        result = re.sub(r'libcpp_pair\[([^,]+),\s*([^\]]+)\]', r'std::pair<\1, \2>', result)
        result = re.sub(r'shared_ptr\[([^\]]+)\]', r'std::shared_ptr<\1>', result)

        # OpenMS type aliases that need namespace
        openms_typedefs = {
            'Int': 'OpenMS::Int',
            'UInt': 'OpenMS::UInt',
            'Size': 'OpenMS::Size',
            'String': 'OpenMS::String',
            'StringList': 'OpenMS::StringList',
            'IntList': 'OpenMS::IntList',
            'DoubleList': 'OpenMS::DoubleList',
            'FloatDataArray': 'OpenMS::DataArrays::FloatDataArray',
            'IntegerDataArray': 'OpenMS::DataArrays::IntegerDataArray',
            'StringDataArray': 'OpenMS::DataArrays::StringDataArray',
            'FloatDataArrays': 'OpenMS::MSSpectrum::FloatDataArrays',
            'IntegerDataArrays': 'OpenMS::MSSpectrum::IntegerDataArrays',
            'StringDataArrays': 'OpenMS::MSSpectrum::StringDataArrays',
            'Param': 'OpenMS::Param',
            'DataValue': 'OpenMS::DataValue',
            'ParamValue': 'OpenMS::ParamValue',
            'DateTime': 'OpenMS::DateTime',
            'CVTerm': 'OpenMS::CVTerm',
            'AASequence': 'OpenMS::AASequence',
            'EmpiricalFormula': 'OpenMS::EmpiricalFormula',
            'ProteinIdentification': 'OpenMS::ProteinIdentification',
            'PeptideIdentification': 'OpenMS::PeptideIdentification',
            'PeptideHit': 'OpenMS::PeptideHit',
            'ProteinHit': 'OpenMS::ProteinHit',
            'FeatureHandle': 'OpenMS::FeatureHandle',
            'Precursor': 'OpenMS::Precursor',
            'Product': 'OpenMS::Product',
            'SourceFile': 'OpenMS::SourceFile',
            'ContactPerson': 'OpenMS::ContactPerson',
            'Sample': 'OpenMS::Sample',
            'Software': 'OpenMS::Software',
            'Instrument': 'OpenMS::Instrument',
            'IonSource': 'OpenMS::IonSource',
            'MassAnalyzer': 'OpenMS::MassAnalyzer',
            'IonDetector': 'OpenMS::IonDetector',
            'HPLC': 'OpenMS::HPLC',
            'Gradient': 'OpenMS::Gradient',
            'Digestion': 'OpenMS::Digestion',
            'Modification': 'OpenMS::Modification',
            'Tagging': 'OpenMS::Tagging',
            'Element': 'OpenMS::Element',
            'Residue': 'OpenMS::Residue',
            'ModificationsDB': 'OpenMS::ModificationsDB',
            'ResidueDB': 'OpenMS::ResidueDB',
            'ElementDB': 'OpenMS::ElementDB',
            'EnzymeDB': 'OpenMS::EnzymeDB',
            'DigestionEnzyme': 'OpenMS::DigestionEnzyme',
            'DigestionEnzymeProtein': 'OpenMS::DigestionEnzymeProtein',
            'DigestionEnzymeRNA': 'OpenMS::DigestionEnzymeRNA',
            'ResidueModification': 'OpenMS::ResidueModification',
            'ModificationDefinition': 'OpenMS::ModificationDefinition',
            'ModificationDefinitionsSet': 'OpenMS::ModificationDefinitionsSet',
            'ModifiedPeptideGenerator': 'OpenMS::ModifiedPeptideGenerator',
            'ProteaseDigestion': 'OpenMS::ProteaseDigestion',
            'RNaseDigestion': 'OpenMS::RNaseDigestion',
            'TheoreticalSpectrumGenerator': 'OpenMS::TheoreticalSpectrumGenerator',
            'PeakPickerHiRes': 'OpenMS::PeakPickerHiRes',
            'Normalizer': 'OpenMS::Normalizer',
            'SpectrumLookup': 'OpenMS::SpectrumLookup',
            'MRMMapping': 'OpenMS::MRMMapping',
            'MRMFeatureFinderScoring': 'OpenMS::MRMFeatureFinderScoring',
            'MRMTransitionGroupPicker': 'OpenMS::MRMTransitionGroupPicker',
            'TransitionTSVFile': 'OpenMS::TransitionTSVFile',
            'TransitionPQPFile': 'OpenMS::TransitionPQPFile',
            'TargetedExperiment': 'OpenMS::TargetedExperiment',
            'LightTargetedExperiment': 'OpenMS::LightTargetedExperiment',
            'ReactionMonitoringTransition': 'OpenMS::ReactionMonitoringTransition',
            'IncludeExcludeTarget': 'OpenMS::IncludeExcludeTarget',
            'TargetedExperimentHelper': 'OpenMS::TargetedExperimentHelper',
            'SwathMap': 'OpenMS::SwathMap',
            'OpenSwathScoring': 'OpenMS::OpenSwathScoring',
            'ChromatogramExtractor': 'OpenMS::ChromatogramExtractor',
            'MzMLFile': 'OpenMS::MzMLFile',
            'MzXMLFile': 'OpenMS::MzXMLFile',
            'MzDataFile': 'OpenMS::MzDataFile',
            'FeatureXMLFile': 'OpenMS::FeatureXMLFile',
            'ConsensusXMLFile': 'OpenMS::ConsensusXMLFile',
            'IdXMLFile': 'OpenMS::IdXMLFile',
            'PepXMLFile': 'OpenMS::PepXMLFile',
            'ProtXMLFile': 'OpenMS::ProtXMLFile',
            'MascotXMLFile': 'OpenMS::MascotXMLFile',
            'FASTAFile': 'OpenMS::FASTAFile',
            'FASTAEntry': 'OpenMS::FASTAFile::FASTAEntry',
            # Constructor parameter types
            'CVMappingFile': 'OpenMS::CVMappingFile',
            'FeatureFinderIdentificationAlgorithm': 'OpenMS::FeatureFinderIdentificationAlgorithm',
            'FeatureFinderAlgorithmMetaboIdent': 'OpenMS::FeatureFinderAlgorithmMetaboIdent',
            'ConsensusMapNormalizerAlgorithmMedian': 'OpenMS::ConsensusMapNormalizerAlgorithmMedian',
            'ConsensusMapNormalizerAlgorithmQuantile': 'OpenMS::ConsensusMapNormalizerAlgorithmQuantile',
            'IonSource': 'OpenMS::IonSource',
            'MassAnalyzer': 'OpenMS::MassAnalyzer',
            'IonDetector': 'OpenMS::IonDetector',
            'ControlledVocabulary': 'OpenMS::ControlledVocabulary',
            'CVMappingRule': 'OpenMS::CVMappingRule',
            'CVMappings': 'OpenMS::CVMappings',
            'ConsensusFeature': 'OpenMS::ConsensusFeature',
            'Peak1D': 'OpenMS::Peak1D',
            'Peak2D': 'OpenMS::Peak2D',
            'ChromatogramPeak': 'OpenMS::ChromatogramPeak',
            'MobilityPeak1D': 'OpenMS::MobilityPeak1D',
            'MSSpectrum': 'OpenMS::MSSpectrum',
            'MSChromatogram': 'OpenMS::MSChromatogram',
            'MSExperiment': 'OpenMS::MSExperiment',
            'PeakMap': 'OpenMS::PeakMap',
            # More constructor types
            'SpectrumAccessOpenMS': 'OpenMS::SpectrumAccessOpenMS',
            'SpectrumAccessOpenMSCached': 'OpenMS::SpectrumAccessOpenMSCached',
            'MapAlignmentEvaluationAlgorithmRecall': 'OpenMS::MapAlignmentEvaluationAlgorithmRecall',
            'UInt64': 'OpenMS::UInt64',
            'LightTransition': 'OpenSwath::LightTransition',
            'LightCompound': 'OpenSwath::LightCompound',
            'LightProtein': 'OpenSwath::LightProtein',
            'LightPeptide': 'OpenSwath::LightPeptide',
            'ChromExtractParams': 'OpenMS::ChromExtractParams',
        }

        # OpenMS classes that need namespace (only add if not already namespaced)
        for typedef, qualified in openms_typedefs.items():
            # Match the type when not already qualified
            pattern = r'\b(?<!OpenMS::)(?<!std::)' + typedef + r'\b'
            result = re.sub(pattern, qualified, result)

        # Convert DPosition2 to DPosition<2> (common alias in .pxd)
        result = re.sub(r'\bDPosition2\b', r'OpenMS::DPosition<2>', result)
        result = re.sub(r'\bDPosition1\b', r'OpenMS::DPosition<1>', result)

        # Handle common nested type aliases (convert to actual types)
        # These are typically typedefs in peak classes like OpenMS::Peak1D::PositionType
        # We need to match the full qualified form first, then unqualified
        nested_types = {
            'PositionType': 'double',           # DPosition<N> reduces to coordinate type
            'CoordinateType': 'double',
            'IntensityType': 'float',           # Peak intensity is float in OpenMS
            'MassType': 'double',
        }
        for nested, actual in nested_types.items():
            # First, replace fully-qualified forms like OpenMS::Peak1D::PositionType -> double
            result = re.sub(r'OpenMS::\w+::' + nested + r'\b', actual, result)
            # Then replace any remaining unqualified forms
            result = re.sub(r'\b' + nested + r'\b', actual, result)

        # Handle enum types that need full qualification
        enum_types = {
            'SpectrumType': 'OpenMS::SpectrumSettings::SpectrumType',
            'SpectrumLevel': 'OpenMS::SpectrumLevel',
            'IMFormat': 'OpenMS::IMFormat',
            'DriftTimeUnit': 'OpenMS::DriftTimeUnit',
            'PeakType': 'OpenMS::Peak1D::PeakType',
            'DataType': 'OpenMS::DataValue::DataType',
            'Polarity': 'OpenMS::IonSource::Polarity',
            'InletType': 'OpenMS::IonSource::InletType',
            'IonizationMethod': 'OpenMS::IonSource::IonizationMethod',
            'AnalyzerType': 'OpenMS::MassAnalyzer::AnalyzerType',
            'ResolutionMethod': 'OpenMS::MassAnalyzer::ResolutionMethod',
            'ResolutionType': 'OpenMS::MassAnalyzer::ResolutionType',
            'ScanDirection': 'OpenMS::MassAnalyzer::ScanDirection',
            'ScanLaw': 'OpenMS::MassAnalyzer::ScanLaw',
            'ReflectronState': 'OpenMS::MassAnalyzer::ReflectronState',
            'Type': 'OpenMS::IonDetector::Type',
            'AcquisitionMode': 'OpenMS::IonDetector::AcquisitionMode',
            'DecoyType': 'OpenMS::TargetedExperimentHelper::RetentionTime::RTType',
            'TermSpecificity': 'OpenMS::ResidueModification::TermSpecificity',
        }
        for enum_type, qualified in enum_types.items():
            pattern = r'\b(?<!OpenMS::)(?<!::)' + enum_type + r'\b'
            result = re.sub(pattern, qualified, result)

        # Strip reference qualifiers unless we're preserving them for static_cast
        if not preserve_reference:
            result = result.replace("const ", "").replace(" const", "")
            result = result.replace("&", "").strip()

        return result

    def _escape_string(self, s: str) -> str:
        """Escape a string for use in C++ code."""
        return s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")

    def _qualify_openms_types(self, type_str: str) -> str:
        """Add OpenMS:: prefix to unqualified OpenMS types.

        Handles both the main type and template parameters, e.g.:
        RangeManagerContainer< RangeMZ, RangeIntensity >
        -> OpenMS::RangeManagerContainer<OpenMS::RangeMZ, OpenMS::RangeIntensity>
        """
        import re

        # OpenMS types that need qualification when unqualified
        openms_types = {
            # Range types
            "RangeManagerContainer", "RangeManager",
            "RangeMZ", "RangeRT", "RangeIntensity", "RangeMobility",
            "RangeBase",
            # Settings types
            "SpectrumSettings", "ChromatogramSettings",
            "ExperimentalSettings", "AcquisitionInfo",
            # Data types
            "MetaInfoInterface", "MetaInfo",
            "Peak1D", "Peak2D", "ChromatogramPeak", "MobilityPeak1D",
            "DefaultParamHandler",
            # Other common types
            "String", "DataValue", "ParamValue",
        }

        result = type_str
        for t in openms_types:
            # Match the type when not already qualified with OpenMS::
            # Use word boundaries and negative lookbehind for ::
            pattern = r'(?<![:\w])' + t + r'(?![:\w])'
            result = re.sub(pattern, f'OpenMS::{t}', result)

        return result

    def _write_module_file(
        self, output_path: Path, module_num: int, content: ModuleContent
    ) -> None:
        """Write a module C++ file as a standalone NB_MODULE."""
        lines = []

        # Header comment
        lines.append("// Auto-generated by pyOpenMS2 generator (v2 - libclang)")
        lines.append(f"// Module {module_num}")
        lines.append("")

        # Includes
        for inc in sorted(content.includes):
            lines.append(f"#include {inc}")
        lines.append("")

        # Namespace
        lines.append("namespace nb = nanobind;")
        lines.append("using namespace nb::literals;")
        lines.append("")

        # Standalone NB_MODULE for this sub-module
        lines.append(f'NB_MODULE(_pyopenms2_{module_num}, m) {{')
        lines.append(f'    m.doc() = "pyOpenMS2 module {module_num}";')
        lines.append("")

        # Class bindings
        for binding in content.class_bindings:
            lines.append(binding)
            lines.append("")

        # Enum bindings
        for binding in content.enum_bindings:
            lines.append(binding)
            lines.append("")

        lines.append("}")

        output_path.write_text("\n".join(lines))

    def _write_main_module(self, output_path: Path, num_modules: int) -> None:
        """Write the main module file.

        Note: With standalone NB_MODULE per sub-module, this main module is
        optional. It's kept for backward compatibility but just creates an
        empty module - the actual classes are in _pyopenms2_N modules.
        """
        lines = []

        lines.append("// Auto-generated by pyOpenMS2 generator (v2 - libclang)")
        lines.append("// Main module - placeholder, actual bindings in _pyopenms2_N modules")
        lines.append("")
        lines.append("#include <nanobind/nanobind.h>")
        lines.append("")
        lines.append("namespace nb = nanobind;")
        lines.append("")

        # Main module definition - empty, just for compatibility
        lines.append('NB_MODULE(_pyopenms2, m) {')
        lines.append('    m.doc() = "pyOpenMS2: Python bindings for OpenMS (nanobind)\\n"')
        lines.append('             "Note: Import pyopenms instead for all classes.";')
        lines.append("}")

        output_path.write_text("\n".join(lines))
