"""
Addon Processor for pyOpenMS2

This module processes existing .pyx addon files from pyOpenMS and determines:
1. Which addons can be converted to pure Python
2. Which addons need C++ implementation for performance
3. Generate Python addon code for injection at import time
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import List, Optional, Tuple


class AddonTier(Enum):
    """Classification of addons by implementation strategy."""

    PURE_PYTHON = auto()  # Can be implemented in pure Python
    CPP_PERFORMANCE = auto()  # Needs C++ for performance
    CPP_REQUIRED = auto()  # Requires C++ (uses cdef, memoryview, etc.)


@dataclass
class AddonMethod:
    """A method defined in an addon file."""

    name: str
    tier: AddonTier
    source_code: str
    parameters: List[str] = field(default_factory=list)
    uses_cdef: bool = False
    uses_memoryview: bool = False
    uses_numpy_buffer: bool = False
    calls_cpp_directly: bool = False
    doc: str = ""


@dataclass
class AddonFile:
    """Parsed content of an addon file."""

    class_name: str
    methods: List[AddonMethod] = field(default_factory=list)
    imports: List[str] = field(default_factory=list)
    cimports: List[str] = field(default_factory=list)


class AddonProcessor:
    """
    Processes .pyx addon files to determine implementation strategy.

    Addon files in pyOpenMS add custom Python methods to wrapped C++ classes.
    This processor analyzes them to determine which can be pure Python
    (injected at import time) vs which need C++ implementation.
    """

    # Methods that should always be in C++ for performance
    PERFORMANCE_CRITICAL_METHODS = {
        "get_peaks",
        "set_peaks",
        "get_data_dict",
        "get_data_arrays",
        "get_df",
        "to_df",
        "get_mz_array",
        "get_intensity_array",
        "get_rt_array",
    }

    # Cython constructs that require C++ implementation
    CPP_REQUIRED_PATTERNS = [
        r"\bcdef\b",
        r"\bcpdef\b",
        r"\bctypedef\b",
        r"\bfrom\s+\w+\s+cimport\b",
        r"\bcimport\b",
        r"\[:\]",  # memoryview syntax
        r"\.inst\.get\(\)",  # Direct access to shared_ptr
        r"deref\(",  # Cython dereference
        r"address\(",  # Cython address-of
    ]

    def __init__(self):
        self._cpp_required_re = [re.compile(p) for p in self.CPP_REQUIRED_PATTERNS]

    def parse_file(self, addon_path: Path) -> Optional[AddonFile]:
        """
        Parse an addon file and classify its methods.

        Parameters
        ----------
        addon_path : Path
            Path to the .pyx addon file.

        Returns
        -------
        AddonFile or None
            Parsed addon information, or None if parsing fails.
        """
        class_name = addon_path.stem

        try:
            content = addon_path.read_text(encoding="utf-8")
        except Exception:
            return None

        addon = AddonFile(class_name=class_name)

        # Extract imports
        addon.imports = self._extract_imports(content)
        addon.cimports = self._extract_cimports(content)

        # Extract and classify methods
        methods = self._extract_methods(content)
        for method_name, method_source in methods:
            method = self._classify_method(method_name, method_source)
            addon.methods.append(method)

        return addon

    def _extract_imports(self, content: str) -> List[str]:
        """Extract Python import statements."""
        imports = []
        for match in re.finditer(r"^(import\s+.+|from\s+\S+\s+import\s+.+)$", content, re.MULTILINE):
            line = match.group(1)
            if "cimport" not in line:
                imports.append(line)
        return imports

    def _extract_cimports(self, content: str) -> List[str]:
        """Extract Cython cimport statements."""
        cimports = []
        for match in re.finditer(r"^(from\s+\S+\s+cimport\s+.+|cimport\s+.+)$", content, re.MULTILINE):
            cimports.append(match.group(1))
        return cimports

    def _extract_methods(self, content: str) -> List[Tuple[str, str]]:
        """Extract method definitions from addon file."""
        methods = []

        # Pattern for def statements at the start of a line or indented
        # Addons don't have class definitions - methods are injected directly
        method_pattern = re.compile(
            r"^([ \t]*)def\s+(\w+)\s*\(([^)]*)\)\s*:",
            re.MULTILINE
        )

        for match in method_pattern.finditer(content):
            indent = match.group(1)
            method_name = match.group(2)
            params = match.group(3)

            # Find the method body
            start = match.end()
            body_indent = len(indent) + 4  # Assume 4-space indentation

            # Find where the method ends
            lines = content[start:].split("\n")
            method_lines = []

            for line in lines:
                if not line.strip():
                    method_lines.append(line)
                    continue

                line_indent = len(line) - len(line.lstrip())
                if line_indent < body_indent and line.strip():
                    break

                method_lines.append(line)

            method_source = f"def {method_name}({params}):\n" + "\n".join(method_lines)
            methods.append((method_name, method_source))

        return methods

    def _classify_method(self, name: str, source: str) -> AddonMethod:
        """Classify a method by its implementation requirements."""
        method = AddonMethod(name=name, source_code=source, tier=AddonTier.PURE_PYTHON)

        # Check for performance-critical methods
        if name in self.PERFORMANCE_CRITICAL_METHODS:
            method.tier = AddonTier.CPP_PERFORMANCE
            return method

        # Check for Cython-specific constructs
        for pattern in self._cpp_required_re:
            if pattern.search(source):
                method.tier = AddonTier.CPP_REQUIRED
                method.uses_cdef = "cdef" in source
                method.uses_memoryview = "[:]" in source
                break

        # Check for direct C++ access patterns
        if ".inst.get()" in source or "deref(" in source:
            method.calls_cpp_directly = True
            if method.tier == AddonTier.PURE_PYTHON:
                method.tier = AddonTier.CPP_REQUIRED

        # Check for numpy buffer protocol usage
        if "np.asarray" in source or "numpy.asarray" in source:
            method.uses_numpy_buffer = True

        # Extract docstring
        doc_match = re.search(r'"""(.+?)"""', source, re.DOTALL)
        if doc_match:
            method.doc = doc_match.group(1).strip()
        else:
            doc_match = re.search(r"'''(.+?)'''", source, re.DOTALL)
            if doc_match:
                method.doc = doc_match.group(1).strip()

        return method

    def generate_python_addon(self, addon: AddonFile) -> str:
        """
        Generate Python addon code for pure Python methods.

        Parameters
        ----------
        addon : AddonFile
            The parsed addon file.

        Returns
        -------
        str
            Python code that can be injected into the module.
        """
        lines = []

        lines.append(f'"""')
        lines.append(f"Addon methods for {addon.class_name}")
        lines.append(f'"""')
        lines.append("")
        lines.append("from __future__ import annotations")
        lines.append("")

        # Add imports (excluding cimports)
        if addon.imports:
            for imp in addon.imports:
                lines.append(imp)
            lines.append("")

        lines.append("from . import addon")
        lines.append("")

        # Generate addon methods
        for method in addon.methods:
            if method.tier == AddonTier.PURE_PYTHON:
                # Convert Cython method to pure Python
                converted = self._convert_to_pure_python(method)
                lines.append(f'@addon("{addon.class_name}")')
                lines.append(converted)
                lines.append("")

        return "\n".join(lines)

    def _convert_to_pure_python(self, method: AddonMethod) -> str:
        """
        Convert a Cython method to pure Python.

        This handles simple cases. Complex methods that use cdef, memoryviews,
        etc. will have been classified as CPP_REQUIRED.
        """
        source = method.source_code

        # Remove type annotations from parameters
        source = re.sub(r"(\w+)\s*:\s*\w+(\s*[,)])", r"\1\2", source)

        # Remove cdef local variable declarations
        source = re.sub(r"^\s*cdef\s+\w+\s+\w+\s*$", "", source, flags=re.MULTILINE)

        # Replace self.inst.get() patterns with self
        # These shouldn't appear in pure Python methods, but just in case
        source = source.replace("self.inst.get()", "self._cpp")

        return source

    def get_cpp_required_methods(self, addon: AddonFile) -> List[AddonMethod]:
        """Get methods that require C++ implementation."""
        return [m for m in addon.methods if m.tier in (AddonTier.CPP_REQUIRED, AddonTier.CPP_PERFORMANCE)]

    def get_pure_python_methods(self, addon: AddonFile) -> List[AddonMethod]:
        """Get methods that can be pure Python."""
        return [m for m in addon.methods if m.tier == AddonTier.PURE_PYTHON]

    def generate_cpp_addon_bindings(self, addon: AddonFile) -> str:
        """
        Generate C++ nanobind code for performance-critical addon methods.

        This is a placeholder - actual implementation would need to translate
        the Cython code to C++ nanobind lambda expressions.
        """
        lines = []

        cpp_methods = self.get_cpp_required_methods(addon)
        if not cpp_methods:
            return ""

        lines.append(f"// C++ addon methods for {addon.class_name}")
        lines.append(f"// These methods require C++ implementation for performance")
        lines.append("")

        for method in cpp_methods:
            lines.append(f"// {method.name}")
            if method.name == "get_peaks":
                lines.append(self._generate_get_peaks(addon.class_name))
            elif method.name == "set_peaks":
                lines.append(self._generate_set_peaks(addon.class_name))
            elif method.name == "__repr__":
                lines.append(self._generate_repr(addon.class_name))
            elif method.name == "__len__":
                lines.append(self._generate_len(addon.class_name))
            elif method.name == "__iter__":
                lines.append(self._generate_iter(addon.class_name))
            else:
                lines.append(f"// TODO: Implement {method.name}")
            lines.append("")

        return "\n".join(lines)

    def _generate_get_peaks(self, class_name: str) -> str:
        """Generate C++ code for get_peaks method."""
        return f'''.def("get_peaks", [](const OpenMS::{class_name}& self) {{
        size_t n = self.size();
        auto mz = nb::ndarray<double, nb::shape<-1>>({{n}});
        auto intensity = nb::ndarray<double, nb::shape<-1>>({{n}});

        double* mz_data = mz.data();
        double* int_data = intensity.data();

        for (size_t i = 0; i < n; ++i) {{
            mz_data[i] = self[i].getMZ();
            int_data[i] = self[i].getIntensity();
        }}

        return nb::make_tuple(mz, intensity);
    }}, "Get peaks as (mz_array, intensity_array) tuple")'''

    def _generate_set_peaks(self, class_name: str) -> str:
        """Generate C++ code for set_peaks method."""
        return f'''.def("set_peaks", [](OpenMS::{class_name}& self,
                         nb::ndarray<double, nb::shape<-1>>& mz,
                         nb::ndarray<double, nb::shape<-1>>& intensity) {{
        size_t n = mz.shape(0);
        if (n != intensity.shape(0)) {{
            throw std::runtime_error("mz and intensity arrays must have same length");
        }}

        self.resize(n);
        const double* mz_data = mz.data();
        const double* int_data = intensity.data();

        for (size_t i = 0; i < n; ++i) {{
            self[i].setMZ(mz_data[i]);
            self[i].setIntensity(int_data[i]);
        }}
    }}, nb::arg("mz"), nb::arg("intensity"), "Set peaks from mz and intensity arrays")'''

    def _generate_repr(self, class_name: str) -> str:
        """Generate C++ code for __repr__ method."""
        return f'''.def("__repr__", [](const OpenMS::{class_name}& self) {{
        std::ostringstream oss;
        oss << "{class_name}(...)";
        return oss.str();
    }})'''

    def _generate_len(self, class_name: str) -> str:
        """Generate C++ code for __len__ method."""
        return f'''.def("__len__", [](const OpenMS::{class_name}& self) {{
        return self.size();
    }})'''

    def _generate_iter(self, class_name: str) -> str:
        """Generate C++ code for __iter__ method."""
        return f'''.def("__iter__", [](const OpenMS::{class_name}& self) {{
        return nb::make_iterator(self.begin(), self.end());
    }}, nb::keep_alive<0, 1>())'''
