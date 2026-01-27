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

    def __init__(self, num_modules: int = 8):
        self.num_modules = num_modules

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

            # Generate class binding
            binding = self._generate_class_binding(merged_class)
            if binding:
                content.class_bindings.append(binding)

        return content

    def _generate_class_binding(self, merged_class: MergedClass) -> Optional[str]:
        """Generate nanobind binding code for a class."""
        lines = []

        class_name = merged_class.name
        qualified_name = merged_class.qualified_name

        # Skip nested classes
        if "::" in qualified_name.replace("OpenMS::", "", 1):
            return None

        # Skip abstract classes (can't instantiate)
        if merged_class.is_abstract:
            return None

        # Handle base classes - nanobind supports single inheritance
        base_spec = ""
        if merged_class.base_classes:
            # Filter to known OpenMS bases
            openms_bases = [b for b in merged_class.base_classes if "OpenMS::" in b]
            if openms_bases:
                # Use first OpenMS base class
                base_spec = f", {openms_bases[0]}"

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

        # Add iterator support
        if merged_class.wrap_iter:
            lines.append(f'        .def("__iter__", [](const {qualified_name}& self) {{ return nb::make_iterator(self.begin(), self.end()); }}, nb::keep_alive<0, 1>())')

        # Close class definition
        lines.append("        ;")

        return "\n".join(lines)

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
        """Generate binding for a regular method."""
        method = merged_method.cpp_method

        # Build parameter types
        param_types = [self._normalize_type(p.type_str) for p in method.parameters]
        params_str = ", ".join(param_types)

        # Get return type
        return_type = self._normalize_type(method.return_type) if method.return_type else "void"

        # Check if method is overloaded
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

        # Build method pointer with overload_cast if needed
        if is_overloaded:
            # Use overload_cast for overloaded methods
            const_suffix = ", nb::const_" if method.is_const else ""
            method_ptr = f"nb::overload_cast<{params_str}>(&{qualified_name}::{method.name}{const_suffix})"
        else:
            method_ptr = f"&{qualified_name}::{method.name}"

        # Build argument names
        param_names = []
        for p in method.parameters:
            name = p.name if p.name and not p.name.startswith("arg") else ""
            if name:
                param_names.append(f'"{name}"_a')

        # Build def call
        args_str = ", ".join(param_names) if param_names else ""

        result = f'.def("{method_name}", {method_ptr}'
        if args_str:
            result += f", {args_str}"
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

    def _normalize_type(self, type_str: str) -> str:
        """Normalize C++ type string for use in bindings."""
        if not type_str:
            return "void"

        # Handle OpenMS namespace types
        # Keep qualified names as-is
        return type_str

    def _escape_string(self, s: str) -> str:
        """Escape a string for use in C++ code."""
        return s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")

    def _write_module_file(
        self, output_path: Path, module_num: int, content: ModuleContent
    ) -> None:
        """Write a module C++ file."""
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

        # Module function
        lines.append(f"void register_module_{module_num}(nb::module_& m) {{")

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
        """Write the main module file."""
        lines = []

        lines.append("// Auto-generated by pyOpenMS2 generator (v2 - libclang)")
        lines.append("// Main module")
        lines.append("")
        lines.append("#include <nanobind/nanobind.h>")
        lines.append("")
        lines.append("namespace nb = nanobind;")
        lines.append("")

        # Declare registration functions
        for i in range(1, num_modules + 1):
            lines.append(f"void register_module_{i}(nb::module_& m);")
        lines.append("")

        # Main module definition
        lines.append('NB_MODULE(_pyopenms2, m) {')
        lines.append('    m.doc() = "pyOpenMS2: Python bindings for OpenMS (nanobind)";')
        lines.append("")

        # Call all registration functions
        for i in range(1, num_modules + 1):
            lines.append(f"    register_module_{i}(m);")

        lines.append("}")

        output_path.write_text("\n".join(lines))
