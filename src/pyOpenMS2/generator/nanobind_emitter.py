"""
Nanobind C++ Code Emitter for pyOpenMS2

This module generates nanobind C++ binding code from parsed .pxd declarations.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .pxd_parser import ClassDecl, EnumDecl, Method, Parameter
from .type_registry import TypeRegistry, TypeInfo, TypeCategory


@dataclass
class ModuleContent:
    """Content for a single C++ module file."""

    includes: Set[str]
    forward_declarations: List[str]
    class_bindings: List[str]
    enum_bindings: List[str]


class NanobindEmitter:
    """
    Generates nanobind C++ binding code from parsed class declarations.
    """

    def __init__(self, type_registry: TypeRegistry, num_modules: int = 8):
        self.type_registry = type_registry
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
            "==": "__eq__",
            "!=": "__ne__",
            "<": "__lt__",
            "<=": "__le__",
            ">": "__gt__",
            ">=": "__ge__",
            "+": "__add__",
            "-": "__sub__",
            "*": "__mul__",
            "/": "__truediv__",
            "+=": "__iadd__",
            "-=": "__isub__",
            "*=": "__imul__",
            "/=": "__itruediv__",
            "[]": "__getitem__",
            "()": "__call__",
        }

        # Cython method names to actual C++ method names
        # (Cython uses iadd/isub/etc. because it can't directly name operator+=)
        self._cython_to_cpp_method = {
            "iadd": "operator+=",
            "isub": "operator-=",
            "imul": "operator*=",
            "idiv": "operator/=",
        }

    def emit(
        self,
        classes: Dict[str, ClassDecl],
        addons: Dict[str, dict],
        output_dir: Path,
    ) -> None:
        """
        Generate nanobind C++ files.

        Parameters
        ----------
        classes : dict
            Parsed class declarations.
        addons : dict
            Addon information (for identifying performance-critical methods).
        output_dir : Path
            Output directory for generated files.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve types in all classes
        for class_decl in classes.values():
            self._resolve_class_types(class_decl)

        # Distribute classes across modules
        modules = self._distribute_classes(classes)

        # Generate module files
        for i, module_classes in enumerate(modules, 1):
            content = self._generate_module_content(module_classes, classes, addons)
            output_file = output_dir / f"module_{i}.cpp"
            self._write_module_file(output_file, i, content)

        # Generate main module file
        self._write_main_module(output_dir / "main_module.cpp", len(modules))

        # Generate type caster forward declarations
        self._write_type_caster_forwards(output_dir / "type_casters_fwd.h", classes)

    def _resolve_class_types(self, class_decl: ClassDecl) -> None:
        """Resolve all type references in a class declaration."""
        for method in class_decl.methods + class_decl.constructors:
            if method.return_type and method.return_type != "void":
                method.return_type_info = self.type_registry.resolve_type(
                    method.return_type
                )
            for param in method.parameters:
                param.type_info = self.type_registry.resolve_type(param.type_str)

    def _distribute_classes(
        self, classes: Dict[str, ClassDecl]
    ) -> List[List[ClassDecl]]:
        """Distribute classes across modules for parallel compilation."""
        # Sort by name for deterministic distribution
        sorted_classes = sorted(classes.values(), key=lambda c: c.name)

        # Use hash-based distribution for even split
        modules: List[List[ClassDecl]] = [[] for _ in range(self.num_modules)]

        for class_decl in sorted_classes:
            # Hash class name to determine module
            hash_val = int(hashlib.md5(class_decl.name.encode()).hexdigest(), 16)
            module_idx = hash_val % self.num_modules
            modules[module_idx].append(class_decl)

        return modules

    def _generate_module_content(
        self,
        module_classes: List[ClassDecl],
        all_classes: Dict[str, ClassDecl],
        addons: Dict[str, dict],
    ) -> ModuleContent:
        """Generate content for a module."""
        content = ModuleContent(
            includes=self._standard_includes.copy(),
            forward_declarations=[],
            class_bindings=[],
            enum_bindings=[],
        )

        for class_decl in module_classes:
            # Add OpenMS header
            if class_decl.header_file:
                header = class_decl.header_file
                # Don't double-wrap with angle brackets
                if not header.startswith("<"):
                    header = f"<{header}>"
                content.includes.add(header)

            # Generate class binding
            binding = self._generate_class_binding(class_decl, addons.get(class_decl.name, {}))
            if binding:
                content.class_bindings.append(binding)

            # Generate enum bindings
            for enum_decl in class_decl.enums:
                enum_binding = self._generate_enum_binding(enum_decl, class_decl.name)
                content.enum_bindings.append(enum_binding)

        return content

    def _generate_class_binding(
        self, class_decl: ClassDecl, addon_info: dict
    ) -> Optional[str]:
        """Generate nanobind binding code for a class."""
        lines = []

        # Class name and qualified name
        class_name = class_decl.name
        qualified_name = f"{class_decl.namespace}::{class_name}"

        # Skip nested classes for now (they have complex namespace issues)
        # A nested class would have namespace like "OpenMS::OuterClass" or the qualified
        # name would have more than 2 parts
        if "::" in qualified_name.replace("OpenMS::", "", 1):
            return None

        # For now, skip classes with base classes entirely (complex inheritance issues)
        # This is a conservative approach to get a working build
        base_spec = ""
        if class_decl.base_classes:
            return None

        # Start class definition
        # Note: nanobind handles shared_ptr automatically, no need to specify holder
        lines.append(f'    nb::class_<{qualified_name}{base_spec}>(m, "{class_name}")')

        # Add docstring
        if class_decl.doc:
            doc = self._escape_string(class_decl.doc)
            lines[-1] = lines[-1][:-1]  # Remove trailing paren
            lines[-1] += f', "{doc}")'

        # Generate constructors
        for ctor in class_decl.constructors:
            if ctor.wrap_ignore:
                continue
            ctor_code = self._generate_constructor(ctor, qualified_name)
            lines.append(f"        {ctor_code}")

        # Generate methods - filter out default parameter variants
        # (Cython treats methods with defaults as separate overloads, but C++ doesn't)
        filtered_methods = self._filter_default_param_overloads(class_decl.methods)
        for method in filtered_methods:
            if method.wrap_ignore:
                continue
            method_code = self._generate_method(method, qualified_name, class_decl)
            if method_code:
                lines.append(f"        {method_code}")

        # Add hash support
        if class_decl.wrap_hash:
            lines.append(f'        .def("__hash__", [](const {qualified_name}& self) {{ return std::hash<{qualified_name}>{{}}(self); }})')

        # Add iterator support
        if class_decl.wrap_iter:
            iter_type = class_decl.wrap_iter
            lines.append(f'        .def("__iter__", [](const {qualified_name}& self) {{ return nb::make_iterator(self.begin(), self.end()); }}, nb::keep_alive<0, 1>())')

        # Close class definition
        lines.append("        ;")

        return "\n".join(lines)

    def _filter_default_param_overloads(self, methods: List[Method]) -> List[Method]:
        """
        Filter out methods that appear to be default parameter variants.

        In Cython .pxd files, methods with default parameters are declared as
        separate overloads. But in C++, they're a single method with defaults.
        This function keeps only the most complete version of each method.
        """
        from collections import defaultdict

        # Group methods by name
        by_name: dict = defaultdict(list)
        for method in methods:
            by_name[method.name].append(method)

        filtered = []
        for name, overloads in by_name.items():
            if len(overloads) == 1:
                filtered.append(overloads[0])
                continue

            # Sort by parameter count descending
            sorted_overloads = sorted(overloads, key=lambda m: len(m.parameters), reverse=True)

            # Check if these look like default parameter variants
            # (each subsequent has fewer params, and param types match prefix)
            is_default_params = True
            max_params = sorted_overloads[0].parameters
            for i, m in enumerate(sorted_overloads[1:], 1):
                # Check if this method's params are a prefix of max_params
                if len(m.parameters) >= len(max_params):
                    is_default_params = False
                    break
                for j, p in enumerate(m.parameters):
                    if j >= len(max_params):
                        is_default_params = False
                        break
                    # Compare type strings (rough check)
                    if p.type_str != max_params[j].type_str:
                        is_default_params = False
                        break
                if not is_default_params:
                    break

            if is_default_params and len(sorted_overloads) > 1:
                # Keep only the most complete version
                filtered.append(sorted_overloads[0])
            else:
                # These are true overloads - keep all
                filtered.extend(overloads)

        return filtered

    def _generate_constructor(self, ctor: Method, qualified_name: str) -> str:
        """Generate constructor binding."""
        if not ctor.parameters:
            return ".def(nb::init<>())"

        param_types = []
        for param in ctor.parameters:
            if param.type_info:
                param_types.append(self.type_registry.get_cpp_type(param.type_info))
            else:
                param_types.append(param.type_str)

        types_str = ", ".join(param_types)
        return f".def(nb::init<{types_str}>())"

    def _generate_method(
        self, method: Method, qualified_name: str, class_decl: ClassDecl
    ) -> Optional[str]:
        """Generate method binding."""
        method_name = method.wrap_as or method.name

        # Handle operators
        if method.is_operator and method.operator_name:
            py_op = self._operator_map.get(method.operator_name)
            if py_op:
                return self._generate_operator_binding(method, qualified_name, py_op)
            return None

        # Static methods
        if method.is_static:
            return self._generate_static_method(method, qualified_name, method_name)

        # Regular methods
        return self._generate_regular_method(method, qualified_name, method_name, class_decl)

    def _generate_regular_method(
        self, method: Method, qualified_name: str, method_name: str, class_decl: ClassDecl
    ) -> str:
        """Generate binding for a regular method."""
        # Build parameter types for method pointer
        param_types = []
        param_names = []
        for param in method.parameters:
            if param.type_info:
                param_types.append(self.type_registry.get_cpp_type(param.type_info))
            else:
                param_types.append(param.type_str)
            param_names.append(f'"{param.name}"_a')

        # Return type
        return_type = "void"
        if method.return_type_info:
            return_type = self.type_registry.get_cpp_type(method.return_type_info)
        elif method.return_type:
            return_type = method.return_type

        # Determine actual C++ method name (may differ from Cython method name)
        cpp_method_name = self._cython_to_cpp_method.get(method.name, method.name)

        # For in-place operators (iadd, isub, etc.), use nanobind's operator syntax
        if method.name in self._cython_to_cpp_method:
            op_map = {
                "iadd": "+=",
                "isub": "-=",
                "imul": "*=",
                "idiv": "/=",
            }
            op_symbol = op_map.get(method.name)
            if op_symbol:
                return f".def(nb::self {op_symbol} nb::self)"

        # Check if method is overloaded in this class
        method_count = sum(1 for m in class_decl.methods if m.name == method.name)
        is_overloaded = method_count > 1

        # Build method pointer
        params_str = ", ".join(param_types)

        if is_overloaded:
            # TODO: overload_cast needs exact C++ signatures including const&
            # For now, skip overloaded methods as we can't reliably determine signatures
            # from .pxd files alone
            return None
        else:
            # Direct method pointer for non-overloaded methods
            # Skip if return type suggests const/non-const overloads exist
            # (references, vectors, maps often have both versions)
            if return_type:
                if "&" in return_type:
                    return None
                # Cython uses libcpp_vector, libcpp_map, etc.
                if "libcpp_" in return_type.lower():
                    return None
                if "vector" in return_type.lower():
                    return None

            # For now, only bind simple methods with no parameters
            # This avoids many edge cases with overloads, const/non-const, static methods
            if method.parameters:
                return None

            method_ptr = f"&{qualified_name}::{cpp_method_name}"

        # Build def call
        args_str = ", ".join(param_names) if param_names else ""

        result = f'.def("{method_name}", {method_ptr}'
        if args_str:
            result += f", {args_str}"
        if method.doc:
            doc = self._escape_string(method.doc)
            result += f', "{doc}"'
        result += ")"

        return result

    def _is_likely_const_method(self, method: Method) -> bool:
        """Heuristic to determine if a method is likely const."""
        name = method.name.lower()
        # Common const method patterns
        const_patterns = [
            'get', 'is', 'has', 'can', 'should', 'empty', 'size',
            'length', 'count', 'begin', 'end', 'front', 'back',
            'find', 'contains', 'check', 'valid', 'to', 'at',
            'operator==', 'operator!=', 'operator<', 'operator>',
            'operator<=', 'operator>=', 'operator[]',
        ]
        for pattern in const_patterns:
            if name.startswith(pattern) or name == pattern:
                return True

        # Methods returning references to internal state might be non-const
        # Methods starting with 'set', 'add', 'remove', 'clear', 'push', 'pop' are non-const
        nonconst_patterns = ['set', 'add', 'remove', 'clear', 'push', 'pop', 'insert', 'erase', 'update', 'modify']
        for pattern in nonconst_patterns:
            if name.startswith(pattern):
                return False

        # Default: assume const if no parameters and not void return
        if not method.parameters and method.return_type and method.return_type != "void":
            return True

        return False

    def _generate_static_method(
        self, method: Method, qualified_name: str, method_name: str
    ) -> str:
        """Generate binding for a static method."""
        # For now, skip static methods - they have the same overload_cast issues
        # TODO: implement proper static method binding with correct signatures
        return None

    def _generate_operator_binding(
        self, method: Method, qualified_name: str, py_op: str
    ) -> str:
        """Generate binding for an operator."""
        op = method.operator_name

        if op == "[]":
            # Subscript operator - need lambda for bounds checking
            return f'.def("__getitem__", []({qualified_name}& self, size_t i) {{ return self[i]; }})'
        elif op in ["==", "!=", "<", "<=", ">", ">="]:
            # Comparison operators
            return f'.def(nb::self {op} nb::self)'
        elif op in ["+", "-", "*", "/"]:
            # Arithmetic operators
            return f'.def(nb::self {op} nb::self)'
        elif op in ["+=", "-=", "*=", "/="]:
            # In-place operators
            return f'.def(nb::self {op} nb::self)'

        return ""

    def _generate_enum_binding(self, enum_decl: EnumDecl, parent_class: str) -> str:
        """Generate binding for an enum."""
        lines = []
        qualified_name = f"{enum_decl.namespace}::{enum_decl.name}"
        parent_qualified = f"{enum_decl.namespace}::{parent_class}"

        lines.append(f'    nb::enum_<{qualified_name}>(cls_{parent_class.lower()}, "{enum_decl.name}")')

        for value in enum_decl.values:
            lines.append(f'        .value("{value.name}", {qualified_name}::{value.name})')

        lines.append("        .export_values();")

        return "\n".join(lines)

    def _write_module_file(
        self, output_path: Path, module_num: int, content: ModuleContent
    ) -> None:
        """Write a module C++ file."""
        lines = []

        # Header comment
        lines.append("// Auto-generated by pyOpenMS2 generator")
        lines.append(f"// Module {module_num}")
        lines.append("")

        # Includes
        for inc in sorted(content.includes):
            if inc.startswith("<"):
                lines.append(f"#include {inc}")
            else:
                lines.append(f'#include {inc}')
        lines.append("")

        # Namespace
        lines.append("namespace nb = nanobind;")
        lines.append("using namespace nb::literals;")
        lines.append("")

        # Module function
        lines.append(f"void bind_module_{module_num}(nb::module_& m) {{")

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
        """Write the main module file that imports all submodules."""
        lines = []

        lines.append("// Auto-generated by pyOpenMS2 generator")
        lines.append("// Main module entry point")
        lines.append("")
        lines.append("#include <nanobind/nanobind.h>")
        lines.append("")
        lines.append("namespace nb = nanobind;")
        lines.append("")

        # Forward declarations
        for i in range(1, num_modules + 1):
            lines.append(f"void bind_module_{i}(nb::module_& m);")
        lines.append("")

        # Main module definition
        lines.append('NB_MODULE(_pyopenms2, m) {')
        lines.append('    m.doc() = "pyOpenMS2: Python bindings for OpenMS using nanobind";')
        lines.append("")

        for i in range(1, num_modules + 1):
            lines.append(f"    bind_module_{i}(m);")

        lines.append("}")

        output_path.write_text("\n".join(lines))

    def _write_type_caster_forwards(
        self, output_path: Path, classes: Dict[str, ClassDecl]
    ) -> None:
        """Write forward declarations for type casters."""
        lines = []

        lines.append("// Auto-generated by pyOpenMS2 generator")
        lines.append("// Forward declarations for nanobind type visibility")
        lines.append("")
        lines.append("#pragma once")
        lines.append("")
        lines.append("#include <nanobind/nanobind.h>")
        lines.append("")

        # Forward declare all wrapped classes in OpenMS namespace
        lines.append("namespace OpenMS {")
        for class_name in sorted(classes.keys()):
            lines.append(f"    class {class_name};")
        lines.append("}")
        lines.append("")

        output_path.write_text("\n".join(lines))

    def _escape_string(self, s: str) -> str:
        """Escape a string for C++ string literal."""
        return (
            s.replace("\\", "\\\\")
            .replace('"', '\\"')
            .replace("\n", "\\n")
            .replace("\t", "\\t")
        )
