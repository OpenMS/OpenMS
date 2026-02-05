"""
PXD Parser for pyOpenMS2 nanobind bindings

This module parses existing Cython .pxd files from pyOpenMS and extracts
class declarations, method signatures, and wrap- directives for generating
nanobind bindings.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .type_registry import TypeRegistry, TypeInfo

logger = logging.getLogger(__name__)


class WrapDirective(Enum):
    """Types of wrap- directives in .pxd files."""

    DOC = auto()  # wrap-doc: Documentation
    AS = auto()  # wrap-as: Rename method
    IGNORE = auto()  # wrap-ignore: Skip wrapping
    INSTANCES = auto()  # wrap-instances: Template specializations
    INHERITS = auto()  # wrap-inherits: Parent classes
    UPPER_LIMIT = auto()  # wrap-upper-limit: Bounds checking
    ATTACH = auto()  # wrap-attach: Enums/static methods
    MANUAL_MEMORY = auto()  # wrap-manual-memory: Custom memory management
    HASH = auto()  # wrap-hash: Hash function
    ITER = auto()  # wrap-iter: Iterator protocol
    CONST = auto()  # wrap-const: Mark method as const (for nanobind)
    STATIC = auto()  # wrap-static: Mark method as static (for nanobind)


@dataclass
class Parameter:
    """A function/method parameter."""

    name: str
    type_str: str
    type_info: Optional[TypeInfo] = None
    default_value: Optional[str] = None
    is_output: bool = False  # For output parameters (passed by reference)


@dataclass
class Method:
    """A class method declaration."""

    name: str
    return_type: str
    return_type_info: Optional[TypeInfo] = None
    parameters: List[Parameter] = field(default_factory=list)
    is_static: bool = False
    is_const: bool = False
    is_constructor: bool = False
    is_destructor: bool = False
    is_operator: bool = False
    operator_name: Optional[str] = None
    doc: str = ""
    wrap_as: Optional[str] = None  # Renamed method
    wrap_ignore: bool = False
    wrap_upper_limit: Optional[str] = None
    nogil: bool = True
    except_clause: bool = True
    cpp_name: Optional[str] = None  # C++ qualified name for standalone functions


@dataclass
class EnumValue:
    """An enum value."""

    name: str
    value: Optional[int] = None
    doc: str = ""


@dataclass
class MemberVariable:
    """A public member variable declaration."""

    name: str
    type_str: str


@dataclass
class EnumDecl:
    """An enum declaration."""

    name: str
    values: List[EnumValue] = field(default_factory=list)
    doc: str = ""
    namespace: str = "OpenMS"
    attached_to: Optional[str] = None  # Class this enum is attached to
    cpp_name: Optional[str] = None  # Full C++ qualified name (e.g., "OpenMS::FileTypes::Type")
    is_scoped: bool = False  # True for 'enum class' (C++11 scoped enums)


@dataclass
class ClassDecl:
    """A class declaration."""

    name: str
    namespace: str = "OpenMS"
    header_file: str = ""
    doc: str = ""
    methods: List[Method] = field(default_factory=list)
    constructors: List[Method] = field(default_factory=list)
    base_classes: List[str] = field(default_factory=list)
    enums: List[EnumDecl] = field(default_factory=list)
    template_args: List[str] = field(default_factory=list)
    template_instances: Dict[str, List[str]] = field(default_factory=dict)
    is_abstract: bool = False
    wrap_hash: Optional[str] = None
    wrap_manual_memory: bool = False
    wrap_iter: Optional[str] = None
    attached_enums: List[str] = field(default_factory=list)
    member_variables: List[MemberVariable] = field(default_factory=list)
    cpp_qualified_name: Optional[str] = None  # Full C++ name from pxd "..." syntax (e.g., "OpenMS::PeakIntegrator::PeakArea")


class PxdParser:
    """
    Parser for Cython .pxd declaration files.

    This parser extracts class declarations, method signatures, and wrap-
    directives from existing pyOpenMS .pxd files to generate nanobind bindings.
    """

    def __init__(self, type_registry: TypeRegistry):
        self.type_registry = type_registry
        self._current_file: Optional[Path] = None
        self._current_class: Optional[str] = None

        # Regex patterns
        self._extern_pattern = re.compile(
            r'cdef\s+extern\s+from\s+"([^"]+)"\s+namespace\s+"([^"]+)":'
        )
        self._class_pattern = re.compile(
            r'cdef\s+cppclass\s+(\w+)(?:\s+"([^"]+)")?(?:\[([^\]]+)\])?\s*(?:\(([^)]+)\))?'
        )
        self._method_pattern = re.compile(
            r"(\w+)\s*\(([^)]*)\)\s*(except\s*\+)?\s*(nogil)?\s*(?:#.*)?$"
        )
        self._static_method_pattern = re.compile(
            r"@staticmethod"
        )
        # Pattern for enum: cdef enum EnumName or cdef enum class EnumName "C++Name"
        # Note: "class" is Cython syntax for scoped enums, not the enum name
        self._enum_pattern = re.compile(
            r'cdef\s+enum\s+(?:class\s+)?(\w+)(?:\s+"([^"]+)")?'
        )
        self._scoped_enum_pattern = re.compile(r'cdef\s+enum\s+class\s+')
        self._wrap_directive_pattern = re.compile(
            r"#?\s*wrap-(\w+)(?::\s*(.*?))?(?=\s+wrap-|\s*$)"
        )
        # Pattern for standalone functions with C++ alias:
        # return_type func_name "OpenMS::Class::func_name" (params) except + nogil
        self._standalone_func_pattern = re.compile(
            r'^\s*(.+?)\s+(\w+)\s+"([^"]+)"\s*\(([^)]*)\)\s*(except\s*\+)?\s*(nogil)?'
        )
        # Pattern for namespace-level functions without C++ alias (using wrap-attach or namespace):
        # return_type func_name(params) except + nogil  # wrap-attach:ClassName
        # Note: return type can include brackets like libcpp_vector[Type]
        self._namespace_func_pattern = re.compile(
            r'^\s*(?!cdef\s|def\s|class\s|ctypedef\s|@|#)([A-Za-z_][\w:<>,\s\*&\[\]]*?)\s+(\w+)\s*\(([^)]*)\)\s*(except\s*\+)?\s*(nogil)?\s*(?:#.*)?$'
        )

    def parse_file(self, pxd_path: Path) -> Dict[str, ClassDecl]:
        """
        Parse a .pxd file and extract class declarations.

        Parameters
        ----------
        pxd_path : Path
            Path to the .pxd file.

        Returns
        -------
        dict
            Dictionary mapping class names to ClassDecl objects.
        """
        self._current_file = pxd_path

        with open(pxd_path, "r", encoding="utf-8") as f:
            content = f.read()

        classes: Dict[str, ClassDecl] = {}
        enums: List[EnumDecl] = []
        # Standalone functions to be attached to classes as static methods
        # List of (target_class_name, method, header_file)
        standalone_funcs: List[Tuple[str, Method, str]] = []

        lines = content.split("\n")
        i = 0

        current_namespace = "OpenMS"
        current_header = ""

        while i < len(lines):
            line = lines[i]
            stripped = line.strip()

            # Skip empty lines and comments (unless wrap directive)
            if not stripped or (stripped.startswith("#") and "wrap-" not in stripped):
                i += 1
                continue

            # Check for extern from declaration
            extern_match = self._extern_pattern.match(stripped)
            if extern_match:
                current_header = extern_match.group(1)
                current_namespace = extern_match.group(2)
                i += 1
                continue

            # Check for class declaration
            class_match = self._class_pattern.match(stripped)
            if class_match:
                class_name = class_match.group(1)
                cpp_qualified_name = class_match.group(2)  # Optional "OpenMS::Foo::Bar"
                template_args_str = class_match.group(3)
                base_classes_str = class_match.group(4)

                # Determine the C++ namespace from qualified name if provided
                if cpp_qualified_name:
                    # Extract namespace from "OpenMS::Foo::Bar" -> "OpenMS::Foo"
                    parts = cpp_qualified_name.rsplit("::", 1)
                    if len(parts) == 2:
                        current_namespace = parts[0]

                template_args = []
                if template_args_str:
                    template_args = [a.strip() for a in template_args_str.split(",")]

                base_classes = []
                if base_classes_str:
                    base_classes = [b.strip() for b in base_classes_str.split(",")]

                # Parse class body
                class_decl, end_idx = self._parse_class_body(
                    lines, i, class_name, current_namespace, current_header,
                    template_args, base_classes, cpp_qualified_name=cpp_qualified_name
                )

                if class_decl:
                    classes[class_name] = class_decl
                    self.type_registry.register_class(class_name, current_namespace)

                i = end_idx
                continue

            # Check for enum declaration
            enum_match = self._enum_pattern.match(stripped)
            if enum_match:
                enum_name = enum_match.group(1)
                cpp_name = enum_match.group(2)  # Optional C++ alias like "OpenMS::FileTypes::Type"
                # Check if this is a scoped enum (cdef enum class)
                is_scoped = bool(self._scoped_enum_pattern.match(stripped))
                enum_decl, end_idx = self._parse_enum_body(
                    lines, i, enum_name, current_namespace
                )
                if enum_decl:
                    if cpp_name:
                        enum_decl.cpp_name = cpp_name
                    enum_decl.is_scoped = is_scoped
                    enums.append(enum_decl)
                    self.type_registry.register_enum(enum_name, current_namespace)
                i = end_idx
                continue

            # Check for standalone function with C++ alias (namespace-level static functions)
            # Pattern: return_type func_name "OpenMS::Class::func_name" (params) except + nogil
            func_match = self._standalone_func_pattern.match(stripped)
            if func_match:
                return_type = func_match.group(1).strip()
                func_name = func_match.group(2)
                cpp_qualified = func_match.group(3)  # e.g., "OpenMS::PEFFFile::isPEFFFile"
                params_str = func_match.group(4)

                # Extract target class from C++ qualified name or namespace
                # "OpenMS::PEFFFile::isPEFFFile" -> "PEFFFile"
                # "OpenMS::PEFFEntry::fromFASTAEntry" -> "PEFFEntry"
                target_class = None
                if cpp_qualified:
                    parts = cpp_qualified.split("::")
                    if len(parts) >= 2:
                        # Second-to-last part is the class name
                        target_class = parts[-2]
                if not target_class and current_namespace.startswith("OpenMS::"):
                    # Fall back to namespace: "OpenMS::PEFFFile" -> "PEFFFile"
                    target_class = current_namespace.split("::")[-1]

                if target_class:
                    params = self._parse_parameters(params_str)
                    method = Method(
                        name=func_name,
                        return_type=return_type,
                        parameters=params,
                        is_static=True,
                        nogil="nogil" in stripped,
                        except_clause="except" in stripped,
                        cpp_name=cpp_qualified,  # Store original C++ name for binding
                    )

                    # Parse wrap directives from this line and following lines
                    line_directives = self._parse_wrap_directives(stripped)
                    if WrapDirective.DOC in line_directives:
                        method.doc = line_directives[WrapDirective.DOC]
                    if WrapDirective.IGNORE in line_directives:
                        method.wrap_ignore = True

                    # Check for doc continuation on following lines
                    j = i + 1
                    doc_lines = []
                    while j < len(lines):
                        next_line = lines[j].strip()
                        if next_line.startswith("#") and "wrap-doc:" not in next_line:
                            # Continue doc from previous wrap-doc
                            doc_content = next_line.lstrip("#").strip()
                            if doc_content:
                                doc_lines.append(doc_content)
                            j += 1
                        else:
                            break
                    if doc_lines:
                        if method.doc:
                            method.doc += "\n" + "\n".join(doc_lines)
                        else:
                            method.doc = "\n".join(doc_lines)

                    standalone_funcs.append((target_class, method, current_header))
                    logger.debug(f"Found standalone function {func_name} -> {target_class}.{func_name}")

                i += 1
                continue

            # Check for namespace-level function without C++ alias (using wrap-attach or namespace)
            # Pattern: return_type func_name(params) except + nogil  # wrap-attach:ClassName
            namespace_func_match = self._namespace_func_pattern.match(stripped)
            if namespace_func_match and current_namespace:
                return_type = namespace_func_match.group(1).strip()
                func_name = namespace_func_match.group(2)
                params_str = namespace_func_match.group(3)

                # Parse wrap directives from this line
                line_directives = self._parse_wrap_directives(stripped)

                # Determine target class from wrap-attach or namespace
                target_class = None
                if WrapDirective.ATTACH in line_directives:
                    target_class = line_directives[WrapDirective.ATTACH].strip()
                if not target_class and current_namespace.startswith("OpenMS::"):
                    # Extract class from namespace: "OpenMS::ProForma" -> "ProForma"
                    target_class = current_namespace.split("::")[-1]

                if target_class:
                    params = self._parse_parameters(params_str)
                    # Construct fully qualified C++ name from namespace
                    cpp_name = f"{current_namespace}::{func_name}"

                    method = Method(
                        name=func_name,
                        return_type=return_type,
                        parameters=params,
                        is_static=True,
                        nogil="nogil" in stripped,
                        except_clause="except" in stripped,
                        cpp_name=cpp_name,
                    )

                    if WrapDirective.DOC in line_directives:
                        method.doc = line_directives[WrapDirective.DOC]
                    if WrapDirective.IGNORE in line_directives:
                        method.wrap_ignore = True

                    # Check for doc continuation on following lines
                    j = i + 1
                    doc_lines = []
                    while j < len(lines):
                        next_line = lines[j].strip()
                        if next_line.startswith("#") and "wrap-doc:" not in next_line:
                            doc_content = next_line.lstrip("#").strip()
                            if doc_content:
                                doc_lines.append(doc_content)
                            j += 1
                        else:
                            break
                    if doc_lines:
                        if method.doc:
                            method.doc += "\n" + "\n".join(doc_lines)
                        else:
                            method.doc = "\n".join(doc_lines)

                    standalone_funcs.append((target_class, method, current_header))
                    logger.debug(
                        f"Found namespace function {func_name} -> {target_class}.{func_name} "
                        f"(from namespace {current_namespace})"
                    )

                i += 1
                continue

            i += 1

        # Attach enums to classes based on wrap-attach directives or file name
        file_stem = pxd_path.stem  # e.g., "IMTypes" from "IMTypes.pxd"
        for enum_decl in enums:
            # First, check explicit attached_to from namespace or wrap-attach
            if enum_decl.attached_to and enum_decl.attached_to in classes:
                classes[enum_decl.attached_to].enums.append(enum_decl)
            # If not attached, try to attach to class matching the file name
            elif not enum_decl.attached_to and file_stem in classes:
                enum_decl.attached_to = file_stem
                classes[file_stem].enums.append(enum_decl)
                logger.debug(f"Auto-attached enum {enum_decl.name} to class {file_stem} (by file name)")

        # Attach standalone functions to classes as static methods
        for target_class, method, header in standalone_funcs:
            if target_class in classes:
                classes[target_class].methods.append(method)
                logger.debug(f"Attached standalone function {method.name} to {target_class}")
            else:
                # Class not found in this file - might be in another file
                # Create a placeholder class if needed, or just log a warning
                logger.warning(
                    f"Target class {target_class} not found for standalone function "
                    f"{method.name} in {pxd_path.name}"
                )

        return classes

    def _parse_class_body(
        self,
        lines: List[str],
        start_idx: int,
        class_name: str,
        namespace: str,
        header_file: str,
        template_args: List[str],
        base_classes: List[str],
        cpp_qualified_name: Optional[str] = None,
    ) -> Tuple[Optional[ClassDecl], int]:
        """Parse the body of a class declaration."""
        self._current_class = class_name

        class_decl = ClassDecl(
            name=class_name,
            namespace=namespace,
            header_file=header_file,
            template_args=template_args,
            base_classes=base_classes,
            cpp_qualified_name=cpp_qualified_name,
        )

        # Get indentation level of class definition
        base_indent = self._get_indent(lines[start_idx])
        i = start_idx + 1

        # Collect wrap directives that appeared before or on class line
        class_line = lines[start_idx]
        class_directives = self._parse_wrap_directives(class_line)

        # Look for directives on previous lines
        for j in range(max(0, start_idx - 10), start_idx):
            prev_line = lines[j].strip()
            if prev_line.startswith("#"):
                directives = self._parse_wrap_directives(prev_line)
                class_directives.update(directives)

        # Apply class-level directives
        if WrapDirective.DOC in class_directives:
            class_decl.doc = class_directives[WrapDirective.DOC]
        if WrapDirective.HASH in class_directives:
            class_decl.wrap_hash = class_directives[WrapDirective.HASH]
        if WrapDirective.MANUAL_MEMORY in class_directives:
            class_decl.wrap_manual_memory = True
        if WrapDirective.INSTANCES in class_directives:
            class_decl.template_instances = self._parse_instances(
                class_directives[WrapDirective.INSTANCES]
            )
        if WrapDirective.INHERITS in class_directives:
            inherits = class_directives[WrapDirective.INHERITS]
            if inherits:
                class_decl.base_classes.extend(
                    [b.strip() for b in inherits.split(",")]
                )
        if WrapDirective.ITER in class_directives:
            class_decl.wrap_iter = class_directives[WrapDirective.ITER]

        # Parse class body
        is_static_next = False
        pending_doc = []
        pending_hash = False  # True when we saw wrap-hash: with no value, expecting continuation
        pending_class_doc = False  # True when we saw wrap-doc: (class-level), expecting continuation
        seen_method = False  # Track whether we've seen any method/constructor yet

        while i < len(lines):
            line = lines[i]
            current_indent = self._get_indent(line)
            stripped = line.strip()

            # Check if we've exited the class body
            if stripped and current_indent <= base_indent and not stripped.startswith("#"):
                break

            # Skip empty lines
            if not stripped:
                i += 1
                continue

            # Check for static method decorator
            if "@staticmethod" in stripped:
                is_static_next = True
                i += 1
                continue

            # Check for wrap directives (including doc continuation)
            if stripped.startswith("#"):
                # Check for ABSTRACT class comment (fallback when libclang can't parse)
                # Only match if the comment refers to THIS class (contains the class name)
                # or is a standalone directive like "# ABSTRACT CLASS"
                if "ABSTRACT" in stripped.upper() and "class" in stripped.lower():
                    upper = stripped.upper()
                    # Match "# ABSTRACT CLASS" pattern or class name mentioned
                    if class_name.upper() in upper or re.match(r'^#\s*ABSTRACT\s+CLASS', upper):
                        class_decl.is_abstract = True
                        logger.debug(f"Class {class_name} marked as abstract from pxd comment")
                directives = self._parse_wrap_directives(stripped)

                # Handle continuation lines for pending directives
                if (pending_hash or pending_class_doc) and WrapDirective.DOC in directives and WrapDirective.HASH not in directives:
                    doc_val = directives[WrapDirective.DOC].strip()
                    if pending_hash:
                        class_decl.wrap_hash = doc_val
                        pending_hash = False
                    elif pending_class_doc:
                        pending_doc.append(doc_val)
                elif WrapDirective.DOC in directives:
                    pending_doc.append(directives[WrapDirective.DOC])

                # Apply class-level directives found inside class body
                if WrapDirective.HASH in directives:
                    # Flush pending class doc before processing hash (which clears pending_class_doc)
                    if pending_class_doc and pending_doc and not seen_method:
                        class_decl.doc = "\n".join(pending_doc)
                        pending_doc = []
                    value = directives[WrapDirective.HASH].strip()
                    if value:
                        class_decl.wrap_hash = value
                    else:
                        pending_hash = True
                    pending_class_doc = False
                if WrapDirective.DOC in directives and not seen_method:
                    # wrap-doc: before any method → class-level doc
                    value = directives[WrapDirective.DOC].strip() if WrapDirective.DOC in directives else ""
                    if not value:
                        pending_class_doc = True
                if WrapDirective.MANUAL_MEMORY in directives:
                    class_decl.wrap_manual_memory = True
                if WrapDirective.ITER in directives:
                    class_decl.wrap_iter = directives[WrapDirective.ITER]
                if WrapDirective.INSTANCES in directives:
                    # wrap-instances: may have value on same line or on continuation lines
                    value = directives[WrapDirective.INSTANCES].strip()
                    if value:
                        class_decl.template_instances.update(self._parse_instances(value))
                    # Collect continuation lines with ":=" format
                    while i + 1 < len(lines):
                        next_line = lines[i + 1].strip()
                        if next_line.startswith("#") and ":=" in next_line:
                            i += 1
                            # Parse "# InstanceName := ClassName[T1, T2]"
                            inst_text = next_line.lstrip("# ")
                            if ":=" in inst_text:
                                inst_name, inst_type = inst_text.split(":=", 1)
                                inst_name = inst_name.strip()
                                inst_type = inst_type.strip()
                                # Extract template args from ClassName[T1, T2]
                                bracket_match = re.match(r"(\w+)\[([^\]]+)\]", inst_type)
                                if bracket_match:
                                    args = [a.strip() for a in bracket_match.group(2).split(",")]
                                    class_decl.template_instances[inst_name] = args
                                    logger.debug(f"Template instance: {inst_name} := {class_name}[{', '.join(args)}]")
                        else:
                            break
                i += 1
                continue

            # Handle multi-line method signatures (parameters spanning multiple lines)
            code_part = stripped.split('#')[0]
            if '(' in code_part and ')' not in code_part:
                # Join continuation lines until we find the closing ')'
                full_line = stripped
                while i + 1 < len(lines):
                    i += 1
                    next_stripped = lines[i].strip()
                    full_line += ' ' + next_stripped
                    if ')' in next_stripped.split('#')[0]:
                        break
                stripped = full_line

            # Parse method/constructor
            method = self._parse_method_line(stripped, class_name)
            if method:
                # On first method, flush any class-level doc from pending_doc
                if not seen_method and pending_doc and pending_class_doc:
                    class_decl.doc = "\n".join(pending_doc)
                    pending_doc = []
                    pending_class_doc = False
                seen_method = True

                method.is_static = is_static_next
                is_static_next = False

                # Apply pending documentation
                if pending_doc:
                    method.doc = "\n".join(pending_doc)
                    pending_doc = []

                # Parse wrap directives from this line
                line_directives = self._parse_wrap_directives(stripped)
                if WrapDirective.IGNORE in line_directives:
                    method.wrap_ignore = True
                if WrapDirective.AS in line_directives:
                    method.wrap_as = line_directives[WrapDirective.AS]
                if WrapDirective.UPPER_LIMIT in line_directives:
                    method.wrap_upper_limit = line_directives[WrapDirective.UPPER_LIMIT]
                if WrapDirective.DOC in line_directives:
                    method.doc = line_directives[WrapDirective.DOC]
                if WrapDirective.CONST in line_directives:
                    method.is_const = True
                if WrapDirective.STATIC in line_directives:
                    method.is_static = True

                # Infer const from method name patterns if not explicitly set
                # and method takes no parameters or only const refs
                if not method.is_const and not method.is_constructor and not method.is_static:
                    method.is_const = self._infer_const(method)

                if method.is_constructor:
                    class_decl.constructors.append(method)
                else:
                    class_decl.methods.append(method)
            else:
                # Try to parse as a member variable: "type_name var_name"
                # Must not contain parentheses (methods), "cdef", "ctypedef", "enum", "pass"
                member_var = self._parse_member_variable(stripped)
                if member_var:
                    class_decl.member_variables.append(member_var)

            i += 1

        self._current_class = None
        return class_decl, i

    def _parse_enum_body(
        self,
        lines: List[str],
        start_idx: int,
        enum_name: str,
        namespace: str,
    ) -> Tuple[Optional[EnumDecl], int]:
        """Parse the body of an enum declaration."""
        enum_decl = EnumDecl(name=enum_name, namespace=namespace)

        base_indent = self._get_indent(lines[start_idx])
        i = start_idx + 1

        # Check for wrap-attach directive
        start_line = lines[start_idx]
        directives = self._parse_wrap_directives(start_line)
        if WrapDirective.ATTACH in directives:
            enum_decl.attached_to = directives[WrapDirective.ATTACH]

        # Also check previous lines for directives
        for j in range(max(0, start_idx - 5), start_idx):
            prev_line = lines[j].strip()
            if prev_line.startswith("#"):
                prev_directives = self._parse_wrap_directives(prev_line)
                if WrapDirective.ATTACH in prev_directives:
                    enum_decl.attached_to = prev_directives[WrapDirective.ATTACH]
                if WrapDirective.DOC in prev_directives:
                    enum_decl.doc = prev_directives[WrapDirective.DOC]

        # Auto-deduce attached_to from namespace if not explicitly set
        # e.g., namespace "OpenMS::FileTypes" -> attached_to "FileTypes"
        if not enum_decl.attached_to and namespace.startswith("OpenMS::"):
            parts = namespace.split("::")
            if len(parts) >= 2:
                # namespace is "OpenMS::ClassName" or "OpenMS::Outer::Inner"
                # Use the first part after "OpenMS" as the class
                potential_class = parts[1]
                # Only set if it looks like a class name (starts uppercase)
                if potential_class and potential_class[0].isupper():
                    enum_decl.attached_to = potential_class

        # Look for wrap-attach in the enum body (multi-line format)
        # Format:
        #   # wrap-attach:
        #   #  ClassName
        pending_attach = False
        while i < len(lines):
            line = lines[i]
            current_indent = self._get_indent(line)
            stripped = line.strip()

            if stripped and current_indent <= base_indent and not stripped.startswith("#"):
                break

            if not stripped or stripped.startswith("#"):
                # Check for wrap-attach in comment lines
                if stripped.startswith("#"):
                    if "wrap-attach:" in stripped:
                        # Check if value is on same line
                        attach_directives = self._parse_wrap_directives(stripped)
                        if WrapDirective.ATTACH in attach_directives:
                            val = attach_directives[WrapDirective.ATTACH].strip()
                            if val:
                                enum_decl.attached_to = val
                            else:
                                pending_attach = True
                    elif pending_attach:
                        # This line should have the class name (e.g., "#  Residue")
                        class_name = stripped.lstrip("#").strip()
                        if class_name and class_name[0].isupper():
                            enum_decl.attached_to = class_name
                            pending_attach = False
                i += 1
                continue

            # Parse enum value - strip trailing comments first
            # e.g., "UNKNOWN,            # < Unknown file extension" -> "UNKNOWN"
            comment_idx = stripped.find("#")
            clean_value = stripped[:comment_idx].strip() if comment_idx >= 0 else stripped

            if not clean_value:
                i += 1
                continue

            # Handle multiple comma-separated values on one line (e.g., "CMD, GUI, NONE")
            values_on_line = [v.strip() for v in clean_value.split(",") if v.strip()]

            for val in values_on_line:
                if not val or val.startswith("cdef"):
                    continue

                # Handle Cython syntax: NAME "fully::qualified::cpp::name"
                # e.g., SUM "OpenMS::MSExperiment::RasterAggregation::SUM"
                # Extract just the Python/enum name before the quoted string
                if '"' in val:
                    val = val.split('"')[0].strip()

                if "=" in val:
                    parts = val.split("=")
                    name = parts[0].strip()
                    value = int(parts[1].strip())
                    enum_decl.values.append(EnumValue(name=name, value=value))
                else:
                    enum_decl.values.append(EnumValue(name=val))

            i += 1

        return enum_decl, i

    def _parse_member_variable(self, line: str) -> Optional[MemberVariable]:
        """Parse a member variable declaration like 'double area' or 'libcpp_string id'."""
        # Remove trailing comments
        comment_idx = line.find("#")
        clean = line[:comment_idx].strip() if comment_idx >= 0 else line.strip()
        if not clean:
            return None
        # Skip keywords and lines with parens (methods), colons (cdef/ctypedef), etc.
        if "(" in clean or ")" in clean:
            return None
        if clean.startswith(("cdef ", "ctypedef ", "enum ", "pass", "cppclass ", "@")):
            return None
        # Match: type_tokens... var_name
        # e.g., "double area", "libcpp_vector[Peak1D] peaks", "Int nn_threshold"
        parts = clean.rsplit(None, 1)
        if len(parts) == 2:
            type_str, var_name = parts
            # var_name must be a simple identifier
            if re.match(r"^[a-zA-Z_]\w*$", var_name) and not re.match(r"^[a-zA-Z_]\w*$", clean):
                return MemberVariable(name=var_name, type_str=type_str)
        return None

    def _parse_method_line(self, line: str, class_name: str) -> Optional[Method]:
        """Parse a method declaration line."""
        # Remove trailing comments (but keep wrap directives for later)
        comment_idx = line.find("#")
        clean_line = line[:comment_idx].strip() if comment_idx >= 0 else line.strip()

        if not clean_line:
            return None

        # Check for constructor
        constructor_match = re.match(
            rf"{class_name}\s*\(([^)]*)\)\s*(except\s*\+)?\s*(nogil)?",
            clean_line
        )
        if constructor_match:
            params = self._parse_parameters(constructor_match.group(1))
            return Method(
                name=class_name,
                return_type="void",
                parameters=params,
                is_constructor=True,
                nogil="nogil" in clean_line,
                except_clause="except" in clean_line,
            )

        # Check for destructor
        if clean_line.startswith("~"):
            return None  # Skip destructors

        # Check for operator
        operator_match = re.match(
            r"(\w+)\s+operator([^\(]+)\(([^)]*)\)\s*(except\s*\+)?\s*(nogil)?",
            clean_line
        )
        if operator_match:
            return_type = operator_match.group(1)
            op_name = operator_match.group(2).strip()
            params = self._parse_parameters(operator_match.group(3))
            return Method(
                name=f"operator{op_name}",
                return_type=return_type,
                parameters=params,
                is_operator=True,
                operator_name=op_name,
                nogil="nogil" in clean_line,
                except_clause="except" in clean_line,
            )

        # Regular method
        # Pattern: return_type method_name(params) except + nogil
        # Return type can be multi-word like "unsigned int", "const char*", etc.
        method_match = re.match(
            r"(.+?)\s+(\w+)\s*\(([^)]*)\)\s*(except\s*\+)?\s*(nogil)?",
            clean_line
        )
        if method_match:
            return_type = method_match.group(1)
            method_name = method_match.group(2)
            params_str = method_match.group(3)

            params = self._parse_parameters(params_str)

            return Method(
                name=method_name,
                return_type=return_type,
                parameters=params,
                nogil="nogil" in clean_line,
                except_clause="except" in clean_line,
            )

        # Void methods (no explicit return type)
        void_match = re.match(
            r"void\s+(\w+)\s*\(([^)]*)\)\s*(except\s*\+)?\s*(nogil)?",
            clean_line
        )
        if void_match:
            method_name = void_match.group(1)
            params_str = void_match.group(2)
            params = self._parse_parameters(params_str)

            return Method(
                name=method_name,
                return_type="void",
                parameters=params,
                nogil="nogil" in clean_line,
                except_clause="except" in clean_line,
            )

        return None

    def _parse_parameters(self, params_str: str) -> List[Parameter]:
        """Parse a parameter list string."""
        if not params_str.strip():
            return []

        params = []
        # Split by comma, but handle nested brackets
        depth = 0
        current = []
        for char in params_str:
            if char in "[<(":
                depth += 1
            elif char in "]>)":
                depth -= 1
            elif char == "," and depth == 0:
                if current:
                    params.append("".join(current).strip())
                current = []
                continue
            current.append(char)

        if current:
            params.append("".join(current).strip())

        result = []
        for param in params:
            param = param.strip()
            if not param:
                continue

            # Check for default value
            default = None
            if "=" in param:
                parts = param.split("=", 1)
                param = parts[0].strip()
                default = parts[1].strip()

            # Parse type and name
            # Pattern: type name or type[template] name
            # Handle reference (&) and pointer (*) that may be space-separated
            # Handle multi-word types like "unsigned int", "long long", etc.

            # Multi-word type prefixes that should not be treated as names
            type_prefixes = {
                'unsigned', 'signed', 'const', 'long', 'short', 'static',
                'volatile', 'mutable', 'constexpr'
            }

            parts = param.rsplit(None, 1)
            if len(parts) == 2:
                type_str = parts[0]
                name = parts[1]
                # If name is just & or *, it's part of the type
                if name in ("&", "*"):
                    type_str = param  # Use whole thing as type
                    name = f"arg{len(result)}"
                # If type_str is a type prefix, the whole thing is a type without name
                elif type_str in type_prefixes:
                    type_str = param  # Use whole thing as type (e.g., "unsigned int")
                    name = f"arg{len(result)}"
                # If name looks like a type modifier (int, long, etc.), it's part of the type
                elif name in ("int", "long", "short", "char", "double", "float"):
                    type_str = param  # Use whole thing as type
                    name = f"arg{len(result)}"
                # If type_str is 'const', the whole thing is a type (const Type&)
                elif type_str == "const":
                    type_str = param  # Use whole thing as type
                    name = f"arg{len(result)}"
                # If name looks like a type (contains & or * or is a known type)
                elif "&" in name or "*" in name:
                    type_str = param  # Use whole thing as type
                    name = f"arg{len(result)}"
            else:
                # Could be just a type with implicit name
                type_str = parts[0]
                name = f"arg{len(result)}"

            # Check for output parameter (reference to non-const)
            is_output = "&" in type_str and "const" not in type_str

            result.append(Parameter(
                name=name,
                type_str=type_str,
                default_value=default,
                is_output=is_output,
            ))

        return result

    def _infer_const(self, method: Method) -> bool:
        """
        Infer whether a method is const based on naming conventions.

        Heuristics:
        - Getters (get*, is*, has*, can*, should*, will*, empty, size, begin, end)
        - Methods returning const ref/ptr typically const
        - Methods taking no params (except self) and returning value typically const
        """
        name = method.name

        # Common const getter patterns
        const_prefixes = ('get', 'is', 'has', 'can', 'should', 'will', 'find', 'count')
        const_names = ('empty', 'size', 'length', 'begin', 'end', 'cbegin', 'cend',
                       'front', 'back', 'data', 'at', 'operator[]', 'min', 'max')

        # Check for explicit non-const patterns
        non_const_prefixes = ('set', 'add', 'remove', 'delete', 'clear', 'push', 'pop',
                              'insert', 'erase', 'swap', 'resize', 'reserve', 'assign',
                              'update', 'reset', 'sort', 'reverse', 'append', 'extend')

        # Non-const names
        if name.startswith(non_const_prefixes):
            return False

        # Const names
        if name.startswith(const_prefixes):
            return True
        if name in const_names:
            return True

        # Methods with no parameters that return a value are often const
        if len(method.parameters) == 0 and method.return_type != 'void':
            return True

        return False

    def _parse_wrap_directives(self, line: str) -> Dict[WrapDirective, str]:
        """Parse wrap- directives from a line."""
        directives = {}

        # Find all wrap-* patterns in the line
        for match in self._wrap_directive_pattern.finditer(line):
            directive_name = match.group(1).upper().replace("-", "_")
            directive_value = match.group(2).strip() if match.group(2) else ""

            try:
                directive = WrapDirective[directive_name]
                directives[directive] = directive_value
            except KeyError:
                pass  # Unknown directive, ignore

        # Special handling for doc: continuation
        if "#  " in line and WrapDirective.DOC not in directives:
            # This might be a doc continuation line
            doc_match = re.search(r"#\s{2,}(.+)$", line)
            if doc_match:
                directives[WrapDirective.DOC] = doc_match.group(1)

        return directives

    def _parse_instances(self, instances_str: str) -> Dict[str, List[str]]:
        """Parse wrap-instances directive value."""
        result = {}

        # Format: Name1[T1, T2], Name2[T3]
        # or: Name:T1,T2
        parts = instances_str.split(",")

        for part in parts:
            part = part.strip()
            if "[" in part:
                match = re.match(r"(\w+)\[([^\]]+)\]", part)
                if match:
                    name = match.group(1)
                    args = [a.strip() for a in match.group(2).split(",")]
                    result[name] = args
            elif ":" in part:
                name, args_str = part.split(":", 1)
                args = [a.strip() for a in args_str.split(",")]
                result[name.strip()] = args

        return result

    def _get_indent(self, line: str) -> int:
        """Get the indentation level of a line."""
        return len(line) - len(line.lstrip())

    def resolve_types(self, classes: Dict[str, ClassDecl]) -> None:
        """
        Resolve all type strings in class declarations to TypeInfo.

        Parameters
        ----------
        classes : dict
            Dictionary of class declarations to process.
        """
        for class_decl in classes.values():
            # Resolve method return types and parameters
            for method in class_decl.methods + class_decl.constructors:
                method.return_type_info = self.type_registry.resolve_type(
                    method.return_type
                )
                for param in method.parameters:
                    param.type_info = self.type_registry.resolve_type(param.type_str)
