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
from typing import Dict, List, Optional, Set, Tuple, Union

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


@dataclass
class EnumValue:
    """An enum value."""

    name: str
    value: Optional[int] = None
    doc: str = ""


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
            r'cdef\s+cppclass\s+(\w+)(?:\s+"([^"]+)")?(?:\[([^\]]+)\])?(?:\(([^)]+)\))?'
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
            r"#\s*wrap-(\w+)(?::\s*(.+))?$"
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
                    template_args, base_classes
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
    ) -> Tuple[Optional[ClassDecl], int]:
        """Parse the body of a class declaration."""
        self._current_class = class_name

        class_decl = ClassDecl(
            name=class_name,
            namespace=namespace,
            header_file=header_file,
            template_args=template_args,
            base_classes=base_classes,
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
                if "ABSTRACT" in stripped.upper() and "class" in stripped.lower():
                    class_decl.is_abstract = True
                    logger.debug(f"Class {class_name} marked as abstract from pxd comment")
                directives = self._parse_wrap_directives(stripped)
                if WrapDirective.DOC in directives:
                    pending_doc.append(directives[WrapDirective.DOC])
                i += 1
                continue

            # Parse method/constructor
            method = self._parse_method_line(stripped, class_name)
            if method:
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

        while i < len(lines):
            line = lines[i]
            current_indent = self._get_indent(line)
            stripped = line.strip()

            if stripped and current_indent <= base_indent and not stripped.startswith("#"):
                break

            if not stripped or stripped.startswith("#"):
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

                if "=" in val:
                    parts = val.split("=")
                    name = parts[0].strip()
                    value = int(parts[1].strip())
                    enum_decl.values.append(EnumValue(name=name, value=value))
                else:
                    enum_decl.values.append(EnumValue(name=val))

            i += 1

        return enum_decl, i

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
