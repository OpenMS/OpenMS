"""
Type Registry for pyOpenMS nanobind bindings

This module provides a registry of type mappings and conversion information
for generating nanobind bindings from pyOpenMS .pxd declarations.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Dict, List, Optional, Set


class TypeCategory(Enum):
    """Categories of types for binding generation."""

    PRIMITIVE = auto()  # int, double, bool, etc.
    STRING = auto()  # OpenMS::String, std::string
    OPAQUE = auto()  # Types passed by reference without conversion
    CONTAINER = auto()  # vector, map, set, etc.
    ENUM = auto()  # Enum types
    CLASS = auto()  # OpenMS classes
    CUSTOM_CASTER = auto()  # Types with custom nanobind type casters


@dataclass
class TypeInfo:
    """Information about a type for binding generation."""

    cpp_type: str
    python_type: str
    category: TypeCategory
    needs_custom_caster: bool = False
    caster_header: Optional[str] = None
    is_template: bool = False
    template_args: List[str] = field(default_factory=list)
    is_const: bool = False
    is_reference: bool = False
    is_pointer: bool = False
    is_shared_ptr: bool = False
    namespace: str = ""


class TypeRegistry:
    """
    Registry of type information for nanobind binding generation.

    This class provides mappings from Cython/C++ type declarations to
    information needed for generating nanobind bindings.
    """

    def __init__(self):
        self._types: Dict[str, TypeInfo] = {}
        self._aliases: Dict[str, str] = {}
        self._known_classes: Set[str] = set()
        self._enums: Set[str] = set()

        # Initialize built-in type mappings
        self._init_primitives()
        self._init_string_types()
        self._init_container_templates()
        self._init_openms_types()

    def _init_primitives(self):
        """Initialize primitive type mappings."""
        primitives = {
            # Integer types
            "int": ("int", "int"),
            "long": ("long", "int"),
            "long long": ("long long", "int"),
            "short": ("short", "int"),
            "char": ("char", "int"),
            "unsigned int": ("unsigned int", "int"),
            "unsigned long": ("unsigned long", "int"),
            "unsigned long long": ("unsigned long long", "int"),
            "unsigned short": ("unsigned short", "int"),
            "unsigned char": ("unsigned char", "int"),
            "size_t": ("size_t", "int"),
            "Size": ("size_t", "int"),
            "UInt": ("unsigned int", "int"),
            "UInt32": ("uint32_t", "int"),
            "UInt64": ("uint64_t", "int"),
            "Int": ("int", "int"),
            "Int32": ("int32_t", "int"),
            "Int64": ("int64_t", "int"),
            "SignedSize": ("ptrdiff_t", "int"),
            # Float types
            "float": ("float", "float"),
            "double": ("double", "float"),
            # Boolean
            "bool": ("bool", "bool"),
            "bint": ("bool", "bool"),  # Cython boolean
            # Void
            "void": ("void", "None"),
        }

        for cython_type, (cpp_type, py_type) in primitives.items():
            self._types[cython_type] = TypeInfo(
                cpp_type=cpp_type,
                python_type=py_type,
                category=TypeCategory.PRIMITIVE,
            )
            self._aliases[cython_type] = cython_type

    def _init_string_types(self):
        """Initialize string type mappings."""
        # OpenMS::String
        self._types["String"] = TypeInfo(
            cpp_type="OpenMS::String",
            python_type="str",
            category=TypeCategory.STRING,
            needs_custom_caster=True,
            caster_header="string.h",
            namespace="OpenMS",
        )

        # std::string
        self._types["std_string"] = TypeInfo(
            cpp_type="std::string",
            python_type="str",
            category=TypeCategory.STRING,
        )

        # Aliases
        self._aliases["String"] = "String"
        self._aliases["OpenMS::String"] = "String"
        self._aliases["libcpp_string"] = "std_string"
        self._aliases["libcpp_utf8_string"] = "String"
        self._aliases["string"] = "std_string"

    def _init_container_templates(self):
        """Initialize container template type patterns."""
        # These are patterns, not concrete types
        # The actual type resolution happens in resolve_type()
        pass

    def _init_openms_types(self):
        """Initialize OpenMS-specific types with custom casters."""
        # DPosition types
        self._types["DPosition1"] = TypeInfo(
            cpp_type="OpenMS::DPosition<1>",
            python_type="float",
            category=TypeCategory.CUSTOM_CASTER,
            needs_custom_caster=True,
            caster_header="dposition.h",
            is_template=True,
            template_args=["1"],
            namespace="OpenMS",
        )

        self._types["DPosition2"] = TypeInfo(
            cpp_type="OpenMS::DPosition<2>",
            python_type="Tuple[float, float]",
            category=TypeCategory.CUSTOM_CASTER,
            needs_custom_caster=True,
            caster_header="dposition.h",
            is_template=True,
            template_args=["2"],
            namespace="OpenMS",
        )

        # DataValue
        self._types["DataValue"] = TypeInfo(
            cpp_type="OpenMS::DataValue",
            python_type="Union[None, int, float, str, List[str], List[int], List[float]]",
            category=TypeCategory.CUSTOM_CASTER,
            needs_custom_caster=True,
            caster_header="datavalue.h",
            namespace="OpenMS",
        )

        # ParamValue
        self._types["ParamValue"] = TypeInfo(
            cpp_type="OpenMS::ParamValue",
            python_type="Union[None, int, float, str, List[str], List[int], List[float]]",
            category=TypeCategory.CUSTOM_CASTER,
            needs_custom_caster=True,
            caster_header="datavalue.h",
            namespace="OpenMS",
        )

        # List types
        self._types["StringList"] = TypeInfo(
            cpp_type="OpenMS::StringList",
            python_type="List[str]",
            category=TypeCategory.CUSTOM_CASTER,
            needs_custom_caster=True,
            caster_header="stl_containers.h",
            namespace="OpenMS",
        )

        self._types["IntList"] = TypeInfo(
            cpp_type="OpenMS::IntList",
            python_type="List[int]",
            category=TypeCategory.CONTAINER,
            namespace="OpenMS",
        )

        self._types["DoubleList"] = TypeInfo(
            cpp_type="OpenMS::DoubleList",
            python_type="List[float]",
            category=TypeCategory.CONTAINER,
            namespace="OpenMS",
        )

        # Aliases
        self._aliases["DPosition[1]"] = "DPosition1"
        self._aliases["DPosition[2]"] = "DPosition2"
        self._aliases["_DPosition1"] = "DPosition1"
        self._aliases["_DPosition2"] = "DPosition2"
        self._aliases["OpenMS::DPosition<1>"] = "DPosition1"
        self._aliases["OpenMS::DPosition<2>"] = "DPosition2"
        self._aliases["DataValue"] = "DataValue"
        self._aliases["OpenMS::DataValue"] = "DataValue"
        self._aliases["ParamValue"] = "ParamValue"
        self._aliases["OpenMS::ParamValue"] = "ParamValue"

    def register_class(self, class_name: str, namespace: str = "OpenMS") -> None:
        """
        Register a class that will be wrapped.

        Parameters
        ----------
        class_name : str
            Name of the class (without namespace).
        namespace : str
            C++ namespace (default: "OpenMS").
        """
        full_name = f"{namespace}::{class_name}" if namespace else class_name
        self._known_classes.add(class_name)
        self._known_classes.add(full_name)

        if class_name not in self._types:
            self._types[class_name] = TypeInfo(
                cpp_type=full_name,
                python_type=class_name,
                category=TypeCategory.CLASS,
                namespace=namespace,
            )
            self._aliases[class_name] = class_name
            self._aliases[full_name] = class_name

    def register_enum(self, enum_name: str, namespace: str = "OpenMS") -> None:
        """
        Register an enum type.

        Parameters
        ----------
        enum_name : str
            Name of the enum.
        namespace : str
            C++ namespace.
        """
        full_name = f"{namespace}::{enum_name}" if namespace else enum_name
        self._enums.add(enum_name)
        self._enums.add(full_name)

        self._types[enum_name] = TypeInfo(
            cpp_type=full_name,
            python_type=enum_name,
            category=TypeCategory.ENUM,
            namespace=namespace,
        )
        self._aliases[enum_name] = enum_name
        self._aliases[full_name] = enum_name

    def resolve_type(self, type_str: str) -> TypeInfo:
        """
        Resolve a type string to TypeInfo.

        Parameters
        ----------
        type_str : str
            Type string from .pxd file (e.g., "libcpp_vector[String]").

        Returns
        -------
        TypeInfo
            Resolved type information.
        """
        # Clean up the type string
        type_str = type_str.strip()

        # Check for const qualifier
        is_const = type_str.startswith("const ")
        if is_const:
            type_str = type_str[6:].strip()

        # Check for reference
        is_reference = type_str.endswith("&")
        if is_reference:
            type_str = type_str[:-1].strip()

        # Check for pointer
        is_pointer = type_str.endswith("*")
        if is_pointer:
            type_str = type_str[:-1].strip()

        # Check for shared_ptr
        is_shared_ptr = False
        shared_ptr_match = re.match(r"shared_ptr\[(.+)\]$", type_str)
        if shared_ptr_match:
            is_shared_ptr = True
            type_str = shared_ptr_match.group(1).strip()

        # Check aliases first
        if type_str in self._aliases:
            base_type = self._types.get(self._aliases[type_str])
            if base_type:
                return TypeInfo(
                    cpp_type=base_type.cpp_type,
                    python_type=base_type.python_type,
                    category=base_type.category,
                    needs_custom_caster=base_type.needs_custom_caster,
                    caster_header=base_type.caster_header,
                    is_template=base_type.is_template,
                    template_args=base_type.template_args.copy(),
                    is_const=is_const,
                    is_reference=is_reference,
                    is_pointer=is_pointer,
                    is_shared_ptr=is_shared_ptr,
                    namespace=base_type.namespace,
                )

        # Check for direct type match
        if type_str in self._types:
            base_type = self._types[type_str]
            return TypeInfo(
                cpp_type=base_type.cpp_type,
                python_type=base_type.python_type,
                category=base_type.category,
                needs_custom_caster=base_type.needs_custom_caster,
                caster_header=base_type.caster_header,
                is_template=base_type.is_template,
                template_args=base_type.template_args.copy(),
                is_const=is_const,
                is_reference=is_reference,
                is_pointer=is_pointer,
                is_shared_ptr=is_shared_ptr,
                namespace=base_type.namespace,
            )

        # Parse container types
        container_info = self._parse_container_type(type_str)
        if container_info:
            container_info.is_const = is_const
            container_info.is_reference = is_reference
            container_info.is_pointer = is_pointer
            container_info.is_shared_ptr = is_shared_ptr
            return container_info

        # Check if it's a known class
        if type_str in self._known_classes or type_str.replace("OpenMS::", "") in self._known_classes:
            clean_name = type_str.replace("OpenMS::", "")
            return TypeInfo(
                cpp_type=f"OpenMS::{clean_name}",
                python_type=clean_name,
                category=TypeCategory.CLASS,
                is_const=is_const,
                is_reference=is_reference,
                is_pointer=is_pointer,
                is_shared_ptr=is_shared_ptr,
                namespace="OpenMS",
            )

        # Unknown type - treat as opaque class
        return TypeInfo(
            cpp_type=type_str,
            python_type=type_str,
            category=TypeCategory.OPAQUE,
            is_const=is_const,
            is_reference=is_reference,
            is_pointer=is_pointer,
            is_shared_ptr=is_shared_ptr,
        )

    def _parse_container_type(self, type_str: str) -> Optional[TypeInfo]:
        """Parse container type strings like libcpp_vector[T], libcpp_map[K, V]."""
        # Vector patterns
        vector_match = re.match(r"(?:libcpp_)?vector\[(.+)\]$", type_str)
        if vector_match:
            inner = vector_match.group(1).strip()
            inner_type = self.resolve_type(inner)

            # Check for special vector types
            if inner_type.cpp_type == "OpenMS::String":
                return TypeInfo(
                    cpp_type="std::vector<OpenMS::String>",
                    python_type="List[str]",
                    category=TypeCategory.CUSTOM_CASTER,
                    needs_custom_caster=True,
                    caster_header="stl_containers.h",
                    is_template=True,
                    template_args=[inner],
                )
            elif inner_type.cpp_type == "OpenMS::DPosition<2>":
                return TypeInfo(
                    cpp_type="std::vector<OpenMS::DPosition<2>>",
                    python_type="numpy.ndarray",
                    category=TypeCategory.CUSTOM_CASTER,
                    needs_custom_caster=True,
                    caster_header="dposition.h",
                    is_template=True,
                    template_args=[inner],
                )
            else:
                return TypeInfo(
                    cpp_type=f"std::vector<{inner_type.cpp_type}>",
                    python_type=f"List[{inner_type.python_type}]",
                    category=TypeCategory.CONTAINER,
                    is_template=True,
                    template_args=[inner],
                )

        # Map patterns
        map_match = re.match(r"(?:libcpp_)?map\[(.+),\s*(.+)\]$", type_str)
        if map_match:
            key_str = map_match.group(1).strip()
            val_str = map_match.group(2).strip()
            key_type = self.resolve_type(key_str)
            val_type = self.resolve_type(val_str)

            needs_caster = key_type.cpp_type == "OpenMS::String"

            return TypeInfo(
                cpp_type=f"std::map<{key_type.cpp_type}, {val_type.cpp_type}>",
                python_type=f"Dict[{key_type.python_type}, {val_type.python_type}]",
                category=TypeCategory.CUSTOM_CASTER if needs_caster else TypeCategory.CONTAINER,
                needs_custom_caster=needs_caster,
                caster_header="stl_containers.h" if needs_caster else None,
                is_template=True,
                template_args=[key_str, val_str],
            )

        # Set patterns
        set_match = re.match(r"(?:libcpp_)?set\[(.+)\]$", type_str)
        if set_match:
            inner = set_match.group(1).strip()
            inner_type = self.resolve_type(inner)

            needs_caster = inner_type.cpp_type == "OpenMS::String"

            return TypeInfo(
                cpp_type=f"std::set<{inner_type.cpp_type}>",
                python_type=f"Set[{inner_type.python_type}]",
                category=TypeCategory.CUSTOM_CASTER if needs_caster else TypeCategory.CONTAINER,
                needs_custom_caster=needs_caster,
                caster_header="stl_containers.h" if needs_caster else None,
                is_template=True,
                template_args=[inner],
            )

        # Pair patterns
        pair_match = re.match(r"(?:libcpp_)?pair\[(.+),\s*(.+)\]$", type_str)
        if pair_match:
            first_str = pair_match.group(1).strip()
            second_str = pair_match.group(2).strip()
            first_type = self.resolve_type(first_str)
            second_type = self.resolve_type(second_str)

            return TypeInfo(
                cpp_type=f"std::pair<{first_type.cpp_type}, {second_type.cpp_type}>",
                python_type=f"Tuple[{first_type.python_type}, {second_type.python_type}]",
                category=TypeCategory.CONTAINER,
                is_template=True,
                template_args=[first_str, second_str],
            )

        return None

    def get_required_caster_headers(self, types: List[TypeInfo]) -> Set[str]:
        """
        Get the set of type caster headers required for the given types.

        Parameters
        ----------
        types : list of TypeInfo
            Types to check.

        Returns
        -------
        set of str
            Required header file names.
        """
        headers = set()
        for t in types:
            if t.needs_custom_caster and t.caster_header:
                headers.add(t.caster_header)
        return headers

    def get_cpp_type(self, type_info: TypeInfo) -> str:
        """
        Get the full C++ type declaration including qualifiers.

        Parameters
        ----------
        type_info : TypeInfo
            The type information.

        Returns
        -------
        str
            Full C++ type string.
        """
        result = type_info.cpp_type

        if type_info.is_shared_ptr:
            result = f"std::shared_ptr<{result}>"

        if type_info.is_const:
            result = f"const {result}"

        if type_info.is_reference:
            result = f"{result}&"
        elif type_info.is_pointer:
            result = f"{result}*"

        return result

    def is_known_class(self, name: str) -> bool:
        """Check if a class name is registered."""
        return name in self._known_classes

    def is_known_enum(self, name: str) -> bool:
        """Check if an enum name is registered."""
        return name in self._enums
