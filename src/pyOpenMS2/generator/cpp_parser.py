"""
C++ header parser using libclang.

This module provides accurate C++ type information by parsing headers
directly, avoiding the limitations of .pxd file parsing.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set

logger = logging.getLogger(__name__)

# Try to import clang bindings
try:
    from clang.cindex import (
        Index,
        CursorKind,
        TypeKind,
        AccessSpecifier,
        TranslationUnit,
        Config,
    )

    # Configure libclang library path if not found automatically
    import ctypes.util
    if not ctypes.util.find_library("clang"):
        # Try common paths - prefer versioned libclang (not libclang-cpp)
        search_paths = [
            "/usr/lib/x86_64-linux-gnu",
            "/usr/lib/llvm-14/lib",
            "/usr/lib/llvm-15/lib",
            "/usr/lib/llvm-16/lib",
            "/usr/lib/llvm-17/lib",
        ]
        for lib_path in search_paths:
            lib_dir = Path(lib_path)
            if not lib_dir.exists():
                continue
            # Look for libclang-XX.so.1 (not libclang-cpp)
            for f in sorted(lib_dir.glob("libclang-[0-9]*.so.[0-9]*"), reverse=True):
                if "cpp" not in f.name and f.is_file():
                    Config.set_library_file(str(f))
                    break
            else:
                # Try libclang.so
                clang_lib = lib_dir / "libclang.so"
                if clang_lib.exists():
                    Config.set_library_file(str(clang_lib))
                continue
            break

    CLANG_AVAILABLE = True
except ImportError:
    CLANG_AVAILABLE = False
    logger.warning("clang Python bindings not available. Install with: pip install clang")


@dataclass
class CppParameter:
    """A C++ function parameter."""
    name: str
    type_str: str  # Full C++ type string
    is_const: bool = False
    is_reference: bool = False
    is_pointer: bool = False
    default_value: Optional[str] = None


@dataclass
class CppMethod:
    """A C++ class method."""
    name: str
    return_type: str
    parameters: List[CppParameter] = field(default_factory=list)
    is_const: bool = False
    is_static: bool = False
    is_virtual: bool = False
    is_pure_virtual: bool = False
    is_constructor: bool = False
    is_destructor: bool = False
    access: str = "public"  # public, protected, private


@dataclass
class CppClass:
    """A C++ class declaration."""
    name: str
    qualified_name: str  # Full namespace::class name
    namespace: str
    header_file: str
    base_classes: List[str] = field(default_factory=list)
    methods: List[CppMethod] = field(default_factory=list)
    constructors: List[CppMethod] = field(default_factory=list)
    is_abstract: bool = False
    is_template: bool = False
    template_params: List[str] = field(default_factory=list)
    doc: str = ""

    def find_method(self, name: str, param_count: Optional[int] = None) -> Optional[CppMethod]:
        """
        Find a method by name, optionally filtering by parameter count.

        Parameters
        ----------
        name : str
            Method name.
        param_count : int, optional
            Expected number of parameters.

        Returns
        -------
        CppMethod or None
            The matching method, or None if not found.
        """
        candidates = [m for m in self.methods if m.name == name]
        if not candidates:
            return None
        if param_count is not None:
            candidates = [m for m in candidates if len(m.parameters) == param_count]
        if len(candidates) == 1:
            return candidates[0]
        # Multiple matches - return the first (caller should use get_overloads)
        return candidates[0] if candidates else None

    def get_overloads(self, name: str) -> List[CppMethod]:
        """Get all overloads of a method by name."""
        return [m for m in self.methods if m.name == name]

    def has_overloads(self, name: str) -> bool:
        """Check if a method has multiple overloads."""
        return len(self.get_overloads(name)) > 1


class CppHeaderParser:
    """
    Parse C++ headers using libclang to extract accurate type information.
    """

    def __init__(self, include_paths: List[Path] = None, compiler_args: List[str] = None):
        """
        Initialize the parser.

        Parameters
        ----------
        include_paths : list of Path
            Additional include directories for header resolution.
        compiler_args : list of str
            Additional compiler arguments (e.g., -std=c++17).
        """
        if not CLANG_AVAILABLE:
            raise RuntimeError("clang Python bindings not available")

        self.index = Index.create()
        self.include_paths = include_paths or []
        self.compiler_args = compiler_args or ["-std=c++17", "-x", "c++"]

        # Build include path arguments
        for path in self.include_paths:
            self.compiler_args.append(f"-I{path}")

        self._classes: Dict[str, CppClass] = {}

    def parse_header(self, header_path: Path) -> Dict[str, CppClass]:
        """
        Parse a single header file.

        Parameters
        ----------
        header_path : Path
            Path to the header file.

        Returns
        -------
        dict
            Dictionary mapping class names to CppClass objects.
        """
        logger.debug(f"Parsing {header_path}")

        tu = self.index.parse(
            str(header_path),
            args=self.compiler_args,
            options=TranslationUnit.PARSE_SKIP_FUNCTION_BODIES,
        )

        # Check for errors
        for diag in tu.diagnostics:
            if diag.severity >= 3:  # Error or fatal
                logger.warning(f"Parse error in {header_path}: {diag.spelling}")

        # Walk the AST
        self._walk_cursor(tu.cursor, str(header_path))

        return self._classes

    def parse_headers(self, header_paths: List[Path]) -> Dict[str, CppClass]:
        """
        Parse multiple header files.

        Parameters
        ----------
        header_paths : list of Path
            Paths to header files.

        Returns
        -------
        dict
            Dictionary mapping class names to CppClass objects.
        """
        for path in header_paths:
            try:
                self.parse_header(path)
            except Exception as e:
                logger.warning(f"Failed to parse {path}: {e}")

        return self._classes

    def _walk_cursor(self, cursor, source_file: str, namespace: str = ""):
        """Recursively walk the AST cursor."""
        for child in cursor.get_children():
            # Only process nodes from our source file
            if child.location.file and str(child.location.file) != source_file:
                continue

            if child.kind == CursorKind.NAMESPACE:
                ns_name = child.spelling
                new_namespace = f"{namespace}::{ns_name}" if namespace else ns_name
                self._walk_cursor(child, source_file, new_namespace)

            elif child.kind in (CursorKind.CLASS_DECL, CursorKind.STRUCT_DECL):
                self._process_class(child, source_file, namespace)

            elif child.kind == CursorKind.CLASS_TEMPLATE:
                self._process_class_template(child, source_file, namespace)

    def _process_class(self, cursor, source_file: str, namespace: str):
        """Process a class declaration."""
        class_name = cursor.spelling
        if not class_name:
            return

        qualified_name = f"{namespace}::{class_name}" if namespace else class_name

        cpp_class = CppClass(
            name=class_name,
            qualified_name=qualified_name,
            namespace=namespace,
            header_file=source_file,
        )

        # Get base classes
        for child in cursor.get_children():
            if child.kind == CursorKind.CXX_BASE_SPECIFIER:
                base_type = child.type.spelling
                cpp_class.base_classes.append(base_type)

        # Process members
        self._process_class_members(cursor, cpp_class)

        # Check if abstract (has any pure virtual methods)
        cpp_class.is_abstract = any(m.is_pure_virtual for m in cpp_class.methods)

        self._classes[class_name] = cpp_class
        logger.debug(f"Parsed class: {qualified_name} ({len(cpp_class.methods)} methods)")

    def _process_class_template(self, cursor, source_file: str, namespace: str):
        """Process a class template declaration."""
        class_name = cursor.spelling
        if not class_name:
            return

        qualified_name = f"{namespace}::{class_name}" if namespace else class_name

        cpp_class = CppClass(
            name=class_name,
            qualified_name=qualified_name,
            namespace=namespace,
            header_file=source_file,
            is_template=True,
        )

        # Get template parameters
        for child in cursor.get_children():
            if child.kind == CursorKind.TEMPLATE_TYPE_PARAMETER:
                cpp_class.template_params.append(child.spelling)
            elif child.kind == CursorKind.TEMPLATE_NON_TYPE_PARAMETER:
                cpp_class.template_params.append(child.spelling)

        # Process the actual class within the template
        for child in cursor.get_children():
            if child.kind in (CursorKind.CLASS_DECL, CursorKind.STRUCT_DECL):
                # Get base classes
                for subchild in child.get_children():
                    if subchild.kind == CursorKind.CXX_BASE_SPECIFIER:
                        cpp_class.base_classes.append(subchild.type.spelling)

                self._process_class_members(child, cpp_class)
                break

        cpp_class.is_abstract = any(m.is_pure_virtual for m in cpp_class.methods)
        self._classes[class_name] = cpp_class

    def _process_class_members(self, cursor, cpp_class: CppClass):
        """Process class members (methods, constructors)."""
        current_access = "private"  # Default for classes

        for child in cursor.get_children():
            # Track access specifier
            if child.kind == CursorKind.CXX_ACCESS_SPEC_DECL:
                if child.access_specifier == AccessSpecifier.PUBLIC:
                    current_access = "public"
                elif child.access_specifier == AccessSpecifier.PROTECTED:
                    current_access = "protected"
                else:
                    current_access = "private"
                continue

            # Skip non-public members for binding
            if current_access != "public":
                continue

            if child.kind == CursorKind.CONSTRUCTOR:
                ctor = self._process_method(child, is_constructor=True)
                ctor.access = current_access
                cpp_class.constructors.append(ctor)

            elif child.kind == CursorKind.DESTRUCTOR:
                # Skip destructors
                pass

            elif child.kind == CursorKind.CXX_METHOD:
                method = self._process_method(child)
                method.access = current_access
                cpp_class.methods.append(method)

    def _process_method(self, cursor, is_constructor: bool = False) -> CppMethod:
        """Process a method declaration."""
        method = CppMethod(
            name=cursor.spelling,
            return_type=cursor.result_type.spelling if not is_constructor else "",
            is_constructor=is_constructor,
            is_const=cursor.is_const_method(),
            is_static=cursor.is_static_method(),
            is_virtual=cursor.is_virtual_method(),
            is_pure_virtual=cursor.is_pure_virtual_method(),
        )

        # Process parameters
        for child in cursor.get_children():
            if child.kind == CursorKind.PARM_DECL:
                param = self._process_parameter(child)
                method.parameters.append(param)

        return method

    def _process_parameter(self, cursor) -> CppParameter:
        """Process a parameter declaration."""
        param_type = cursor.type
        type_str = param_type.spelling

        return CppParameter(
            name=cursor.spelling or f"arg{cursor.hash}",
            type_str=type_str,
            is_const="const" in type_str,
            is_reference=param_type.kind == TypeKind.LVALUEREFERENCE,
            is_pointer=param_type.kind == TypeKind.POINTER,
        )


@dataclass
class MergedMethod:
    """Method info merged from C++ AST and .pxd declarations."""
    cpp_method: CppMethod
    wrap_as: Optional[str] = None  # Python name if different from C++ name
    wrap_ignore: bool = False
    doc: str = ""


@dataclass
class MergedClass:
    """Class info merged from C++ AST and .pxd declarations."""
    cpp_class: CppClass
    methods: List[MergedMethod]
    constructors: List[CppMethod]
    wrap_hash: bool = False
    wrap_iter: Optional[str] = None
    wrap_manual_memory: bool = False
    doc: str = ""

    @property
    def name(self) -> str:
        return self.cpp_class.name

    @property
    def qualified_name(self) -> str:
        return self.cpp_class.qualified_name

    @property
    def namespace(self) -> str:
        return self.cpp_class.namespace

    @property
    def header_file(self) -> str:
        return self.cpp_class.header_file

    @property
    def base_classes(self) -> List[str]:
        return self.cpp_class.base_classes

    @property
    def is_abstract(self) -> bool:
        return self.cpp_class.is_abstract


def merge_with_pxd(
    cpp_classes: Dict[str, CppClass],
    pxd_classes: Dict[str, "ClassDecl"],  # From pxd_parser
) -> Dict[str, MergedClass]:
    """
    Merge C++ AST information with .pxd declarations.

    The .pxd files serve as an allowlist - only classes/methods declared
    in .pxd files will be wrapped, but with accurate C++ signatures.

    Parameters
    ----------
    cpp_classes : dict
        Classes parsed from C++ headers.
    pxd_classes : dict
        Classes declared in .pxd files.

    Returns
    -------
    dict
        Merged class information with accurate C++ types.
    """
    merged = {}

    for class_name, pxd_class in pxd_classes.items():
        if class_name not in cpp_classes:
            logger.warning(f"Class {class_name} in .pxd but not found in C++ headers")
            continue

        cpp_class = cpp_classes[class_name]

        # Build method allowlist from .pxd with directives
        pxd_methods_by_name = {}
        for m in pxd_class.methods:
            if m.name not in pxd_methods_by_name:
                pxd_methods_by_name[m.name] = []
            pxd_methods_by_name[m.name].append(m)

        # Match C++ methods with .pxd declarations
        merged_methods = []
        for cpp_method in cpp_class.methods:
            if cpp_method.name not in pxd_methods_by_name:
                continue
            if cpp_method.access != "public":
                continue

            # Find the corresponding .pxd declaration
            pxd_overloads = pxd_methods_by_name[cpp_method.name]

            # Try to match by parameter count
            matching_pxd = None
            for pxd_m in pxd_overloads:
                if len(pxd_m.parameters) == len(cpp_method.parameters):
                    matching_pxd = pxd_m
                    break

            # If no exact match, use the first one
            if matching_pxd is None and pxd_overloads:
                matching_pxd = pxd_overloads[0]

            merged_method = MergedMethod(
                cpp_method=cpp_method,
                wrap_as=matching_pxd.wrap_as if matching_pxd else None,
                wrap_ignore=matching_pxd.wrap_ignore if matching_pxd else False,
                doc=matching_pxd.doc if matching_pxd else "",
            )
            merged_methods.append(merged_method)

        # Filter constructors
        constructors = [c for c in cpp_class.constructors if c.access == "public"]

        merged_class = MergedClass(
            cpp_class=cpp_class,
            methods=merged_methods,
            constructors=constructors,
            wrap_hash=getattr(pxd_class, 'wrap_hash', False),
            wrap_iter=getattr(pxd_class, 'wrap_iter', None),
            wrap_manual_memory=getattr(pxd_class, 'wrap_manual_memory', False),
            doc=getattr(pxd_class, 'doc', ''),
        )
        merged[class_name] = merged_class

    return merged


# Example usage and testing
if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.DEBUG)

    if not CLANG_AVAILABLE:
        print("clang bindings not available. Install with: pip install clang")
        sys.exit(1)

    # Test with a simple header
    test_header = Path("/tmp/test_header.h")
    test_header.write_text("""
namespace OpenMS {

class TestClass {
public:
    TestClass();
    TestClass(const TestClass& other);

    int getValue() const;
    void setValue(int value);

    static TestClass create();

    virtual void doSomething() = 0;

private:
    int value_;
};

}  // namespace OpenMS
""")

    parser = CppHeaderParser()
    classes = parser.parse_header(test_header)

    for name, cls in classes.items():
        print(f"\nClass: {cls.qualified_name}")
        print(f"  Abstract: {cls.is_abstract}")
        print(f"  Base classes: {cls.base_classes}")
        print(f"  Constructors: {len(cls.constructors)}")
        for ctor in cls.constructors:
            params = ", ".join(p.type_str for p in ctor.parameters)
            print(f"    {ctor.name}({params})")
        print(f"  Methods: {len(cls.methods)}")
        for method in cls.methods:
            params = ", ".join(p.type_str for p in method.parameters)
            const = " const" if method.is_const else ""
            static = "static " if method.is_static else ""
            virtual = "virtual " if method.is_virtual else ""
            pure = " = 0" if method.is_pure_virtual else ""
            print(f"    {static}{virtual}{method.return_type} {method.name}({params}){const}{pure}")
