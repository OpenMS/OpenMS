"""
C++ header parser using libclang.

This module provides accurate C++ type information by parsing headers
directly, avoiding the limitations of .pxd file parsing.

Supports caching parsed results to JSON for faster subsequent runs.
"""

import json
import hashlib
import logging
import os
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Any

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
    type_str: str  # Full C++ type string (original spelling)
    canonical_type: str = ""  # Canonical type with typedefs resolved
    is_const: bool = False
    is_reference: bool = False
    is_pointer: bool = False
    default_value: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary for JSON caching."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "CppParameter":
        """Deserialize from dictionary."""
        return cls(**d)


@dataclass
class CppMethod:
    """A C++ class method."""
    name: str
    return_type: str  # Original spelling
    canonical_return_type: str = ""  # Canonical type with typedefs resolved
    parameters: List[CppParameter] = field(default_factory=list)
    is_const: bool = False
    is_static: bool = False
    is_virtual: bool = False
    is_pure_virtual: bool = False
    is_constructor: bool = False
    is_destructor: bool = False
    access: str = "public"  # public, protected, private

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary for JSON caching."""
        d = asdict(self)
        d['parameters'] = [p.to_dict() if hasattr(p, 'to_dict') else p for p in self.parameters]
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "CppMethod":
        """Deserialize from dictionary."""
        d = d.copy()
        d['parameters'] = [CppParameter.from_dict(p) if isinstance(p, dict) else p for p in d.get('parameters', [])]
        return cls(**d)


@dataclass
class CppClass:
    """A C++ class declaration."""
    name: str
    qualified_name: str  # Full namespace::class name
    namespace: str
    header_file: str
    base_classes: List[str] = field(default_factory=list)  # Original spelling
    canonical_base_classes: List[str] = field(default_factory=list)  # Canonical types
    methods: List[CppMethod] = field(default_factory=list)
    constructors: List[CppMethod] = field(default_factory=list)
    is_abstract: bool = False
    is_template: bool = False
    template_params: List[str] = field(default_factory=list)
    doc: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary for JSON caching."""
        return {
            'name': self.name,
            'qualified_name': self.qualified_name,
            'namespace': self.namespace,
            'header_file': self.header_file,
            'base_classes': self.base_classes,
            'canonical_base_classes': self.canonical_base_classes,
            'methods': [m.to_dict() for m in self.methods],
            'constructors': [c.to_dict() for c in self.constructors],
            'is_abstract': self.is_abstract,
            'is_template': self.is_template,
            'template_params': self.template_params,
            'doc': self.doc,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "CppClass":
        """Deserialize from dictionary."""
        return cls(
            name=d['name'],
            qualified_name=d['qualified_name'],
            namespace=d['namespace'],
            header_file=d['header_file'],
            base_classes=d.get('base_classes', []),
            canonical_base_classes=d.get('canonical_base_classes', []),
            methods=[CppMethod.from_dict(m) for m in d.get('methods', [])],
            constructors=[CppMethod.from_dict(c) for c in d.get('constructors', [])],
            is_abstract=d.get('is_abstract', False),
            is_template=d.get('is_template', False),
            template_params=d.get('template_params', []),
            doc=d.get('doc', ''),
        )

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

    Supports caching parsed results to JSON for faster subsequent runs.
    """

    # Cache format version - bump this when dataclass structure changes
    CACHE_VERSION = 5  # Bumped for struct default public access fix

    def __init__(
        self,
        include_paths: List[Path] = None,
        compiler_args: List[str] = None,
        cache_dir: Optional[Path] = None,
    ):
        """
        Initialize the parser.

        Parameters
        ----------
        include_paths : list of Path
            Additional include directories for header resolution.
        compiler_args : list of str
            Additional compiler arguments (e.g., -std=c++17).
        cache_dir : Path, optional
            Directory to store cached parse results. If None, caching is disabled.
        """
        if not CLANG_AVAILABLE:
            raise RuntimeError("clang Python bindings not available")

        self.index = Index.create()
        self.include_paths = include_paths or []
        self.compiler_args = compiler_args or ["-std=c++20", "-x", "c++"]
        self.cache_dir = Path(cache_dir) if cache_dir else None

        # Build include path arguments
        for path in self.include_paths:
            self.compiler_args.append(f"-I{path}")

        # Add Qt include paths if available (needed for proper type resolution)
        qt_include_paths = [
            Path("/usr/include/x86_64-linux-gnu/qt6"),
            Path("/usr/include/x86_64-linux-gnu/qt6/QtCore"),
            Path("/usr/include/x86_64-linux-gnu/qt6/QtNetwork"),
            Path("/usr/include/qt6"),
            Path("/usr/include/qt6/QtCore"),
            Path("/usr/include/qt6/QtNetwork"),
        ]
        for qt_path in qt_include_paths:
            if qt_path.exists():
                self.compiler_args.append(f"-I{qt_path}")

        # Add OpenSwathAlgo include paths (relative to OpenMS includes)
        for inc_path in self.include_paths:
            # Try to find openswathalgo relative to openms
            openswath_path = inc_path.parent.parent / "openswathalgo" / "include"
            if openswath_path.exists():
                self.compiler_args.append(f"-I{openswath_path}")

        self._classes: Dict[str, CppClass] = {}

        # Create cache directory if specified
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"libclang cache directory: {self.cache_dir}")

    @staticmethod
    def _normalize_filename(filename: str) -> str:
        """Normalize a filename for stable comparisons across platforms."""
        return os.path.normpath(filename)

    def _get_cache_key(self, header_path: Path) -> str:
        """Generate a cache key for a header file based on path and mtime."""
        path_str = str(header_path.resolve())
        mtime = header_path.stat().st_mtime
        # Include compiler args in key for reproducibility
        args_str = "|".join(sorted(self.compiler_args))
        key_data = f"{path_str}|{mtime}|{args_str}|v{self.CACHE_VERSION}"
        return hashlib.md5(key_data.encode()).hexdigest()

    def _get_cache_path(self, header_path: Path) -> Path:
        """Get the cache file path for a header."""
        cache_key = self._get_cache_key(header_path)
        return self.cache_dir / f"{cache_key}.json"

    def _load_from_cache(self, header_path: Path) -> Optional[Dict[str, CppClass]]:
        """Load parsed classes from cache if valid."""
        if not self.cache_dir:
            return None

        cache_path = self._get_cache_path(header_path)
        if not cache_path.exists():
            return None

        try:
            with open(cache_path, 'r') as f:
                data = json.load(f)

            # Validate cache version
            if data.get('version') != self.CACHE_VERSION:
                logger.debug(f"Cache version mismatch for {header_path}")
                return None

            # Deserialize classes
            classes = {}
            for name, class_dict in data.get('classes', {}).items():
                classes[name] = CppClass.from_dict(class_dict)

            logger.debug(f"Loaded {len(classes)} classes from cache for {header_path.name}")
            return classes

        except (json.JSONDecodeError, KeyError, TypeError) as e:
            logger.debug(f"Cache load failed for {header_path}: {e}")
            return None

    def _save_to_cache(self, header_path: Path, classes: Dict[str, CppClass]):
        """Save parsed classes to cache."""
        if not self.cache_dir:
            return

        cache_path = self._get_cache_path(header_path)
        try:
            data = {
                'version': self.CACHE_VERSION,
                'header': str(header_path),
                'classes': {name: cls.to_dict() for name, cls in classes.items()},
            }
            with open(cache_path, 'w') as f:
                json.dump(data, f, indent=2)
            logger.debug(f"Saved {len(classes)} classes to cache for {header_path.name}")
        except Exception as e:
            logger.warning(f"Failed to save cache for {header_path}: {e}")

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
        header_path = Path(header_path).resolve()

        # Try loading from cache first
        cached = self._load_from_cache(header_path)
        if cached is not None:
            self._classes.update(cached)
            return self._classes

        logger.debug(f"Parsing {header_path} with libclang")

        tu = self.index.parse(
            str(header_path),
            args=self.compiler_args,
            options=TranslationUnit.PARSE_SKIP_FUNCTION_BODIES,
        )

        # Check for errors
        for diag in tu.diagnostics:
            if diag.severity >= 3:  # Error or fatal
                logger.warning(f"Parse error in {header_path}: {diag.spelling}")

        parsed_classes: Dict[str, CppClass] = {}
        allowed_files = {self._normalize_filename(str(header_path))}
        self._walk_cursor(tu.cursor, allowed_files, parsed_classes)
        if parsed_classes:
            self._classes.update(parsed_classes)
            self._save_to_cache(header_path, parsed_classes)

        return self._classes

    def parse_headers(
        self,
        header_paths: List[Path],
        *,
        batch_mode: str = "auto",
        batch_size: int = 50,
        auto_batch_min_headers: int = 25,
    ) -> Dict[str, CppClass]:
        """
        Parse multiple header files.

        Parameters
        ----------
        header_paths : list of Path
            Paths to header files.
        batch_mode : {"auto", "on", "off"}
            Whether to batch-parse multiple headers per translation unit to avoid
            repeatedly parsing common transitive includes.
        batch_size : int
            Number of headers per batch when batching is enabled.
        auto_batch_min_headers : int
            Minimum number of cache misses required to enable batching when
            batch_mode is "auto".

        Returns
        -------
        dict
            Dictionary mapping class names to CppClass objects.
        """
        batch_mode = (batch_mode or "auto").lower()
        if batch_mode not in {"auto", "on", "off"}:
            raise ValueError(f"Invalid batch_mode: {batch_mode!r} (expected 'auto', 'on', or 'off')")
        if batch_size <= 0:
            raise ValueError("batch_size must be > 0")
        if auto_batch_min_headers <= 0:
            raise ValueError("auto_batch_min_headers must be > 0")

        # Deduplicate while preserving order
        unique_headers: List[Path] = []
        seen: Set[str] = set()
        for path in header_paths:
            p = Path(path).resolve()
            key = self._normalize_filename(str(p))
            if key in seen:
                continue
            unique_headers.append(p)
            seen.add(key)

        cache_hits = 0
        cache_misses = 0
        to_parse: List[Path] = []

        for path in unique_headers:
            try:
                cached = self._load_from_cache(path)
                if cached is not None:
                    cache_hits += 1
                    self._classes.update(cached)
                else:
                    cache_misses += 1
                    to_parse.append(path)
            except Exception as e:
                logger.warning(f"Cache check failed for {path}: {e}")
                cache_misses += 1
                to_parse.append(path)

        use_batch = False
        if batch_mode == "on":
            use_batch = True
        elif batch_mode == "off":
            use_batch = False
        else:
            use_batch = len(to_parse) >= auto_batch_min_headers

        if use_batch and to_parse:
            for i in range(0, len(to_parse), batch_size):
                batch = to_parse[i:i + batch_size]
                try:
                    self._parse_headers_in_batch(batch)
                except Exception as e:
                    logger.warning(f"Failed to batch-parse {len(batch)} headers: {e}")
                    # Fall back to individual parsing for this batch
                    for path in batch:
                        try:
                            self.parse_header(path)
                        except Exception as inner_e:
                            logger.warning(f"Failed to parse {path}: {inner_e}")
        else:
            for path in to_parse:
                try:
                    self.parse_header(path)
                except Exception as e:
                    logger.warning(f"Failed to parse {path}: {e}")

        if self.cache_dir:
            logger.info(f"libclang cache: {cache_hits} hits, {cache_misses} misses")

        return self._classes

    def _parse_headers_in_batch(self, header_paths: List[Path]) -> Dict[str, CppClass]:
        """
        Parse a batch of headers in a single translation unit.

        This avoids repeatedly parsing common transitive includes (e.g., STL),
        which is a major bottleneck when parsing hundreds of headers one-by-one.
        """
        header_paths = [Path(p).resolve() for p in header_paths]
        allowed_files = {self._normalize_filename(str(p)) for p in header_paths}

        include_lines = "\n".join(f'#include "{p.as_posix()}"' for p in header_paths)
        batch_source = (
            "// Auto-generated by pyOpenMS2 CppHeaderParser for faster libclang parsing\n"
            f"{include_lines}\n"
        )

        base_dir = self.cache_dir or Path(tempfile.gettempdir())
        batch_path = str((base_dir / "__pyopenms2_libclang_batch.cpp").resolve())

        tu = self.index.parse(
            batch_path,
            args=self.compiler_args,
            unsaved_files=[(batch_path, batch_source)],
            options=TranslationUnit.PARSE_SKIP_FUNCTION_BODIES,
        )

        # Check for errors
        for diag in tu.diagnostics:
            if diag.severity >= 3:  # Error or fatal
                logger.warning(f"Parse error in batch: {diag.spelling}")

        parsed_classes: Dict[str, CppClass] = {}
        self._walk_cursor(tu.cursor, allowed_files, parsed_classes)
        if not parsed_classes:
            return self._classes

        self._classes.update(parsed_classes)

        # Save per-header caches
        classes_by_header: Dict[str, Dict[str, CppClass]] = {}
        for name, cls in parsed_classes.items():
            classes_by_header.setdefault(cls.header_file, {})[name] = cls

        for header_path in header_paths:
            header_key = self._normalize_filename(str(header_path))
            classes_for_header = classes_by_header.get(header_key)
            if classes_for_header:
                self._save_to_cache(header_path, classes_for_header)

        return self._classes

    def clear_cache(self):
        """Clear all cached parse results."""
        if self.cache_dir and self.cache_dir.exists():
            for cache_file in self.cache_dir.glob("*.json"):
                cache_file.unlink()
            logger.info("Cleared libclang cache")

    def _walk_cursor(
        self,
        cursor,
        allowed_files: Set[str],
        out_classes: Dict[str, CppClass],
        namespace: str = "",
    ):
        """Recursively walk the AST cursor, extracting declarations from allowed files."""
        for child in cursor.get_children():
            if not child.location.file:
                continue

            child_file = self._normalize_filename(str(child.location.file))
            if child_file not in allowed_files:
                continue

            if child.kind == CursorKind.NAMESPACE:
                ns_name = child.spelling
                new_namespace = f"{namespace}::{ns_name}" if namespace else ns_name
                self._walk_cursor(child, allowed_files, out_classes, new_namespace)

            elif child.kind in (CursorKind.CLASS_DECL, CursorKind.STRUCT_DECL):
                cpp_class = self._process_class(child, namespace, out_classes=out_classes)
                if cpp_class:
                    out_classes[cpp_class.name] = cpp_class

            elif child.kind == CursorKind.CLASS_TEMPLATE:
                cpp_class = self._process_class_template(child, namespace, out_classes=out_classes)
                if cpp_class:
                    out_classes[cpp_class.name] = cpp_class

    def _process_class(
        self,
        cursor,
        namespace: str,
        out_classes: Optional[Dict[str, CppClass]] = None,
        parent_class_name: str = "",
    ) -> Optional[CppClass]:
        """Process a class declaration.

        Parameters
        ----------
        cursor : clang.cindex.Cursor
            The class cursor to process.
        namespace : str
            The enclosing namespace.
        out_classes : dict, optional
            Dictionary to collect nested classes into. If provided, nested classes
            will be added with flattened names (e.g., ParentClass_NestedClass).
        parent_class_name : str
            For nested classes, the name of the parent class (used for name flattening).
        """
        class_name = cursor.spelling
        if not class_name:
            return None
        if not cursor.is_definition():
            return None

        # For nested classes, create a flattened name
        if parent_class_name:
            flattened_name = f"{parent_class_name}_{class_name}"
            # Qualified name includes the parent class
            qualified_name = f"{namespace}::{parent_class_name}::{class_name}" if namespace else f"{parent_class_name}::{class_name}"
        else:
            flattened_name = class_name
            qualified_name = f"{namespace}::{class_name}" if namespace else class_name

        cpp_class = CppClass(
            name=flattened_name,  # Use flattened name for matching with pxd
            qualified_name=qualified_name,
            namespace=namespace,
            header_file=self._normalize_filename(str(cursor.location.file)),
        )

        # Get base classes (both original and canonical spelling)
        for child in cursor.get_children():
            if child.kind == CursorKind.CXX_BASE_SPECIFIER:
                base_type = child.type
                cpp_class.base_classes.append(base_type.spelling)
                cpp_class.canonical_base_classes.append(base_type.get_canonical().spelling)

        # Process members (including nested classes)
        self._process_class_members(cursor, cpp_class, out_classes, flattened_name)

        # Check if abstract (has any pure virtual methods)
        cpp_class.is_abstract = any(m.is_pure_virtual for m in cpp_class.methods)

        logger.debug(f"Parsed class: {qualified_name} (flattened: {flattened_name}, {len(cpp_class.methods)} methods)")
        return cpp_class

    def _process_class_template(
        self,
        cursor,
        namespace: str,
        out_classes: Optional[Dict[str, CppClass]] = None,
        parent_class_name: str = "",
    ) -> Optional[CppClass]:
        """Process a class template declaration."""
        class_name = cursor.spelling
        if not class_name:
            return None
        if not cursor.is_definition():
            return None

        # For nested templates, create a flattened name
        if parent_class_name:
            flattened_name = f"{parent_class_name}_{class_name}"
            qualified_name = f"{namespace}::{parent_class_name}::{class_name}" if namespace else f"{parent_class_name}::{class_name}"
        else:
            flattened_name = class_name
            qualified_name = f"{namespace}::{class_name}" if namespace else class_name

        cpp_class = CppClass(
            name=flattened_name,
            qualified_name=qualified_name,
            namespace=namespace,
            header_file=self._normalize_filename(str(cursor.location.file)),
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
                # Get base classes (both original and canonical)
                for subchild in child.get_children():
                    if subchild.kind == CursorKind.CXX_BASE_SPECIFIER:
                        base_type = subchild.type
                        cpp_class.base_classes.append(base_type.spelling)
                        cpp_class.canonical_base_classes.append(base_type.get_canonical().spelling)

                self._process_class_members(child, cpp_class, out_classes, flattened_name)
                break

        cpp_class.is_abstract = any(m.is_pure_virtual for m in cpp_class.methods)
        return cpp_class

    def _process_class_members(
        self,
        cursor,
        cpp_class: CppClass,
        out_classes: Optional[Dict[str, CppClass]] = None,
        parent_class_name: str = "",
    ):
        """Process class members (methods, constructors, nested classes).

        Parameters
        ----------
        cursor : clang.cindex.Cursor
            The class cursor.
        cpp_class : CppClass
            The CppClass object to populate.
        out_classes : dict, optional
            Dictionary to collect nested classes into.
        parent_class_name : str
            The flattened name of the current class (for nested class name generation).
        """
        # Default access is public for structs, private for classes
        if cursor.kind == CursorKind.STRUCT_DECL:
            current_access = "public"
        else:
            current_access = "private"

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

            # Handle nested classes and structs
            elif child.kind in (CursorKind.CLASS_DECL, CursorKind.STRUCT_DECL):
                if out_classes is not None and child.is_definition():
                    nested_class = self._process_class(
                        child,
                        cpp_class.namespace,
                        out_classes=out_classes,
                        parent_class_name=parent_class_name,
                    )
                    if nested_class:
                        out_classes[nested_class.name] = nested_class
                        logger.debug(f"Found nested class: {nested_class.name} in {parent_class_name}")

    def _process_method(self, cursor, is_constructor: bool = False) -> CppMethod:
        """Process a method declaration."""
        result_type = cursor.result_type
        return_type = result_type.spelling if not is_constructor else ""
        # Get canonical type (resolves typedefs like Peak1D::IntensityType -> float)
        canonical_return_type = result_type.get_canonical().spelling if not is_constructor else ""

        method = CppMethod(
            name=cursor.spelling,
            return_type=return_type,
            canonical_return_type=canonical_return_type,
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
        # Get canonical type (resolves typedefs)
        canonical_type = param_type.get_canonical().spelling

        return CppParameter(
            name=cursor.spelling or f"arg{cursor.hash}",
            type_str=type_str,
            canonical_type=canonical_type,
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


def _create_fallback_merged_class(class_name: str, pxd_class: "ClassDecl") -> Optional[MergedClass]:
    """
    Create a MergedClass from .pxd info alone (fallback for when libclang fails).

    This creates a minimal binding using only .pxd declarations. Methods will
    need to be added via SPECIAL_METHODS in the emitter.
    """
    # Create a minimal CppClass from .pxd info
    cpp_class = CppClass(
        name=class_name,
        qualified_name=f"{pxd_class.namespace}::{class_name}",
        namespace=pxd_class.namespace,
        header_file=pxd_class.header_file.strip("<>") if pxd_class.header_file else "",
        base_classes=pxd_class.base_classes,
        is_abstract=pxd_class.is_abstract,
        is_template=bool(pxd_class.template_args),
        template_params=pxd_class.template_args,
        doc=pxd_class.doc,
    )

    # Convert .pxd methods to CppMethod and MergedMethod
    # Note: Type info from .pxd is less accurate, but we use it as best effort
    merged_methods = []
    for pxd_m in pxd_class.methods:
        # Skip ignored methods
        if pxd_m.wrap_ignore:
            continue
        # Create basic CppMethod from .pxd
        cpp_method = CppMethod(
            name=pxd_m.name,
            return_type=pxd_m.return_type or "void",
            parameters=[
                CppParameter(
                    name=p.name or f"arg{i}",
                    type_str=getattr(p, 'type_str', getattr(p, 'type', 'void')),
                    default_value=getattr(p, 'default_value', getattr(p, 'default', None))
                )
                for i, p in enumerate(pxd_m.parameters)
            ],
            is_const=pxd_m.is_const,
            is_static=pxd_m.is_static,
            access="public",
        )
        merged_method = MergedMethod(
            cpp_method=cpp_method,
            wrap_as=pxd_m.wrap_as,
            wrap_ignore=pxd_m.wrap_ignore,
            doc=pxd_m.doc or "",
        )
        merged_methods.append(merged_method)

    # Create constructors from pxd
    constructors = []
    for pxd_ctor in pxd_class.constructors:
        cpp_ctor = CppMethod(
            name=pxd_ctor.name,
            return_type="",
            parameters=[
                CppParameter(
                    name=p.name or f"arg{i}",
                    type_str=p.type_str,
                    canonical_type=p.type_str,
                    is_const="const" in p.type_str,
                    is_reference="&" in p.type_str,
                    is_pointer="*" in p.type_str,
                    default_value=p.default_value,
                )
                for i, p in enumerate(pxd_ctor.parameters)
            ],
            access="public",
            is_constructor=True,
        )
        constructors.append(cpp_ctor)

    # Add default constructor if no constructors declared but the class looks instantiable
    if not constructors and not pxd_class.is_abstract:
        constructors.append(CppMethod(
            name=class_name,
            return_type="",
            parameters=[],
            access="public",
            is_constructor=True,
        ))

    return MergedClass(
        cpp_class=cpp_class,
        methods=merged_methods,
        constructors=constructors,
        wrap_hash=getattr(pxd_class, 'wrap_hash', False),
        wrap_iter=getattr(pxd_class, 'wrap_iter', None),
        wrap_manual_memory=getattr(pxd_class, 'wrap_manual_memory', False),
        doc=getattr(pxd_class, 'doc', ''),
    )


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

    # Classes that should fall back to .pxd info when libclang can't parse them
    # (complex template inheritance, missing includes, etc.)
    # NOTE: Template aliases (DRange1 -> DRange<1>) need special handling and
    # cannot use simple .pxd fallback - they need typedef/using mappings.
    FALLBACK_TO_PXD = {
        "Mobilogram",  # Complex RangeManagerContainer template inheritance
    }

    for class_name, pxd_class in pxd_classes.items():
        if class_name not in cpp_classes:
            # Check if this class should use .pxd fallback
            if class_name in FALLBACK_TO_PXD:
                logger.info(f"Class {class_name} using .pxd fallback (not found in C++ headers)")
                # Use pxd_to_merged for single class
                fallback = _create_fallback_merged_class(class_name, pxd_class)
                if fallback:
                    merged[class_name] = fallback
                continue
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

        # Filter constructors from C++ info
        constructors = [c for c in cpp_class.constructors if c.access == "public"]

        # If no constructors from libclang, fall back to .pxd constructors
        if not constructors and pxd_class.constructors:
            for pxd_ctor in pxd_class.constructors:
                # Skip constructors with wrap_ignore=True
                if getattr(pxd_ctor, 'wrap_ignore', False):
                    continue
                cpp_ctor = CppMethod(
                    name=pxd_ctor.name,
                    return_type="",
                    parameters=[
                        CppParameter(
                            name=p.name,
                            type_str=p.type_str,
                            canonical_type=p.type_str,
                            is_const="const" in p.type_str,
                            is_reference="&" in p.type_str,
                            is_pointer="*" in p.type_str,
                            default_value=p.default_value,
                        )
                        for p in pxd_ctor.parameters
                    ],
                    access="public",
                    is_constructor=True,
                )
                constructors.append(cpp_ctor)

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


def pxd_to_merged(pxd_classes: Dict[str, "ClassDecl"]) -> Dict[str, MergedClass]:
    """
    Convert .pxd class declarations to MergedClass format without C++ parsing.

    This allows using NanobindEmitterV2 even without libclang, using type
    information from the .pxd files directly.

    Parameters
    ----------
    pxd_classes : dict
        Classes declared in .pxd files (from PxdParser).

    Returns
    -------
    dict
        MergedClass objects suitable for NanobindEmitterV2.
    """
    merged = {}

    for class_name, pxd_class in pxd_classes.items():
        # Create a minimal CppClass from .pxd info
        cpp_class = CppClass(
            name=class_name,
            qualified_name=f"{pxd_class.namespace}::{class_name}",
            namespace=pxd_class.namespace,
            header_file=pxd_class.header_file.strip("<>") if pxd_class.header_file else "",
            base_classes=pxd_class.base_classes,
            is_abstract=pxd_class.is_abstract,
            is_template=bool(pxd_class.template_args),
            template_params=pxd_class.template_args,
            doc=pxd_class.doc,
        )

        # Convert methods
        merged_methods = []
        for m in pxd_class.methods:
            # Create CppMethod from pxd Method
            cpp_method = CppMethod(
                name=m.name,
                return_type=m.return_type,
                parameters=[
                    CppParameter(
                        name=p.name,
                        type_str=p.type_str,
                        default_value=p.default_value,
                    )
                    for p in m.parameters
                ],
                is_const=m.is_const,
                is_static=m.is_static,
                is_constructor=m.is_constructor,
                is_destructor=m.is_destructor,
                access="public",
            )
            merged_methods.append(MergedMethod(
                cpp_method=cpp_method,
                wrap_as=m.wrap_as,
                wrap_ignore=m.wrap_ignore,
                doc=m.doc,
            ))

        # Convert constructors
        constructors = []
        for c in pxd_class.constructors:
            # Skip constructors with wrap_ignore=True
            if getattr(c, 'wrap_ignore', False):
                continue
            cpp_ctor = CppMethod(
                name=class_name,
                return_type="",
                parameters=[
                    CppParameter(
                        name=p.name,
                        type_str=p.type_str,
                        default_value=p.default_value,
                    )
                    for p in c.parameters
                ],
                is_constructor=True,
                access="public",
            )
            constructors.append(cpp_ctor)

        merged_class = MergedClass(
            cpp_class=cpp_class,
            methods=merged_methods,
            constructors=constructors,
            wrap_hash=bool(pxd_class.wrap_hash),
            wrap_iter=pxd_class.wrap_iter,
            wrap_manual_memory=pxd_class.wrap_manual_memory,
            doc=pxd_class.doc,
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
