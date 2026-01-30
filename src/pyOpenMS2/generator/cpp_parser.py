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
    is_deleted: bool = False  # = delete
    is_defaulted: bool = False  # = default
    access: str = "public"  # public, protected, private
    uses_incomplete_type: bool = False  # Uses forward-declared/incomplete types
    incomplete_types: List[str] = field(default_factory=list)  # Names of incomplete types used

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
    # Overload tracking: method names that have multiple signatures
    overloaded_methods: Set[str] = field(default_factory=set)
    # Methods with const/non-const overloads (common binding issue)
    const_overloaded_methods: Set[str] = field(default_factory=set)
    # Constructor status
    has_deleted_default_constructor: bool = False
    has_deleted_copy_constructor: bool = False
    has_private_constructor: bool = False  # All constructors are private/protected

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
            'overloaded_methods': list(self.overloaded_methods),
            'const_overloaded_methods': list(self.const_overloaded_methods),
            'has_deleted_default_constructor': self.has_deleted_default_constructor,
            'has_deleted_copy_constructor': self.has_deleted_copy_constructor,
            'has_private_constructor': self.has_private_constructor,
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
            overloaded_methods=set(d.get('overloaded_methods', [])),
            const_overloaded_methods=set(d.get('const_overloaded_methods', [])),
            has_deleted_default_constructor=d.get('has_deleted_default_constructor', False),
            has_deleted_copy_constructor=d.get('has_deleted_copy_constructor', False),
            has_private_constructor=d.get('has_private_constructor', False),
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
        self._scoped_enums: Set[str] = set()  # Qualified names of C++11 scoped enums

        # Create cache directory if specified
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"libclang cache directory: {self.cache_dir}")

    @staticmethod
    def _normalize_filename(filename: str) -> str:
        """Normalize a filename for stable comparisons across platforms."""
        return os.path.normpath(os.path.abspath(filename))

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

    def get_scoped_enums(self) -> Set[str]:
        """Get the set of qualified names of scoped enums (C++11 enum class)."""
        return self._scoped_enums

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

            elif child.kind == CursorKind.ENUM_DECL:
                # Detect scoped enums (C++11 enum class)
                if child.is_scoped_enum():
                    enum_name = child.spelling
                    qualified_name = f"{namespace}::{enum_name}" if namespace else enum_name
                    self._scoped_enums.add(qualified_name)

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

        # Collect inherited public methods from base classes
        self._collect_inherited_methods(cursor, cpp_class)

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

        # Track all constructors (including non-public) for deleted/private detection
        all_constructors: List[CppMethod] = []
        has_any_public_constructor = False

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

            if child.kind == CursorKind.CONSTRUCTOR:
                ctor = self._process_method(child, is_constructor=True)
                ctor.access = current_access
                all_constructors.append(ctor)
                if current_access == "public":
                    has_any_public_constructor = True
                    if not ctor.is_deleted:
                        cpp_class.constructors.append(ctor)

            elif child.kind == CursorKind.DESTRUCTOR:
                # Check for pure virtual destructor (makes class abstract)
                if child.is_pure_virtual_method():
                    cpp_class.is_abstract = True

            elif child.kind == CursorKind.CXX_METHOD:
                method = self._process_method(child)
                method.access = current_access
                # Check for pure virtual (makes class abstract)
                if method.is_pure_virtual:
                    cpp_class.is_abstract = True
                # Only add public methods to binding list
                if current_access == "public":
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

        # Analyze constructor status
        for ctor in all_constructors:
            # Check for deleted default constructor (no parameters)
            if len(ctor.parameters) == 0 and ctor.is_deleted:
                cpp_class.has_deleted_default_constructor = True
            # Check for deleted copy constructor (single const ref parameter of same type)
            if len(ctor.parameters) == 1 and ctor.is_deleted:
                param_type = ctor.parameters[0].type_str
                if "const" in param_type and "&" in param_type:
                    cpp_class.has_deleted_copy_constructor = True

        # Check if all constructors are non-public
        if all_constructors and not has_any_public_constructor:
            cpp_class.has_private_constructor = True

        # Detect method overloads
        self._detect_overloads(cpp_class)

    def _collect_inherited_methods(self, cursor, cpp_class: CppClass):
        """Collect public methods from base classes (recursively).

        This ensures that classes with skipped base class specification
        (e.g., MSSpectrum with private std::vector base) still expose
        inherited methods like setNativeID from SpectrumSettings.
        """
        existing_methods = {m.name for m in cpp_class.methods}
        visited = set()

        def _walk_bases(cls_cursor):
            for child in cls_cursor.get_children():
                if child.kind != CursorKind.CXX_BASE_SPECIFIER:
                    continue
                # Resolve the base class definition
                base_type = child.type
                base_decl = base_type.get_declaration()
                if not base_decl or not base_decl.location or not base_decl.location.file:
                    continue
                base_id = (base_decl.spelling, str(base_decl.location.file))
                if base_id in visited:
                    continue
                visited.add(base_id)

                # Collect public methods from this base
                for member in base_decl.get_children():
                    if member.kind == CursorKind.CXX_METHOD:
                        if member.access_specifier == AccessSpecifier.PUBLIC:
                            if member.spelling not in existing_methods:
                                try:
                                    method = self._process_method(member)
                                    method.access = "public"
                                    cpp_class.methods.append(method)
                                    existing_methods.add(method.name)
                                except Exception:
                                    pass  # Skip methods that fail to parse

                # Recurse into base's bases
                _walk_bases(base_decl)

        _walk_bases(cursor)

    def _process_method(self, cursor, is_constructor: bool = False) -> CppMethod:
        """Process a method declaration."""
        result_type = cursor.result_type
        return_type = result_type.spelling if not is_constructor else ""
        # Get canonical type (resolves typedefs like Peak1D::IntensityType -> float)
        canonical_return_type = result_type.get_canonical().spelling if not is_constructor else ""

        # Detect deleted and defaulted methods
        # libclang doesn't have a direct is_deleted() method, so we check tokens
        is_deleted = False
        is_defaulted = False
        try:
            tokens = list(cursor.get_tokens())
            token_spellings = [t.spelling for t in tokens]
            # Look for "= delete" or "= default" pattern
            for i, spelling in enumerate(token_spellings):
                if spelling == "=" and i + 1 < len(token_spellings):
                    next_token = token_spellings[i + 1]
                    if next_token == "delete":
                        is_deleted = True
                    elif next_token == "default":
                        is_defaulted = True
        except Exception:
            pass  # Token access may fail in some cases

        # Detect incomplete (forward-declared) types
        incomplete_types = []

        # Check return type for incomplete types
        if not is_constructor:
            incomplete = self._check_type_incomplete(result_type)
            if incomplete:
                incomplete_types.append(incomplete)

        method = CppMethod(
            name=cursor.spelling,
            return_type=return_type,
            canonical_return_type=canonical_return_type,
            is_constructor=is_constructor,
            is_const=cursor.is_const_method(),
            is_static=cursor.is_static_method(),
            is_virtual=cursor.is_virtual_method(),
            is_pure_virtual=cursor.is_pure_virtual_method(),
            is_deleted=is_deleted,
            is_defaulted=is_defaulted,
        )

        # Process parameters and check for incomplete types
        for child in cursor.get_children():
            if child.kind == CursorKind.PARM_DECL:
                param = self._process_parameter(child)
                method.parameters.append(param)
                # Check parameter type for incomplete types
                incomplete = self._check_type_incomplete(child.type)
                if incomplete:
                    incomplete_types.append(incomplete)

        # Set incomplete type info
        if incomplete_types:
            method.uses_incomplete_type = True
            method.incomplete_types = incomplete_types

        return method

    def _check_type_incomplete(self, clang_type) -> Optional[str]:
        """Check if a type is incomplete (forward-declared).

        Returns the type name if incomplete, None otherwise.
        Template specializations of known wrapped types (e.g., DBoundingBox<2>)
        are not considered incomplete.
        """
        # Get the underlying type for pointers/references
        pointee = clang_type
        if clang_type.kind in (TypeKind.POINTER, TypeKind.LVALUEREFERENCE, TypeKind.RVALUEREFERENCE):
            pointee = clang_type.get_pointee()

        # Get the declaration of the type
        decl = pointee.get_declaration()
        if decl and decl.kind in (CursorKind.CLASS_DECL, CursorKind.STRUCT_DECL):
            # Check if it's a forward declaration (not a definition)
            if not decl.is_definition():
                # Template specializations of known wrapped types appear
                # "incomplete" to libclang but are fully usable since the
                # template itself is defined.  Only whitelist types we know
                # have complete definitions and nanobind type casters.
                if pointee.get_num_template_arguments() > 0:
                    name = decl.spelling
                    if name in _KNOWN_COMPLETE_TEMPLATES:
                        return None
                return decl.spelling

        return None

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

    def _detect_overloads(self, cpp_class: CppClass):
        """Detect overloaded methods and const/non-const overloads.

        Populates cpp_class.overloaded_methods and cpp_class.const_overloaded_methods.
        """
        from collections import defaultdict

        # Group methods by name
        methods_by_name: Dict[str, List[CppMethod]] = defaultdict(list)
        for method in cpp_class.methods:
            methods_by_name[method.name].append(method)

        for name, methods in methods_by_name.items():
            if len(methods) > 1:
                # This method name has multiple overloads
                cpp_class.overloaded_methods.add(name)

                # Check for const/non-const overload pattern
                const_methods = [m for m in methods if m.is_const]
                nonconst_methods = [m for m in methods if not m.is_const]

                if const_methods and nonconst_methods:
                    # Has both const and non-const versions
                    cpp_class.const_overloaded_methods.add(name)


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
    enums: List[Any] = field(default_factory=list)  # EnumDecl from pxd_parser
    template_instances: Dict[str, List[str]] = field(default_factory=dict)  # instance_name -> [template_args]
    pxd_template_params: List[str] = field(default_factory=list)  # Template param names from .pxd
    member_variables: List[Any] = field(default_factory=list)  # MemberVariable from pxd_parser

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

    @property
    def overloaded_methods(self) -> Set[str]:
        """Method names that have multiple overloads."""
        return self.cpp_class.overloaded_methods

    @property
    def const_overloaded_methods(self) -> Set[str]:
        """Method names that have const/non-const overloads."""
        return self.cpp_class.const_overloaded_methods

    @property
    def has_deleted_default_constructor(self) -> bool:
        return self.cpp_class.has_deleted_default_constructor

    @property
    def has_deleted_copy_constructor(self) -> bool:
        return self.cpp_class.has_deleted_copy_constructor

    @property
    def has_private_constructor(self) -> bool:
        return self.cpp_class.has_private_constructor


# Template classes that libclang marks as "incomplete" (forward-declared)
# for implicit specializations, but are fully defined and have nanobind
# type casters or bindings.  Methods using these as return/parameter types
# should NOT be skipped.
_KNOWN_COMPLETE_TEMPLATES = {
    "DBoundingBox",   # DBoundingBox<2> wrapped as DBoundingBox2
    "DPosition",      # DPosition<1>, DPosition<2> have type casters
    "DRange",         # DRange<1>, DRange<2>
    "Matrix",         # Matrix<double> wrapped as MatrixDouble
}

# Mapping from Cython/pxd type names to C++ types for fallback classes
_PXD_TO_CPP_TYPE = {
    "libcpp_utf8_string": "std::string",
    "libcpp_string": "std::string",
    "libcpp_utf8_output_string": "std::string",
    "bool": "bool",
    "int": "int",
    "long": "long",
    "float": "float",
    "double": "double",
    "size_t": "size_t",
    "UInt": "OpenMS::UInt",
    "UInt32": "OpenMS::UInt32",
    "UInt64": "OpenMS::UInt64",
    "Int": "OpenMS::Int",
    "Size": "size_t",
    "SignedSize": "ptrdiff_t",
}


def _convert_pxd_type(type_str: str) -> str:
    """Convert a pxd type string to a C++ type string for fallback classes."""
    if not type_str:
        return type_str

    result = type_str.strip()

    # Handle const and reference modifiers
    is_const = result.startswith("const ")
    if is_const:
        result = result[6:]
    has_ref = result.endswith("&")
    if has_ref:
        result = result[:-1].strip()
    has_ptr = result.endswith("*")
    if has_ptr:
        result = result[:-1].strip()

    # Direct type mapping
    if result in _PXD_TO_CPP_TYPE:
        result = _PXD_TO_CPP_TYPE[result]
    # Handle libcpp_vector[T] -> std::vector<T>
    elif result.startswith("libcpp_vector["):
        inner = result[len("libcpp_vector["):-1]
        inner = _convert_pxd_type(inner)
        result = f"std::vector<{inner}>"
    elif result.startswith("libcpp_pair["):
        inner = result[len("libcpp_pair["):-1]
        parts = [_convert_pxd_type(p.strip()) for p in inner.split(",", 1)]
        result = f"std::pair<{', '.join(parts)}>"
    elif result.startswith("libcpp_set["):
        inner = result[len("libcpp_set["):-1]
        inner = _convert_pxd_type(inner)
        result = f"std::set<{inner}>"
    elif result.startswith("libcpp_map["):
        inner = result[len("libcpp_map["):-1]
        parts = [_convert_pxd_type(p.strip()) for p in inner.split(",", 1)]
        result = f"std::map<{', '.join(parts)}>"
    elif result.startswith("shared_ptr["):
        inner = result[len("shared_ptr["):-1]
        inner = _convert_pxd_type(inner)
        result = f"std::shared_ptr<{inner}>"
    # Note: Don't auto-qualify unrecognized types here — the emitter's
    # _normalize_type handles namespace resolution with full context.

    # Reconstruct with modifiers
    if is_const:
        result = f"const {result}"
    if has_ref:
        result = f"{result} &"
    if has_ptr:
        result = f"{result} *"

    return result


def _create_fallback_merged_class(class_name: str, pxd_class: "ClassDecl") -> Optional[MergedClass]:
    """
    Create a MergedClass from .pxd info alone (fallback for when libclang fails).

    This creates a minimal binding using only .pxd declarations. Methods will
    need to be added via SPECIAL_METHODS in the emitter.
    """
    # Create a minimal CppClass from .pxd info
    # Use the cpp_qualified_name from pxd if available (e.g., "OpenMS::PeakIntegrator::PeakArea")
    # otherwise fall back to namespace::class_name
    if pxd_class.cpp_qualified_name:
        qualified_name = pxd_class.cpp_qualified_name
    else:
        qualified_name = f"{pxd_class.namespace}::{class_name}"
    cpp_class = CppClass(
        name=class_name,
        qualified_name=qualified_name,
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
        # Create basic CppMethod from .pxd (converting pxd types to C++ types)
        cpp_method = CppMethod(
            name=pxd_m.name,
            return_type=_convert_pxd_type(pxd_m.return_type or "void"),
            parameters=[
                CppParameter(
                    name=p.name or f"arg{i}",
                    type_str=_convert_pxd_type(getattr(p, 'type_str', getattr(p, 'type', 'void'))),
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
                    type_str=_convert_pxd_type(p.type_str),
                    canonical_type=_convert_pxd_type(p.type_str),
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
        enums=getattr(pxd_class, 'enums', []),
        template_instances=getattr(pxd_class, 'template_instances', {}),
        pxd_template_params=getattr(pxd_class, 'template_args', []),
        member_variables=[
            type('MemberVariable', (), {'name': mv.name, 'type_str': _convert_pxd_type(mv.type_str)})()
            for mv in getattr(pxd_class, 'member_variables', [])
        ],
    )


def merge_with_pxd(
    cpp_classes: Dict[str, CppClass],
    pxd_classes: Dict[str, "ClassDecl"],  # From pxd_parser
    scoped_enums: Optional[Set[str]] = None,  # Qualified names of C++11 enum class
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
    scoped_enums : set, optional
        Qualified names of C++11 scoped enums (enum class) from libclang.

    Returns
    -------
    dict
        Merged class information with accurate C++ types.
    """
    scoped_enums = scoped_enums or set()
    merged = {}

    # Classes that should fall back to .pxd info when libclang can't parse them
    # (complex template inheritance, missing includes, etc.)
    # NOTE: Template aliases (DRange1 -> DRange<1>) need special handling and
    # cannot use simple .pxd fallback - they need typedef/using mappings.
    FALLBACK_TO_PXD = {
        "Mobilogram",  # Complex RangeManagerContainer template inheritance
        # Nested classes (pxd name != C++ flattened name)
        "ColumnHeader",                  # ConsensusMap::ColumnHeader
        "ExtractionCoordinates",         # ChromatogramExtractorAlgorithm::ExtractionCoordinates
        "MRMFQC_ComponentQCs",           # MRMFeatureQC::ComponentQCs
        "MRMFQC_ComponentGroupQCs",      # MRMFeatureQC::ComponentGroupQCs
        "MRMFQC_ComponentGroupPairQCs",  # MRMFeatureQC::ComponentGroupPairQCs
        "SelectorParameters",            # MRMFeatureSelector::SelectorParameters
        "PI_PeakArea",                   # PeakIntegrator::PeakArea
        "PI_PeakBackground",             # PeakIntegrator::PeakBackground
        "PI_PeakShapeMetrics",           # PeakIntegrator::PeakShapeMetrics
        # Template classes (may not match in libclang by name)
        "Matrix",                        # OpenMS::Matrix<T> (template)
        "BilinearInterpolation",         # OpenMS::Math::BilinearInterpolation<K,V> (template)
        "MRMTransitionGroup",            # OpenMS::MRMTransitionGroup<C,T> (template)
        # Template specialization with cpp name alias
        "DBoundingBox2",                 # OpenMS::DBoundingBox<2>
        # Non-OpenMS namespace classes (OpenSwath, FLASHHelperClasses)
        "LightTransition",               # OpenSwath::LightTransition
        "LightModification",             # OpenSwath::LightModification
        "LightCompound",                 # OpenSwath::LightCompound
        "LightProtein",                  # OpenSwath::LightProtein
        "LightTargetedExperiment",       # OpenSwath::LightTargetedExperiment
        "LogMzPeak",                     # FLASHHelperClasses::LogMzPeak
        "PrecalAveragine",               # FLASHHelperClasses::PrecalculatedAveragine
        "MassFeature_FDHS",              # FLASHHelperClasses::MassFeature
        "IsobaricQuantities",            # FLASHHelperClasses::IsobaricQuantities
    }

    for class_name, pxd_class in pxd_classes.items():
        if class_name not in cpp_classes:
            if class_name in FALLBACK_TO_PXD:
                logger.info(f"Class {class_name} using .pxd fallback (not found in C++ headers)")
                fallback = _create_fallback_merged_class(class_name, pxd_class)
                if fallback:
                    merged[class_name] = fallback
            else:
                logger.debug(f"Class {class_name} in .pxd but not found in C++ headers")
            continue

        cpp_class = cpp_classes[class_name]

        # For FALLBACK_TO_PXD classes that libclang found but with no methods
        # (e.g., template classes where libclang can't enumerate methods),
        # use the pxd fallback which has complete method declarations.
        if class_name in FALLBACK_TO_PXD and not cpp_class.methods:
            logger.info(f"Class {class_name} using .pxd fallback (libclang found 0 methods)")
            fallback = _create_fallback_merged_class(class_name, pxd_class)
            if fallback:
                merged[class_name] = fallback
            continue

        # Build method allowlist from .pxd with directives
        # Include methods from base class pxds (wrap-inherits)
        pxd_methods_by_name = {}
        for m in pxd_class.methods:
            if m.name not in pxd_methods_by_name:
                pxd_methods_by_name[m.name] = []
            pxd_methods_by_name[m.name].append(m)

        # Add inherited methods from base class pxds (recursively)
        def _add_base_pxd_methods(base_names):
            for base_name in base_names:
                if base_name in pxd_classes:
                    for m in pxd_classes[base_name].methods:
                        if m.name not in pxd_methods_by_name:
                            pxd_methods_by_name[m.name] = []
                            pxd_methods_by_name[m.name].append(m)
                    # Recurse into base's bases
                    _add_base_pxd_methods(getattr(pxd_classes[base_name], 'base_classes', []))
        _add_base_pxd_methods(getattr(pxd_class, 'base_classes', []))

        if class_name == "MSSpectrum":
            logger.info(f"MSSpectrum pxd base_classes: {getattr(pxd_class, 'base_classes', [])}")
            logger.info(f"MSSpectrum pxd allowlist has setNativeID: {'setNativeID' in pxd_methods_by_name}")
            logger.info(f"MSSpectrum pxd allowlist has setPrecursors: {'setPrecursors' in pxd_methods_by_name}")
            logger.info(f"MSSpectrum cpp methods count: {len(cpp_class.methods)}")
            logger.info(f"MSSpectrum cpp has setNativeID: {any(m.name == 'setNativeID' for m in cpp_class.methods)}")

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

        # Update enum is_scoped based on libclang detection
        pxd_enums = getattr(pxd_class, 'enums', [])
        for enum in pxd_enums:
            # Check if this enum is a scoped enum (C++11 enum class)
            # Use cpp_name if available, otherwise construct from namespace
            qualified = enum.cpp_name if enum.cpp_name else f"{enum.namespace}::{enum.name}"
            if qualified in scoped_enums:
                enum.is_scoped = True

        merged_class = MergedClass(
            cpp_class=cpp_class,
            methods=merged_methods,
            constructors=constructors,
            wrap_hash=getattr(pxd_class, 'wrap_hash', False),
            wrap_iter=getattr(pxd_class, 'wrap_iter', None),
            wrap_manual_memory=getattr(pxd_class, 'wrap_manual_memory', False),
            doc=getattr(pxd_class, 'doc', ''),
            enums=pxd_enums,
            template_instances=getattr(pxd_class, 'template_instances', {}),
            pxd_template_params=getattr(pxd_class, 'template_args', []),
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
            enums=getattr(pxd_class, 'enums', []),
            template_instances=getattr(pxd_class, 'template_instances', {}),
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
