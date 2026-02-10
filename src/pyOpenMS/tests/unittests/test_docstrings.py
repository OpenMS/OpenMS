"""
Tests for docstring presence, type, and RST validity across all wrapped classes.

Validates that all pyopenms classes and their public methods have properly
formatted docstrings that are compatible with Sphinx documentation generation.
"""

import inspect
import re
import warnings
import pytest

try:
    from docutils.parsers.rst import Parser
    from docutils.utils import Reporter, new_document
    from docutils.frontend import OptionParser
    HAS_DOCUTILS = True
except ImportError:
    HAS_DOCUTILS = False

import pyopenms


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_all_classes():
    """Return all public classes exported by pyopenms."""
    return sorted(
        name for name in dir(pyopenms)
        if not name.startswith("_")
        and inspect.isclass(getattr(pyopenms, name, None))
    )


def _get_public_methods(cls):
    """Return public method names for a class (including __init__)."""
    import enum
    skip = set()
    if issubclass(cls, enum.IntEnum):
        skip = set(dir(int)) - set(dir(enum.Enum))
    methods = []
    for name in sorted(dir(cls)):
        if name in skip:
            continue
        if name.startswith("_") and name != "__init__":
            continue
        obj = getattr(cls, name, None)
        if callable(obj):
            methods.append(name)
    return methods


def _check_rst_code_blocks(docstring):
    """Check for improperly formatted RST code blocks."""
    errors = []
    lines = docstring.split("\n")
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith(".. code-block::"):
            lang = stripped[len(".. code-block::"):].strip()
            if not lang:
                errors.append(f"Line {i+1}: code-block missing language specifier")
        if stripped == "::":
            # Literal block — next non-empty line must be indented
            for j in range(i + 1, len(lines)):
                next_line = lines[j]
                if next_line.strip():
                    if not next_line.startswith(" ") and not next_line.startswith("\t"):
                        errors.append(
                            f"Line {j+1}: literal block content not indented"
                        )
                    break
    return errors


def _check_common_issues(docstring):
    """Check for common docstring problems."""
    errors = []
    # Unbalanced backticks (single or double)
    single = docstring.count("`") - 2 * docstring.count("``")
    if single % 2 != 0:
        errors.append("Unbalanced single backticks")
    # Tabs (RST uses spaces)
    if "\t" in docstring:
        errors.append("Contains tab characters (use spaces for RST)")
    return errors


def _get_docutils_errors(docstring):
    """Parse docstring with docutils and return (level, message) pairs."""
    if not HAS_DOCUTILS:
        return []
    parser = Parser()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        option_parser = OptionParser(components=(Parser,))
        settings = option_parser.get_default_values()
    settings.report_level = Reporter.WARNING_LEVEL
    settings.halt_level = Reporter.SEVERE_LEVEL + 1
    doc = new_document("<docstring>", settings)
    parser.parse(docstring, doc)
    return [
        (node["level"], node.astext())
        for node in doc.findall()
        if node.tagname == "system_message"
    ]


def _is_acceptable_error(error, docstring):
    """Filter RST errors that are expected in generated docstrings."""
    err = error.lower()
    if "hyperlink" in err and "here" in docstring.lower():
        return True
    if "field list" in err:
        return True
    if 'unknown interpreted text role "py:' in err:
        return True
    if 'unknown directive type "py:' in err:
        return True
    if 'unknown directive type "deprecated"' in err:
        return True
    # NumPy-style docstring section headers (Parameters, Returns, Raises, etc.)
    # are valid for Sphinx napoleon but not for raw docutils RST parsing
    if "unexpected section title" in err:
        return True
    # nanobind auto-generated signatures may contain C++ pointer syntax (*)
    if "inline emphasis start-string without end-string" in err:
        return True
    if "inline strong start-string without end-string" in err:
        return True
    return False


def validate_docstring(docstring, name=""):
    """Validate a docstring for RST compatibility.

    Returns (is_valid, errors_list).
    """
    if docstring is None:
        return True, []
    if not isinstance(docstring, str):
        return False, [f"{name}: docstring is {type(docstring).__name__}, not str"]
    errors = []
    errors.extend(_check_rst_code_blocks(docstring))
    errors.extend(_check_common_issues(docstring))
    if HAS_DOCUTILS:
        for level, msg in _get_docutils_errors(docstring):
            if level >= 2:
                sev = {2: "WARNING", 3: "ERROR", 4: "SEVERE"}.get(level, "?")
                errors.append(f"RST {sev}: {msg}")
    return len(errors) == 0, errors


# ---------------------------------------------------------------------------
# Data: collect all classes once at module level for parametrize
# ---------------------------------------------------------------------------

ALL_CLASSES = _get_all_classes()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestAllClassDocstrings:
    """Validate docstrings across every wrapped class."""

    @pytest.mark.parametrize("class_name", ALL_CLASSES)
    def test_class_docstring_is_str(self, class_name):
        """Every class docstring (if present) must be str, not bytes."""
        cls = getattr(pyopenms, class_name)
        doc = cls.__doc__
        if doc is not None:
            assert isinstance(doc, str), (
                f"{class_name}.__doc__ is {type(doc).__name__}, expected str"
            )

    @pytest.mark.parametrize("class_name", ALL_CLASSES)
    def test_class_docstring_valid_rst(self, class_name):
        """Every class docstring must be valid RST."""
        cls = getattr(pyopenms, class_name)
        doc = cls.__doc__
        if not doc:
            pytest.skip("no docstring")
        ok, errors = validate_docstring(doc, class_name)
        if not ok:
            real = [e for e in errors if not _is_acceptable_error(e, doc)]
            assert not real, f"{class_name} docstring RST errors:\n" + "\n".join(real)

    @pytest.mark.parametrize("class_name", ALL_CLASSES)
    def test_method_docstrings_are_str(self, class_name):
        """All method docstrings must be str, not bytes."""
        cls = getattr(pyopenms, class_name)
        binary = []
        for mname in _get_public_methods(cls):
            doc = getattr(getattr(cls, mname, None), "__doc__", None)
            if doc is not None and isinstance(doc, bytes):
                binary.append(mname)
        assert not binary, (
            f"{class_name}: methods with binary docstrings: {binary}"
        )

    @pytest.mark.parametrize("class_name", ALL_CLASSES)
    def test_method_docstrings_valid_rst(self, class_name):
        """All method docstrings must be valid RST."""
        cls = getattr(pyopenms, class_name)
        all_errors = []
        for mname in _get_public_methods(cls):
            doc = getattr(getattr(cls, mname, None), "__doc__", None)
            if not doc:
                continue
            ok, errors = validate_docstring(doc, f"{class_name}.{mname}")
            if not ok:
                real = [e for e in errors if not _is_acceptable_error(e, doc)]
                if real:
                    all_errors.extend(real)
        assert not all_errors, (
            f"{class_name} method RST errors:\n" + "\n".join(all_errors[:20])
        )


class TestDocstringCoverage:
    """Summary statistics on docstring coverage (informational)."""

    def test_class_docstring_coverage(self):
        """Report how many classes have docstrings."""
        total = len(ALL_CLASSES)
        with_doc = sum(
            1 for name in ALL_CLASSES
            if getattr(pyopenms, name).__doc__
        )
        pct = 100.0 * with_doc / total if total else 0
        # Informational — print coverage but don't fail
        print(f"\nClass docstring coverage: {with_doc}/{total} ({pct:.1f}%)")
        # Soft threshold: at least 30% of classes should have docstrings
        assert with_doc >= total * 0.3, (
            f"Only {with_doc}/{total} classes ({pct:.1f}%) have docstrings. "
            f"Expected at least 30%."
        )

    def test_method_docstring_coverage(self):
        """Report how many methods have docstrings."""
        total = 0
        with_doc = 0
        for name in ALL_CLASSES:
            cls = getattr(pyopenms, name)
            for mname in _get_public_methods(cls):
                total += 1
                doc = getattr(getattr(cls, mname, None), "__doc__", None)
                if doc:
                    with_doc += 1
        pct = 100.0 * with_doc / total if total else 0
        print(f"\nMethod docstring coverage: {with_doc}/{total} ({pct:.1f}%)")
        # Most nanobind methods get auto-docstrings from wrap-doc
        assert with_doc >= total * 0.8, (
            f"Only {with_doc}/{total} methods ({pct:.1f}%) have docstrings. "
            f"Expected at least 80%."
        )
