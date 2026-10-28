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

def _resolve(path):
    """The class behind a name from _get_all_classes()."""
    obj = pyopenms
    for part in path.split("."):
        obj = getattr(obj, part)
    # the DB singletons (ElementDB, ...) are exported as callable wrappers
    return getattr(obj, "_singleton_class", obj)


def _get_all_classes():
    """Return all public classes exported by pyopenms.

    Includes the classes nested in them ("ConsensusFeature.Ratio") and the DB
    singletons behind their callable wrappers.
    """
    names = []
    seen = set()

    def visit(path, cls):
        if id(cls) in seen:
            return
        seen.add(id(cls))
        names.append(path)
        for attr, value in vars(cls).items():
            # only classes defined here, not aliases of other classes
            if inspect.isclass(value) and value.__qualname__ == f"{cls.__qualname__}.{attr}":
                visit(f"{path}.{attr}", value)

    for name in dir(pyopenms):
        if name.startswith("_"):
            continue
        obj = _resolve(name)
        if inspect.isclass(obj):
            visit(name, obj)
    return sorted(names)


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


_FIELD_RE = re.compile(r"^:(param|type|returns?|rtype|raises?)\b")
_BULLET_RE = re.compile(r"^[-*] \S")
_DIRECTIVE_OPTION_RE = re.compile(r"^\s+:[\w-]+:")


def _is_paragraph_text(line):
    return (
        bool(line.strip())
        and not line[0].isspace()
        and not _FIELD_RE.match(line)
        and not _BULLET_RE.match(line)
        and not line.startswith("..")
    )


def _check_rst_block_separation(docstring):
    """Check that block markup is separated from the paragraph before it.

    Without the blank line reST folds a directive, a bullet list or a field
    list into the preceding paragraph, so ``.. code-block::`` or ``:return:``
    end up as plain text (#10599). docutils accepts this silently.
    """
    errors = []
    lines = inspect.cleandoc(docstring).split("\n")
    seen_field = False
    for i, line in enumerate(lines):
        prev = lines[i - 1] if i else ""
        if line.startswith(".. ") and prev.strip():
            errors.append(f"Line {i+1}: directive not preceded by a blank line")
        if line.startswith(".. code-block::"):
            j = i + 1
            while j < len(lines) and _DIRECTIVE_OPTION_RE.match(lines[j]):
                j += 1
            if j + 1 >= len(lines) or lines[j].strip() or not lines[j + 1][:1].isspace():
                errors.append(
                    f"Line {i+1}: code-block content must follow a blank line and be indented"
                )
        if _FIELD_RE.match(line):
            if not seen_field and _is_paragraph_text(prev):
                # wording avoids "field list", which _is_acceptable_error waives
                errors.append(f"Line {i+1}: :param:/:return: block not preceded by a blank line")
            seen_field = True
        elif _BULLET_RE.match(line) and not seen_field and _is_paragraph_text(prev):
            errors.append(f"Line {i+1}: bullet list not preceded by a blank line")
        if line.rstrip().endswith("\\n"):
            errors.append(f"Line {i+1}: literal '\\n' instead of a blank line")
    return errors


_GETTER_RE = re.compile(r"^(get|has|is)[A-Z_]")
_SETTER_RE = re.compile(r"^set[A-Z_]")


def _overload_docs(func):
    """The user docstrings of all overloads (without nanobind's signatures)."""
    sigs = getattr(func, "__nb_signature__", None)
    if sigs is not None:
        return [doc for _sig, doc, _defaults in sigs if doc]
    doc = inspect.getdoc(func)
    return [doc] if doc else []


def _accessor_doc_mismatch(name, doc):
    """First line of ``doc`` if it documents the opposite accessor of ``name``."""
    first = next((line.strip() for line in doc.split("\n") if line.strip()), "")
    if _GETTER_RE.match(name) and re.match(r"(?i)sets?\b", first):
        return first
    if _SETTER_RE.match(name) and re.match(r"(?i)(returns?|gets?)\b", first):
        return first
    return None


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


# Sphinx cross-reference roles. Python is the default domain, so these are
# normally written without the ``py:`` prefix (``:meth:`foo```), which the plain
# docutils parser used here does not know. Listed explicitly rather than
# accepting every unknown role, so that a misspelt role is still reported.
_SPHINX_ROLES = frozenset({
    # Python domain
    "mod", "func", "data", "const", "class", "meth", "attr", "exc", "obj",
    # standard cross-reference roles
    "ref", "doc", "term", "envvar", "option", "keyword", "token",
})

_UNKNOWN_ROLE_RE = re.compile(r'unknown interpreted text role "([^"]+)"')


def _is_acceptable_error(error, docstring):
    """Filter RST errors that are expected in generated docstrings."""
    err = error.lower()
    if "hyperlink" in err and "here" in docstring.lower():
        return True
    if "field list" in err:
        return True
    role = _UNKNOWN_ROLE_RE.search(err)
    if role:
        role_name = role.group(1)
        if role_name.startswith("py:"):
            role_name = role_name[3:]
        if role_name in _SPHINX_ROLES:
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
    errors.extend(_check_rst_block_separation(docstring))
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
        cls = _resolve(class_name)
        doc = cls.__doc__
        if doc is not None:
            assert isinstance(doc, str), (
                f"{class_name}.__doc__ is {type(doc).__name__}, expected str"
            )

    @pytest.mark.parametrize("class_name", ALL_CLASSES)
    def test_class_docstring_valid_rst(self, class_name):
        """Every class docstring must be valid RST."""
        cls = _resolve(class_name)
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
        cls = _resolve(class_name)
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
        cls = _resolve(class_name)
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

    @pytest.mark.parametrize("class_name", ALL_CLASSES)
    def test_accessor_docstrings_match_method(self, class_name):
        """A getter must not be documented as a setter, nor vice versa.

        The former pxd-to-nanobind generator attached each docstring to the
        method *following* it, so getRT() read "Checks if retention time is
        set" and setRT() "Returns the retention time" (#10599).
        """
        cls = _resolve(class_name)
        mismatches = []
        for mname in _get_public_methods(cls):
            for doc in _overload_docs(getattr(cls, mname, None)):
                first = _accessor_doc_mismatch(mname, doc)
                if first:
                    mismatches.append(f"{class_name}.{mname}: {first!r}")
        assert not mismatches, (
            "Docstrings describing a different method:\n" + "\n".join(mismatches)
        )


class TestDocstringRegressions:
    """Docstrings reported broken in #10599."""

    def test_peptideidentification_class_example_is_code_block(self):
        doc = inspect.cleandoc(pyopenms.PeptideIdentification.__doc__)
        assert "\n\n.. code-block:: python\n\n   pep_id = oms.PeptideIdentification()" in doc
        assert "\n\n- A list of peptide hits" in doc

    @pytest.mark.parametrize("method, summary", [
        ("getRT", "Returns the retention time of the precursor"),
        ("setRT", "Sets the retention time of the precursor"),
        ("hasRT", "Checks if retention time is set"),
        ("getMZ", "Returns the precursor m/z value"),
        ("setMZ", "Sets the precursor m/z value"),
        ("hasMZ", "Checks if m/z value is set"),
        ("getHits", "Returns all peptide hits"),
        ("insertHit", "Appends a peptide hit to the list"),
        ("setHits", "Sets all peptide hits at once"),
        ("getScoreType", "Returns the type of score"),
        ("setScoreType", "Sets the score type"),
        ("getSpectrumReference", "Get the spectrum reference"),
        ("setSpectrumReference", "Set the spectrum reference"),
    ])
    def test_peptideidentification_method_docstrings(self, method, summary):
        docs = _overload_docs(getattr(pyopenms.PeptideIdentification, method))
        assert docs and all(d.lstrip().startswith(summary) for d in docs), docs
        # summary and field list are separate paragraphs
        assert all("\n\n:" in d for d in docs if ":param" in d or ":return" in d), docs


class TestDocstringCoverage:
    """Summary statistics on docstring coverage (informational)."""

    def test_class_docstring_coverage(self):
        """Report how many classes have docstrings."""
        total = len(ALL_CLASSES)
        with_doc = sum(
            1 for name in ALL_CLASSES
            if _resolve(name).__doc__
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
            cls = _resolve(name)
            for mname in _get_public_methods(cls):
                total += 1
                doc = getattr(getattr(cls, mname, None), "__doc__", None)
                if doc:
                    with_doc += 1
        pct = 100.0 * with_doc / total if total else 0
        print(f"\nMethod docstring coverage: {with_doc}/{total} ({pct:.1f}%)")
        # Most nanobind methods get auto-generated signature docstrings
        assert with_doc >= total * 0.8, (
            f"Only {with_doc}/{total} methods ({pct:.1f}%) have docstrings. "
            f"Expected at least 80%."
        )
