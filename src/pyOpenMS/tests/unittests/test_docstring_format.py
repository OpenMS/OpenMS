"""
Tests for pyOpenMS docstring formatting and Sphinx RST compatibility.

Validates that runtime docstrings are:
1. Present on important classes and methods
2. Properly formatted RST (reStructuredText)
3. Compatible with Sphinx documentation generation

Common issues this catches:
- Missing docstrings on public APIs
- Malformed RST code blocks (missing indentation, wrong directive syntax)
- Broken cross-references
- Invalid field list syntax (:param:, :return:, etc.)
"""
import re
import warnings
import pytest


# Try to import docutils for full RST parsing
try:
    from docutils.parsers.rst import Parser
    from docutils.utils import Reporter, new_document
    from docutils.frontend import OptionParser
    HAS_DOCUTILS = True
except ImportError:
    HAS_DOCUTILS = False

# Try to import pyopenms
try:
    import pyopenms
    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False
    pyopenms = None


# =============================================================================
# Validation helpers
# =============================================================================

def get_docutils_errors(docstring):
    """
    Parse docstring with docutils and return any errors/warnings.

    Returns list of (level, message) tuples where level is:
    - 1: INFO
    - 2: WARNING
    - 3: ERROR
    - 4: SEVERE
    """
    if not HAS_DOCUTILS or not docstring:
        return []

    errors = []

    class ErrorCollector:
        def __init__(self):
            self.errors = []

        def __call__(self, message):
            level = message.get('level', 1)
            text = message.astext()
            self.errors.append((level, text))

    try:
        parser = Parser()
        option_parser = OptionParser(components=(Parser,))
        settings = option_parser.get_default_values()
        settings.report_level = Reporter.WARNING_LEVEL
        settings.halt_level = Reporter.SEVERE_LEVEL + 1

        document = new_document('<docstring>', settings)
        collector = ErrorCollector()
        document.reporter.attach_observer(collector)
        parser.parse(docstring, document)
        errors = collector.errors
    except Exception as e:
        errors.append((3, f"RST parsing failed: {e}"))

    return errors


def check_rst_code_blocks(docstring):
    """
    Check for common RST code block issues.

    Returns list of error messages.
    """
    if not docstring:
        return []

    errors = []
    lines = docstring.split('\n')

    for i, line in enumerate(lines, 1):
        stripped = line.strip()

        # Check for code-block directive without language
        if '.. code-block::' in stripped:
            match = re.search(r'\.\. code-block::(\s*)(\S*)', stripped)
            if match:
                _whitespace, language = match.groups()
                if not language:
                    errors.append(
                        f"Line {i}: code-block directive missing language specifier "
                        "(should be '.. code-block:: python')"
                    )

    return errors


def check_rst_field_lists(docstring):
    """
    Check for properly formatted RST field lists (:param:, :return:, etc).

    Returns list of error messages.
    """
    if not docstring:
        return []

    errors = []

    # Check for malformed fields using negative lookahead
    malformed_patterns = [
        (r':param\s+\w+(?!:)\s', 'Missing colon after :param name'),
        (r':type\s+\w+(?!:)\s', 'Missing colon after :type name'),
        (r':returns(?!:)', ':returns should be :returns: or :return:'),
        (r'@param\s+', 'Use :param: instead of @param (RST not Javadoc)'),
        (r'@return\s+', 'Use :return: instead of @return (RST not Javadoc)'),
    ]

    lines = docstring.split('\n')
    for i, line in enumerate(lines, 1):
        for pattern, message in malformed_patterns:
            if re.search(pattern, line):
                errors.append(f"Line {i}: {message}")

    return errors


def check_common_issues(docstring):
    """
    Check for common docstring issues.

    Returns list of error messages.
    """
    if not docstring:
        return []

    errors = []

    # Check for broken cross-references
    xref_pattern = r':(class|meth|func|attr|mod|obj):`([^`]*)`'
    for match in re.finditer(xref_pattern, docstring):
        ref_type, ref_name = match.groups()
        if not ref_name.strip():
            errors.append(f"Empty cross-reference :{ref_type}:``")
        if ref_name.count('`') > 0:
            errors.append(f"Unbalanced backticks in :{ref_type}:`{ref_name}`")

    # Check for unbalanced backticks in inline code
    lines = docstring.split('\n')
    for i, line in enumerate(lines, 1):
        # Count single backticks (not double)
        backticks = re.findall(r'(?<!`)`(?!`)', line)
        if len(backticks) % 2 != 0:
            errors.append(f"Line {i}: Unbalanced backticks (inline code)")

    # Check for tabs
    if '\t' in docstring:
        errors.append("Contains tabs - use spaces for RST indentation")

    return errors


def validate_docstring(docstring, name=""):
    """
    Validate a docstring for RST compatibility.

    Returns tuple of (is_valid, errors_list).
    """
    if docstring is None:
        return True, []

    if not isinstance(docstring, str):
        prefix = f"{name}: " if name else ""
        return False, [f"{prefix}Docstring is not a string: {type(docstring)}"]

    all_errors = []
    all_errors.extend(check_rst_code_blocks(docstring))
    all_errors.extend(check_rst_field_lists(docstring))
    all_errors.extend(check_common_issues(docstring))

    if HAS_DOCUTILS:
        docutils_errors = get_docutils_errors(docstring)
        for level, msg in docutils_errors:
            if level >= 2:
                severity = {2: 'WARNING', 3: 'ERROR', 4: 'SEVERE'}.get(level, 'UNKNOWN')
                all_errors.append(f"RST {severity}: {msg}")

    return len(all_errors) == 0, all_errors


def _is_acceptable_error(error, docstring):
    """
    Check if an RST validation error is acceptable for nanobind docstrings.

    Some patterns in generated docstrings are technically RST warnings but
    are acceptable in practice.
    """
    error_lower = error.lower()

    if 'hyperlink' in error_lower and 'here' in docstring.lower():
        return True
    if 'field list' in error_lower:
        return True
    if 'unknown interpreted text role "py:' in error_lower:
        return True
    if 'unknown directive type "py:' in error_lower:
        return True
    if 'unknown directive type "deprecated"' in error_lower:
        return True
    # NumPy-style docstring section headers (Parameters, Returns, Raises, etc.)
    # are valid for Sphinx napoleon but not for raw docutils RST parsing
    if "unexpected section title" in error_lower:
        return True
    # nanobind auto-generated signatures and NumPy-style **kwargs may contain
    # unmatched * or ** that RST interprets as emphasis/bold markers
    if "inline emphasis start-string without end-string" in error_lower:
        return True
    if "inline strong start-string without end-string" in error_lower:
        return True

    return False


# =============================================================================
# Tests that require built pyopenms (validate runtime docstrings)
# =============================================================================

@pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not available")
class TestRuntimeDocstrings:
    """Tests that validate runtime docstrings (requires built pyopenms)."""

    def test_pyopenms_import(self):
        """Test that pyopenms can be imported."""
        assert pyopenms is not None

    def test_core_classes_have_docstrings(self):
        """Test that core classes have docstrings."""
        core_classes = [
            'MSSpectrum',
            'MSChromatogram',
            'MSExperiment',
            'Feature',
            'FeatureMap',
            'ConsensusMap',
            'ConsensusFeature',
            'AASequence',
            'Param',
        ]

        missing_docstrings = []
        for class_name in core_classes:
            if not hasattr(pyopenms, class_name):
                continue
            cls = getattr(pyopenms, class_name)
            if not cls.__doc__:
                missing_docstrings.append(class_name)

        assert not missing_docstrings, (
            f"Core classes missing docstrings: {missing_docstrings}"
        )

    def test_docstrings_are_strings(self):
        """Test that docstrings are str, not bytes."""
        classes_to_check = ['MSSpectrum', 'String', 'AASequence']

        binary_docstrings = []
        for class_name in classes_to_check:
            if not hasattr(pyopenms, class_name):
                continue
            cls = getattr(pyopenms, class_name)
            if cls.__doc__ is not None and isinstance(cls.__doc__, bytes):
                binary_docstrings.append(class_name)

        assert not binary_docstrings, (
            f"Classes with binary docstrings: {binary_docstrings}"
        )

    def test_docstrings_valid_rst(self):
        """Test that runtime docstrings are valid RST."""
        classes_to_check = ['MSSpectrum', 'MSExperiment', 'AASequence', 'Feature']

        all_errors = []
        for class_name in classes_to_check:
            if not hasattr(pyopenms, class_name):
                continue
            cls = getattr(pyopenms, class_name)
            if cls.__doc__:
                is_valid, errors = validate_docstring(cls.__doc__, class_name)
                if not is_valid:
                    # Filter out acceptable errors (e.g., Sphinx-specific roles)
                    real_errors = [
                        e for e in errors
                        if not _is_acceptable_error(e, cls.__doc__)
                    ]
                    if real_errors:
                        all_errors.append(f"{class_name}: {real_errors}")

        assert not all_errors, (
            "Docstrings with RST issues:\n" + "\n".join(all_errors[:10])
        )


# =============================================================================
# Unit tests for validation functions
# =============================================================================

class TestDocstringValidation:
    """Unit tests for the docstring validation functions."""

    def test_valid_rst_code_block(self):
        """Test that valid RST code blocks pass validation."""
        docstring = '''
        Example usage:

        .. code-block:: python

            from pyopenms import MSSpectrum
            s = MSSpectrum()
        '''
        _is_valid, errors = validate_docstring(docstring)
        code_errors = [e for e in errors if 'code-block' in e]
        assert not code_errors, f"Unexpected code block errors: {code_errors}"

    def test_code_block_missing_language(self):
        """Test that code-block without language is detected."""
        docstring = '''
        Example:

        .. code-block::

            some code
        '''
        _is_valid, errors = validate_docstring(docstring)
        assert any('missing language' in e for e in errors), (
            f"Should detect missing language: {errors}"
        )

    def test_valid_field_list(self):
        """Test that valid field lists pass validation."""
        docstring = '''
        Do something.

        :param name: The name to use
        :type name: str
        :return: The result
        :rtype: int
        '''
        _is_valid, errors = validate_docstring(docstring)
        field_errors = [e for e in errors if ':param' in e or ':type' in e]
        assert not field_errors, f"Valid field list should pass: {field_errors}"

    def test_javadoc_style_detected(self):
        """Test that Javadoc-style params are detected."""
        docstring = '''
        Do something.

        @param name The name to use
        @return The result
        '''
        _is_valid, errors = validate_docstring(docstring)
        assert any('@param' in e or '@return' in e for e in errors), (
            f"Should detect Javadoc style: {errors}"
        )

    def test_unbalanced_backticks(self):
        """Test that unbalanced backticks are detected."""
        docstring = '''
        Use `get_peaks method to get data.
        '''
        _is_valid, errors = validate_docstring(docstring)
        assert any('backtick' in e.lower() for e in errors), (
            f"Should detect unbalanced backticks: {errors}"
        )

    def test_tabs_detected(self):
        """Test that tabs are detected."""
        docstring = "Some text\twith tabs"
        _is_valid, errors = validate_docstring(docstring)
        assert any('tab' in e.lower() for e in errors), (
            f"Should detect tabs: {errors}"
        )


def test_docutils_availability():
    """Report whether docutils is available for full RST validation."""
    if HAS_DOCUTILS:
        print("docutils is available - full RST parsing enabled")
    else:
        warnings.warn(
            "docutils not available - using regex-based RST validation only. "
            "Install docutils for more thorough validation: pip install docutils",
            stacklevel=2
        )
    assert True
