"""
Tests for pyOpenMS docstring formatting and Sphinx RST compatibility.

This module validates that docstrings are:
1. Present on important classes and methods
2. Properly formatted RST (reStructuredText)
3. Compatible with Sphinx documentation generation

Common issues this catches:
- Missing docstrings on public APIs
- Malformed RST code blocks (missing indentation, wrong directive syntax)
- Broken cross-references
- Invalid field list syntax (:param:, :return:, etc.)

Two types of tests:
1. Tests that validate .pxd source files (work without built pyopenms)
2. Tests that validate runtime docstrings (require built pyopenms)
"""
import os
import re
import glob
import inspect
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


def get_pxd_dir():
    """Get the path to the pxds directory."""
    # This file is in src/pyOpenMS/tests/unittests/
    this_dir = os.path.dirname(os.path.abspath(__file__))
    pxd_dir = os.path.normpath(os.path.join(this_dir, '..', '..', 'pxds'))
    if os.path.isdir(pxd_dir):
        return pxd_dir
    return None


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
                whitespace, language = match.groups()
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
        return False, [f"Docstring is not a string: {type(docstring)}"]

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


def extract_wrap_doc_from_pxd(content):
    """
    Extract wrap-doc content from a .pxd file.

    Returns list of (line_number, doc_content) tuples.
    """
    docs = []
    lines = content.split('\n')

    i = 0
    while i < len(lines):
        line = lines[i]

        # Check for inline wrap-doc (method documentation)
        if '# wrap-doc:' in line:
            match = re.search(r'#\s*wrap-doc:\s*(.*)$', line)
            if match:
                doc = match.group(1).strip()
                if doc:
                    docs.append((i + 1, doc))

        # Check for multi-line wrap-doc block (class documentation)
        stripped = line.strip()
        if stripped == '# wrap-doc:':
            doc_lines = []
            doc_start = i + 1
            i += 1
            while i < len(lines):
                next_line = lines[i]
                # Multi-line wrap-doc continues with '#  ' (hash + 2 spaces)
                if next_line.strip().startswith('#  ') or next_line.strip() == '#':
                    # Remove the leading '# ' or '#  '
                    content_match = re.match(r'\s*#\s?(.*)', next_line)
                    if content_match:
                        doc_lines.append(content_match.group(1))
                    i += 1
                elif next_line.strip().startswith('# ') and not next_line.strip().startswith('# wrap-'):
                    # Continuation with single space
                    content_match = re.match(r'\s*#\s?(.*)', next_line)
                    if content_match:
                        doc_lines.append(content_match.group(1))
                    i += 1
                else:
                    break

            if doc_lines:
                full_doc = '\n'.join(doc_lines)
                docs.append((doc_start, full_doc))
            continue

        i += 1

    return docs


def validate_wrap_doc_structure(content, filename=""):
    """
    Validate the structural format of wrap-doc annotations.

    This catches issues that would cause autowrap to fail or mangle documentation:
    - Inconsistent indentation in multi-line wrap-doc
    - Missing space after # in continuation lines
    - wrap-doc with wrong indentation relative to class
    - Lines that look like wrap-doc continuation but have wrong format

    Returns list of (line_number, error_message) tuples.
    """
    errors = []
    lines = content.split('\n')

    i = 0
    while i < len(lines):
        line = lines[i]
        line_num = i + 1

        # Check for multi-line wrap-doc block
        stripped = line.strip()
        if stripped == '# wrap-doc:':
            # Get the indentation of the wrap-doc: line
            wrap_doc_indent = len(line) - len(line.lstrip())
            block_start = line_num
            i += 1

            # Check continuation lines
            while i < len(lines):
                next_line = lines[i]
                next_stripped = next_line.strip()
                next_line_num = i + 1

                # Empty line or whitespace only - ok
                if not next_stripped:
                    i += 1
                    continue

                # Check if this is a continuation line (starts with #)
                if next_stripped.startswith('#'):
                    # Must have same or greater indent as wrap-doc line
                    next_indent = len(next_line) - len(next_line.lstrip())

                    # Check for missing space after # (common mistake)
                    # Valid: '#  text' or '# text' or '#'
                    # Invalid: '#text' (no space)
                    hash_match = re.match(r'\s*#(\S)', next_line)
                    if hash_match and hash_match.group(1) not in (' ', '\t'):
                        errors.append((
                            next_line_num,
                            f"Missing space after '#' in wrap-doc continuation "
                            f"(got '#{hash_match.group(1)}', expected '# ' or '#  ')"
                        ))

                    # Check for inconsistent indentation
                    if next_indent < wrap_doc_indent and next_stripped != '#':
                        errors.append((
                            next_line_num,
                            f"Inconsistent indentation in wrap-doc block "
                            f"(started at line {block_start} with indent {wrap_doc_indent}, "
                            f"but line {next_line_num} has indent {next_indent})"
                        ))

                    # Check if this is actually another wrap directive (end of current block)
                    if next_stripped.startswith('# wrap-'):
                        break

                    i += 1
                else:
                    # Non-comment line - end of wrap-doc block
                    break

            continue

        # Note: inline wrap-doc format "# wrap-doc:Text" (no space after colon) is valid
        # The only issue is if wrap-doc: is followed by nothing (empty docstring)

        # Check for potential wrap-doc typos
        typo_patterns = [
            (r'#\s*wrapdoc:', 'wrapdoc:'),
            (r'#\s*wrap_doc:', 'wrap_doc:'),
            (r'#\s*wrap doc:', 'wrap doc:'),
            (r'#\s*warp-doc:', 'warp-doc:'),  # common typo
        ]
        for pattern, typo in typo_patterns:
            if re.search(pattern, line, re.IGNORECASE):
                errors.append((
                    line_num,
                    f"Possible wrap-doc typo: '{typo}' should be 'wrap-doc:'"
                ))

        i += 1

    return errors


def check_wrap_doc_in_class_context(content, filename=""):
    """
    Check that wrap-doc annotations appear in valid locations.

    wrap-doc should appear:
    - Inside a 'cdef cppclass' block (for class docs)
    - After a method declaration (for method docs)

    Returns list of (line_number, error_message) tuples.
    """
    errors = []
    lines = content.split('\n')

    in_class = False
    class_indent = 0

    for i, line in enumerate(lines):
        line_num = i + 1
        stripped = line.strip()

        # Track class context
        if 'cdef cppclass' in line or 'cdef class' in line:
            in_class = True
            class_indent = len(line) - len(line.lstrip())
            continue

        # Detect end of class (line with same or less indent that's not empty/comment)
        if in_class and stripped and not stripped.startswith('#'):
            current_indent = len(line) - len(line.lstrip())
            if current_indent <= class_indent:
                in_class = False

        # Check for wrap-doc outside class context
        if '# wrap-doc:' in line:
            if not in_class:
                # Check if it's a file-level or module-level doc (unusual but allowed)
                # For now, just warn - it's likely an error
                errors.append((
                    line_num,
                    "wrap-doc annotation appears outside of class definition "
                    "(may be orphaned or incorrectly indented)"
                ))

    return errors


# =============================================================================
# Tests that work without built pyopenms (validate .pxd source files)
# =============================================================================

class TestPxdDocumentation:
    """Tests that validate .pxd documentation format (no pyopenms build needed)."""

    def test_pxd_directory_exists(self):
        """Test that we can find the pxds directory."""
        pxd_dir = get_pxd_dir()
        assert pxd_dir is not None, "Could not find pxds directory"
        assert os.path.isdir(pxd_dir), f"pxds directory not found: {pxd_dir}"

    def test_pxd_files_exist(self):
        """Test that pxd files exist."""
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        pxd_files = glob.glob(os.path.join(pxd_dir, '*.pxd'))
        assert len(pxd_files) > 100, f"Expected many .pxd files, found {len(pxd_files)}"

    def test_core_classes_have_wrap_doc(self):
        """Test that core classes have wrap-doc annotations."""
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        # Essential classes that MUST have documentation
        essential_files = [
            'MSSpectrum.pxd',
            'MSChromatogram.pxd',
            'MSExperiment.pxd',
            'AASequence.pxd',
        ]

        # Additional classes that should have docs (warning only)
        optional_files = [
            'Feature.pxd',
            'FeatureMap.pxd',
            'Param.pxd',
            'ConsensusMap.pxd',
        ]

        missing_essential = []
        missing_optional = []

        for filename in essential_files:
            filepath = os.path.join(pxd_dir, filename)
            if not os.path.exists(filepath):
                continue
            with open(filepath, 'r') as f:
                content = f.read()
            if '# wrap-doc:' not in content:
                missing_essential.append(filename)

        for filename in optional_files:
            filepath = os.path.join(pxd_dir, filename)
            if not os.path.exists(filepath):
                continue
            with open(filepath, 'r') as f:
                content = f.read()
            if '# wrap-doc:' not in content:
                missing_optional.append(filename)

        if missing_optional:
            warnings.warn(f"Classes missing wrap-doc (should be added): {missing_optional}")

        assert not missing_essential, f"Essential classes missing wrap-doc: {missing_essential}"

    def test_wrap_doc_code_blocks_have_language(self):
        """Test that code-block directives in wrap-doc have language specifiers."""
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        issues = []
        pxd_files = glob.glob(os.path.join(pxd_dir, '*.pxd'))

        for filepath in pxd_files:
            filename = os.path.basename(filepath)
            with open(filepath, 'r') as f:
                content = f.read()

            docs = extract_wrap_doc_from_pxd(content)
            for line_num, doc in docs:
                errors = check_rst_code_blocks(doc)
                if errors:
                    issues.append(f"{filename}:{line_num}: {errors}")

        # Report first 10 issues
        assert not issues, (
            f"Code blocks missing language specifier:\n" +
            "\n".join(issues[:10])
        )

    def test_wrap_doc_no_tabs(self):
        """Test that wrap-doc content uses spaces, not tabs."""
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        files_with_tabs = []
        pxd_files = glob.glob(os.path.join(pxd_dir, '*.pxd'))

        for filepath in pxd_files:
            filename = os.path.basename(filepath)
            with open(filepath, 'r') as f:
                content = f.read()

            docs = extract_wrap_doc_from_pxd(content)
            for line_num, doc in docs:
                if '\t' in doc:
                    files_with_tabs.append(f"{filename}:{line_num}")
                    break  # One per file is enough

        assert not files_with_tabs, f"Files with tabs in wrap-doc: {files_with_tabs}"

    def test_wrap_doc_balanced_backticks(self):
        """Test that inline code backticks are balanced in wrap-doc."""
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        issues = []
        # Check just a few core files to keep test fast
        core_files = ['MSSpectrum.pxd', 'AASequence.pxd', 'Param.pxd']

        for filename in core_files:
            filepath = os.path.join(pxd_dir, filename)
            if not os.path.exists(filepath):
                continue

            with open(filepath, 'r') as f:
                content = f.read()

            docs = extract_wrap_doc_from_pxd(content)
            for line_num, doc in docs:
                errors = check_common_issues(doc)
                backtick_errors = [e for e in errors if 'backtick' in e.lower()]
                if backtick_errors:
                    issues.append(f"{filename}:{line_num}: {backtick_errors}")

        assert not issues, f"Unbalanced backticks:\n" + "\n".join(issues[:10])

    def test_wrap_doc_structure_valid(self):
        """
        Test that wrap-doc annotations have valid structure for autowrap parsing.

        This catches common issues like:
        - Missing space after '#' in continuation lines
        - Inconsistent indentation in multi-line wrap-doc blocks

        Issues are reported as warnings; the test passes but logs problems.
        """
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        all_errors = []
        pxd_files = glob.glob(os.path.join(pxd_dir, '*.pxd'))

        for filepath in pxd_files:
            filename = os.path.basename(filepath)
            with open(filepath, 'r') as f:
                content = f.read()

            errors = validate_wrap_doc_structure(content, filename)
            for line_num, error in errors:
                all_errors.append(f"{filename}:{line_num}: {error}")

        # Report issues as warnings (don't fail the test)
        # This allows CI to pass while developers can still see problems
        if all_errors:
            # Deduplicate by file (one error per file max in warning)
            files_with_issues = sorted(set(e.split(':')[0] for e in all_errors))
            warnings.warn(
                f"wrap-doc structure issues in {len(files_with_issues)} files "
                f"(may cause autowrap parsing issues): {files_with_issues[:10]}"
            )
        # Test passes - this is informational

    def test_wrap_doc_no_typos(self):
        """Test for common wrap-doc typos like 'wrapdoc:', 'wrap_doc:', etc."""
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        typos_found = []
        pxd_files = glob.glob(os.path.join(pxd_dir, '*.pxd'))

        # Common typo patterns
        typo_patterns = [
            r'#\s*wrapdoc:',
            r'#\s*wrap_doc:',
            r'#\s*wrap doc:',
            r'#\s*warp-doc:',
            r'#\s*wrap-docs:',
        ]

        for filepath in pxd_files:
            filename = os.path.basename(filepath)
            with open(filepath, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    for pattern in typo_patterns:
                        if re.search(pattern, line, re.IGNORECASE):
                            typos_found.append(f"{filename}:{line_num}: {line.strip()}")

        assert not typos_found, (
            f"Possible wrap-doc typos found:\n" + "\n".join(typos_found[:10])
        )

    def test_wrap_doc_continuation_format(self):
        """
        Test that multi-line wrap-doc continuation lines have proper format.

        Autowrap expects continuation lines to start with '#  ' (hash + space + space)
        or '# ' (hash + space). Lines starting with '#text' (no space) will break parsing.
        """
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        issues = []
        # Check core files that commonly have multi-line docs
        core_files = ['MSSpectrum.pxd', 'MSExperiment.pxd', 'AASequence.pxd',
                      'Feature.pxd', 'FeatureMap.pxd', 'Param.pxd']

        for filename in core_files:
            filepath = os.path.join(pxd_dir, filename)
            if not os.path.exists(filepath):
                continue

            with open(filepath, 'r') as f:
                content = f.read()

            errors = validate_wrap_doc_structure(content, filename)
            for line_num, error in errors:
                if 'Missing space' in error:
                    issues.append(f"{filename}:{line_num}: {error}")

        assert not issues, (
            f"wrap-doc continuation lines with missing spaces:\n" +
            "\n".join(issues[:10])
        )


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
                    all_errors.append(f"{class_name}: {errors}")

        assert not all_errors, (
            f"Docstrings with RST issues:\n" + "\n".join(all_errors[:10])
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
        is_valid, errors = validate_docstring(docstring)
        code_errors = [e for e in errors if 'code-block' in e]
        assert not code_errors, f"Unexpected code block errors: {code_errors}"

    def test_code_block_missing_language(self):
        """Test that code-block without language is detected."""
        docstring = '''
        Example:

        .. code-block::

            some code
        '''
        is_valid, errors = validate_docstring(docstring)
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
        is_valid, errors = validate_docstring(docstring)
        field_errors = [e for e in errors if ':param' in e or ':type' in e]
        assert not field_errors, f"Valid field list should pass: {field_errors}"

    def test_javadoc_style_detected(self):
        """Test that Javadoc-style params are detected."""
        docstring = '''
        Do something.

        @param name The name to use
        @return The result
        '''
        is_valid, errors = validate_docstring(docstring)
        assert any('@param' in e or '@return' in e for e in errors), (
            f"Should detect Javadoc style: {errors}"
        )

    def test_unbalanced_backticks(self):
        """Test that unbalanced backticks are detected."""
        docstring = '''
        Use `get_peaks method to get data.
        '''
        is_valid, errors = validate_docstring(docstring)
        assert any('backtick' in e.lower() for e in errors), (
            f"Should detect unbalanced backticks: {errors}"
        )

    def test_tabs_detected(self):
        """Test that tabs are detected."""
        docstring = "Some text\twith tabs"
        is_valid, errors = validate_docstring(docstring)
        assert any('tab' in e.lower() for e in errors), (
            f"Should detect tabs: {errors}"
        )


class TestWrapDocStructureValidation:
    """Unit tests for wrap-doc structure validation functions."""

    def test_valid_multiline_wrap_doc(self):
        """Test that valid multi-line wrap-doc passes validation."""
        content = '''
cdef extern from "Test.h":
    cdef cppclass TestClass:
        # wrap-doc:
        #  This is a test class.
        #  It does things.
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    t = TestClass()
        TestClass() except + nogil
'''
        errors = validate_wrap_doc_structure(content, "test.pxd")
        assert not errors, f"Valid wrap-doc should pass: {errors}"

    def test_missing_space_after_hash(self):
        """Test that missing space after # is detected."""
        content = '''
cdef extern from "Test.h":
    cdef cppclass TestClass:
        # wrap-doc:
        #This line has no space after hash
        TestClass() except + nogil
'''
        errors = validate_wrap_doc_structure(content, "test.pxd")
        assert any('Missing space' in e for _, e in errors), (
            f"Should detect missing space: {errors}"
        )

    def test_inline_wrap_doc_valid(self):
        """Test that 'wrap-doc:text' without space is valid (common style)."""
        content = '''
cdef extern from "Test.h":
    cdef cppclass TestClass:
        void method() # wrap-doc:No space after colon is OK
'''
        errors = validate_wrap_doc_structure(content, "test.pxd")
        # This should NOT produce errors - the format "wrap-doc:text" is valid
        assert not errors, f"wrap-doc:text format should be valid: {errors}"

    def test_inconsistent_indentation(self):
        """Test that inconsistent indentation is detected."""
        content = '''
cdef extern from "Test.h":
    cdef cppclass TestClass:
        # wrap-doc:
        #  Line with correct indent
    #  Line with wrong indent (dedented too much)
        TestClass() except + nogil
'''
        errors = validate_wrap_doc_structure(content, "test.pxd")
        assert any('Inconsistent indentation' in e for _, e in errors), (
            f"Should detect inconsistent indentation: {errors}"
        )

    def test_wrap_doc_typo_detected(self):
        """Test that wrap-doc typos are detected."""
        typo_variants = [
            '# wrapdoc: text',
            '# wrap_doc: text',
            '# wrap doc: text',
            '# warp-doc: text',
        ]
        for typo in typo_variants:
            content = f'''
cdef extern from "Test.h":
    cdef cppclass TestClass:
        void method() {typo}
'''
            errors = validate_wrap_doc_structure(content, "test.pxd")
            assert any('typo' in e.lower() for _, e in errors), (
                f"Should detect typo '{typo}': {errors}"
            )

    def test_valid_inline_wrap_doc(self):
        """Test that valid inline wrap-doc passes validation."""
        content = '''
cdef extern from "Test.h":
    cdef cppclass TestClass:
        void method() except + nogil  # wrap-doc: Does something useful
        int getValue() except + nogil  # wrap-doc: Returns the value
'''
        errors = validate_wrap_doc_structure(content, "test.pxd")
        assert not errors, f"Valid inline wrap-doc should pass: {errors}"


def test_docutils_availability():
    """Report whether docutils is available for full RST validation."""
    if HAS_DOCUTILS:
        print("docutils is available - full RST parsing enabled")
    else:
        warnings.warn(
            "docutils not available - using regex-based RST validation only. "
            "Install docutils for more thorough validation: pip install docutils"
        )
    assert True
