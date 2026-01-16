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

Three types of tests:
1. Tests that validate .pxd source files (work without built pyopenms)
2. Tests that validate autowrap-generated .pyx files (primary CI validation)
3. Tests that validate runtime docstrings (require built pyopenms)

The generated .pyx file tests are the primary CI validation because:
- Generated files ARE in the build directory (unlike source .pxd files)
- Tests what users actually see (post-autowrap processing)
- Catches autowrap bugs that might mangle wrap-doc content
"""
import os
import re
import glob
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


def get_generated_pyx_dir():
    """
    Get the directory containing autowrap-generated .pyx files.

    These files are generated in the build directory during the pyopenms build
    and contain the actual docstrings that end up in the compiled module.

    Returns the directory path if found, None otherwise.
    """
    if not HAS_PYOPENMS:
        return None

    # The generated .pyx files are in the same directory as the pyopenms package
    # They're named _pyopenms_0.pyx, _pyopenms_1.pyx, etc. (or _pyopenms.pyx)
    try:
        pkg_dir = os.path.dirname(pyopenms.__file__)
        # Check for split files first (more common in production builds)
        pyx_files = glob.glob(os.path.join(pkg_dir, '_pyopenms_*.pyx'))
        if pyx_files:
            return pkg_dir
        # Check for single file
        single_pyx = os.path.join(pkg_dir, '_pyopenms.pyx')
        if os.path.exists(single_pyx):
            return pkg_dir
    except (AttributeError, TypeError):
        # Some pyopenms packaging configurations may not provide __file__
        # (e.g., frozen executables, zip imports). Treat as "no generated pyx dir".
        pass

    return None


def extract_docstrings_from_pyx(filepath):
    """
    Extract class and method docstrings from an autowrap-generated .pyx file.

    Returns a list of tuples: (name, docstring_type, docstring, line_number)
    where docstring_type is 'class' or 'method'.
    """
    docstrings = []

    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
        lines = content.split('\n')

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Match class definitions: "cdef class ClassName:" or "cdef class ClassName(Base):"
        class_match = re.match(r'^cdef\s+class\s+(\w+)(?:\s*\([^)]*\))?\s*:', stripped)
        if class_match:
            class_name = class_match.group(1)
            # Look for docstring on next non-empty line
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and '"""' in lines[j]:
                docstring, end_line = _extract_triple_quoted_string(lines, j)
                if docstring:
                    docstrings.append((class_name, 'class', docstring, i + 1))
                i = end_line
                continue

        # Match method definitions: "def method_name(self" or "def method_name(cls"
        method_match = re.match(r'^\s*def\s+(\w+)\s*\(\s*(?:self|cls)', stripped)
        if method_match:
            method_name = method_match.group(1)
            # Look for docstring on next non-empty line
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and '"""' in lines[j]:
                docstring, end_line = _extract_triple_quoted_string(lines, j)
                if docstring:
                    docstrings.append((method_name, 'method', docstring, i + 1))
                i = end_line
                continue

        i += 1

    return docstrings


def _extract_triple_quoted_string(lines, start_line):
    """
    Extract a triple-quoted string starting at the given line.

    Returns (docstring_content, end_line_index) or (None, start_line) if not found.
    """
    line = lines[start_line]

    # Find the opening """
    open_pos = line.find('"""')
    if open_pos == -1:
        return None, start_line

    # Check if it's a single-line docstring: """text"""
    after_open = line[open_pos + 3:]
    close_pos = after_open.find('"""')
    if close_pos != -1:
        # Single line docstring
        return after_open[:close_pos].strip(), start_line

    # Multi-line docstring
    doc_lines = [after_open]
    i = start_line + 1

    while i < len(lines):
        current_line = lines[i]
        close_pos = current_line.find('"""')
        if close_pos != -1:
            # Found closing quotes
            doc_lines.append(current_line[:close_pos])
            break
        doc_lines.append(current_line)
        i += 1

    # Join and clean up the docstring
    docstring = '\n'.join(doc_lines)
    # Remove common leading whitespace (dedent)
    docstring_lines = docstring.split('\n')
    if len(docstring_lines) > 1:
        # Find minimum indentation (ignoring empty lines)
        min_indent = float('inf')
        for dline in docstring_lines[1:]:  # Skip first line
            if dline.strip():
                indent = len(dline) - len(dline.lstrip())
                min_indent = min(min_indent, indent)
        if min_indent < float('inf'):
            docstring_lines = [docstring_lines[0]] + [
                dline[min_indent:] if len(dline) > min_indent else dline
                for dline in docstring_lines[1:]
            ]
        docstring = '\n'.join(docstring_lines)

    return docstring.strip(), i


def get_generated_pyx_files():
    """
    Get list of all autowrap-generated .pyx files.

    Returns list of file paths, or empty list if not found.
    """
    pyx_dir = get_generated_pyx_dir()
    if pyx_dir is None:
        return []

    # Get split files
    pyx_files = sorted(glob.glob(os.path.join(pyx_dir, '_pyopenms_*.pyx')))
    if pyx_files:
        return pyx_files

    # Fall back to single file
    single_pyx = os.path.join(pyx_dir, '_pyopenms.pyx')
    if os.path.exists(single_pyx):
        return [single_pyx]

    return []


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


def validate_wrap_doc_structure(content, _filename=""):
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

                # Check if this is a comment line (starts with #)
                if next_stripped.startswith('#'):
                    next_indent = len(next_line) - len(next_line.lstrip())

                    # '##' comments are not wrap-doc continuations
                    if next_stripped.startswith('##'):
                        break

                    # If indentation drops significantly (>4 spaces), this is likely
                    # commented-out code, not a wrap-doc continuation. End the block.
                    if next_indent < wrap_doc_indent - 4 and next_stripped != '#':
                        break

                    # Check for valid continuation formats
                    is_valid_continuation = (
                        next_stripped == '#' or
                        next_stripped.startswith('#  ') or
                        next_stripped.startswith('# wrap-')
                    )

                    is_single_space = (
                        next_stripped.startswith('# ') and
                        not next_stripped.startswith('#  ') and
                        not next_stripped.startswith('# wrap-')
                    )

                    is_no_space = (
                        not next_stripped.startswith('# ') and
                        next_stripped != '#'
                    )

                    # Flag format issues
                    if is_single_space or is_no_space:
                        errors.append((
                            next_line_num,
                            "wrap-doc continuation lines must start with "
                            "'#  ' (hash + two spaces)"
                        ))

                    # Check for slight indentation inconsistency (1-4 spaces)
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


# =============================================================================
# Tests that work without built pyopenms (validate .pxd source files)
# =============================================================================

class TestPxdDocumentation:
    """Tests that validate .pxd documentation format (no pyopenms build needed)."""

    def test_pxd_directory_exists(self):
        """Test that we can find the pxds directory."""
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found (test run from build directory?)")
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
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            if '# wrap-doc:' not in content:
                missing_essential.append(filename)

        for filename in optional_files:
            filepath = os.path.join(pxd_dir, filename)
            if not os.path.exists(filepath):
                continue
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            if '# wrap-doc:' not in content:
                missing_optional.append(filename)

        if missing_optional:
            warnings.warn(
                f"Classes missing wrap-doc (should be added): {missing_optional}",
                stacklevel=2
            )

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
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()

            docs = extract_wrap_doc_from_pxd(content)
            for line_num, doc in docs:
                errors = check_rst_code_blocks(doc)
                if errors:
                    issues.append(f"{filename}:{line_num}: {errors}")

        # Report first 10 issues
        assert not issues, (
            "Code blocks missing language specifier:\n" +
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
            with open(filepath, 'r', encoding='utf-8') as f:
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

            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()

            docs = extract_wrap_doc_from_pxd(content)
            for line_num, doc in docs:
                errors = check_common_issues(doc)
                backtick_errors = [e for e in errors if 'backtick' in e.lower()]
                if backtick_errors:
                    issues.append(f"{filename}:{line_num}: {backtick_errors}")

        assert not issues, "Unbalanced backticks:\n" + "\n".join(issues[:10])

    def test_wrap_doc_structure_valid(self):
        """
        Test that wrap-doc annotations have valid structure for autowrap parsing.

        This catches common issues like:
        - Missing space after '#' in continuation lines
        - Inconsistent indentation in multi-line wrap-doc blocks

        Test fails if any structural issues are found.
        """
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        all_errors = []
        pxd_files = glob.glob(os.path.join(pxd_dir, '*.pxd'))

        for filepath in pxd_files:
            filename = os.path.basename(filepath)
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()

            errors = validate_wrap_doc_structure(content, filename)
            for line_num, error in errors:
                all_errors.append(f"{filename}:{line_num}: {error}")

        # Fail if any structural issues found
        assert not all_errors, (
            f"wrap-doc structure issues in {len(all_errors)} locations:\n" +
            "\n".join(all_errors[:50])
        )

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
            with open(filepath, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    for pattern in typo_patterns:
                        if re.search(pattern, line, re.IGNORECASE):
                            typos_found.append(f"{filename}:{line_num}: {line.strip()}")

        assert not typos_found, (
            "Possible wrap-doc typos found:\n" + "\n".join(typos_found[:10])
        )

    def test_wrap_doc_continuation_format(self):
        """
        Test that multi-line wrap-doc continuation lines have proper format.

        Autowrap expects continuation lines to start with '#  ' (hash + two spaces).
        Lines starting with '# ' (single space) or '#text' (no space) will break parsing.
        """
        pxd_dir = get_pxd_dir()
        if pxd_dir is None:
            pytest.skip("pxds directory not found")

        issues = []
        pxd_files = glob.glob(os.path.join(pxd_dir, '*.pxd'))

        for filepath in pxd_files:
            filename = os.path.basename(filepath)
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()

            errors = validate_wrap_doc_structure(content, filename)
            for line_num, error in errors:
                if 'continuation lines must start' in error:
                    issues.append(f"{filename}:{line_num}: {error}")

        assert not issues, (
            f"wrap-doc continuation lines with incorrect format ({len(issues)} issues):\n" +
            "\n".join(issues[:50])
        )


# =============================================================================
# Tests that validate autowrap-generated .pyx files (requires built pyopenms)
# =============================================================================

@pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not available")
class TestGeneratedDocstrings:
    """
    Tests that validate docstrings in autowrap-generated .pyx files.

    These tests check the actual generated Cython code that autowrap produces,
    which contains the docstrings that end up in the compiled module.

    This is the primary validation for CI because:
    1. Generated files ARE in the build directory (unlike source .pxd files)
    2. Tests what users actually see (post-autowrap processing)
    3. Catches autowrap bugs that might mangle wrap-doc content
    4. Validates docstrings from ALL sources (wrap-doc, addons, etc.)
    """

    def test_generated_pyx_files_exist(self):
        """Test that autowrap-generated .pyx files can be found."""
        pyx_files = get_generated_pyx_files()
        if not pyx_files:
            pytest.skip(
                "Generated .pyx files not found - they may be cleaned after compilation"
            )
        assert len(pyx_files) > 0, "Expected at least one generated .pyx file"

    def test_generated_files_have_docstrings(self):
        """Test that generated .pyx files contain docstrings."""
        pyx_files = get_generated_pyx_files()
        if not pyx_files:
            pytest.skip("Generated .pyx files not found")

        total_docstrings = 0
        for pyx_file in pyx_files:
            docstrings = extract_docstrings_from_pyx(pyx_file)
            total_docstrings += len(docstrings)

        # We expect many docstrings across all generated files
        assert total_docstrings > 100, (
            f"Expected many docstrings in generated files, found {total_docstrings}"
        )

    def test_core_classes_in_generated_files(self):
        """Test that core classes have docstrings in generated files."""
        pyx_files = get_generated_pyx_files()
        if not pyx_files:
            pytest.skip("Generated .pyx files not found")

        # Core classes that must have docstrings
        core_classes = {
            'MSSpectrum', 'MSChromatogram', 'MSExperiment',
            'Feature', 'FeatureMap', 'AASequence', 'Param',
        }

        found_classes = set()
        for pyx_file in pyx_files:
            docstrings = extract_docstrings_from_pyx(pyx_file)
            for name, doc_type, _docstring, _ in docstrings:
                if doc_type == 'class' and name in core_classes:
                    found_classes.add(name)

        missing = core_classes - found_classes
        assert not missing, f"Core classes missing in generated files: {missing}"

    def test_generated_docstrings_valid_rst(self):
        """Test that all docstrings in generated files are valid RST."""
        pyx_files = get_generated_pyx_files()
        if not pyx_files:
            pytest.skip("Generated .pyx files not found")

        all_errors = []

        # Check all docstrings from each file
        for pyx_file in pyx_files:
            filename = os.path.basename(pyx_file)
            docstrings = extract_docstrings_from_pyx(pyx_file)

            for name, doc_type, docstring, line_num in docstrings:
                is_valid, errors = validate_docstring(docstring, name)
                if not is_valid:
                    # Filter out some known acceptable patterns
                    real_errors = [
                        e for e in errors
                        if not _is_acceptable_error(e, docstring)
                    ]
                    if real_errors:
                        all_errors.append(
                            f"{filename}:{line_num} ({doc_type} {name}): {real_errors[:3]}"
                        )

        # Report first 100 issues (may be many if there's a systemic problem)
        assert not all_errors, (
            f"Generated docstrings with RST issues ({len(all_errors)} total):\n" +
            "\n".join(all_errors[:100])
        )

    def test_generated_docstrings_no_tabs(self):
        """Test that generated docstrings use spaces, not tabs."""
        pyx_files = get_generated_pyx_files()
        if not pyx_files:
            pytest.skip("Generated .pyx files not found")

        files_with_tabs = []
        for pyx_file in pyx_files:
            filename = os.path.basename(pyx_file)
            docstrings = extract_docstrings_from_pyx(pyx_file)

            for name, doc_type, docstring, line_num in docstrings:
                if '\t' in docstring:
                    files_with_tabs.append(f"{filename}:{line_num} ({doc_type} {name})")
                    break  # One per file is enough

        assert not files_with_tabs, f"Files with tabs in docstrings: {files_with_tabs}"

    def test_generated_docstrings_are_strings(self):
        """Test that extracted docstrings are proper strings, not bytes."""
        pyx_files = get_generated_pyx_files()
        if not pyx_files:
            pytest.skip("Generated .pyx files not found")

        for pyx_file in pyx_files:
            docstrings = extract_docstrings_from_pyx(pyx_file)
            for name, doc_type, docstring, line_num in docstrings:
                assert isinstance(docstring, str), (
                    f"Docstring for {doc_type} {name} at line {line_num} "
                    f"is {type(docstring)}, expected str"
                )

    def test_generated_docstrings_not_empty_placeholders(self):
        """Test that class docstrings have meaningful content, not just placeholders."""
        pyx_files = get_generated_pyx_files()
        if not pyx_files:
            pytest.skip("Generated .pyx files not found")

        empty_docs = []
        for pyx_file in pyx_files:
            filename = os.path.basename(pyx_file)
            docstrings = extract_docstrings_from_pyx(pyx_file)

            for name, doc_type, docstring, line_num in docstrings:
                if doc_type == 'class':
                    # Check if docstring is essentially empty or just whitespace
                    content = docstring.strip()
                    if not content or content == '""':
                        empty_docs.append(f"{filename}:{line_num} ({name})")

        # Allow some empty docs (internal classes) but not too many
        assert len(empty_docs) <= 50, (
            f"Too many classes with empty docstrings ({len(empty_docs)}): "
            f"{empty_docs[:10]}"
        )


def _is_acceptable_error(error, docstring):
    """
    Check if an RST validation error is acceptable for generated docstrings.

    Some patterns in generated docstrings are technically RST warnings but
    are acceptable in practice (e.g., autowrap's standard format).
    """
    # Autowrap generates "Original C++ documentation is available `here <url>`_"
    # which is valid RST but may trigger warnings in some validators
    if 'hyperlink' in error.lower() and 'here' in docstring.lower():
        return True

    # Method signatures as first line (e.g., "method(self, arg) -> Type")
    # may trigger field list warnings
    if 'field list' in error.lower():
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
                    all_errors.append(f"{class_name}: {errors}")

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
        """Test that incorrect continuation format is detected."""
        content = '''
cdef extern from "Test.h":
    cdef cppclass TestClass:
        # wrap-doc:
        #This line has no space after hash
        TestClass() except + nogil
'''
        errors = validate_wrap_doc_structure(content, "test.pxd")
        assert any('continuation lines must start' in e for _, e in errors), (
            f"Should detect incorrect continuation format: {errors}"
        )

    def test_single_space_after_hash(self):
        """Test that single space after # (instead of two) is detected."""
        content = '''
cdef extern from "Test.h":
    cdef cppclass TestClass:
        # wrap-doc:
        # This line has only one space after hash
        TestClass() except + nogil
'''
        errors = validate_wrap_doc_structure(content, "test.pxd")
        assert any('continuation lines must start' in e for _, e in errors), (
            f"Should detect single-space continuation: {errors}"
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


class TestPyxDocstringExtraction:
    """Unit tests for the .pyx docstring extraction functions."""

    def test_extract_class_docstring(self):
        """Test extraction of class docstrings from Cython code."""
        pyx_content = '''
cdef class MSSpectrum:
    """
    Cython implementation of MSSpectrum

    Original C++ documentation is available `here <http://example.com>`_
    """

    def __init__(self):
        pass
'''
        # Write to temp file and extract
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pyx', delete=False) as f:
            f.write(pyx_content)
            temp_path = f.name

        try:
            docstrings = extract_docstrings_from_pyx(temp_path)
            class_docs = [(n, d) for n, t, d, _ in docstrings if t == 'class']
            assert len(class_docs) == 1, f"Expected 1 class docstring, got {len(class_docs)}"
            name, doc = class_docs[0]
            assert name == 'MSSpectrum'
            assert 'Cython implementation' in doc
        finally:
            os.unlink(temp_path)

    def test_extract_method_docstring(self):
        """Test extraction of method docstrings from Cython code."""
        pyx_content = '''
cdef class TestClass:
    """Class docstring."""

    def get_value(self):
        """
        get_value(self) -> int

        Returns the current value.
        """
        return self._value

    def set_value(self, val):
        """
        set_value(self, val: int) -> None

        Sets the value.
        """
        self._value = val
'''
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pyx', delete=False) as f:
            f.write(pyx_content)
            temp_path = f.name

        try:
            docstrings = extract_docstrings_from_pyx(temp_path)
            method_docs = [(n, d) for n, t, d, _ in docstrings if t == 'method']
            assert len(method_docs) == 2, f"Expected 2 method docstrings, got {len(method_docs)}"
            names = [n for n, _ in method_docs]
            assert 'get_value' in names
            assert 'set_value' in names
        finally:
            os.unlink(temp_path)

    def test_extract_single_line_docstring(self):
        """Test extraction of single-line docstrings."""
        pyx_content = '''
cdef class Simple:
    """A simple class."""

    def method(self):
        """Does something."""
        pass
'''
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pyx', delete=False) as f:
            f.write(pyx_content)
            temp_path = f.name

        try:
            docstrings = extract_docstrings_from_pyx(temp_path)
            assert len(docstrings) >= 1
            class_doc = next((d for n, t, d, _ in docstrings if t == 'class'), None)
            assert class_doc == "A simple class."
        finally:
            os.unlink(temp_path)

    def test_extract_multiline_docstring_dedent(self):
        """Test that multi-line docstrings are properly dedented."""
        pyx_content = '''
cdef class TestClass:
    """
    First line.

    Second paragraph with more detail.
    And another line.
    """
    pass
'''
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pyx', delete=False) as f:
            f.write(pyx_content)
            temp_path = f.name

        try:
            docstrings = extract_docstrings_from_pyx(temp_path)
            class_docs = [d for n, t, d, _ in docstrings if t == 'class']
            assert len(class_docs) == 1
            doc = class_docs[0]
            # Check that it's dedented (no leading spaces on content lines)
            lines = doc.split('\n')
            for line in lines:
                if line.strip():
                    # Content lines should start at position 0 or have consistent indent
                    assert not line.startswith('        '), (
                        f"Docstring not properly dedented: {repr(doc)}"
                    )
        finally:
            os.unlink(temp_path)

    def test_extract_class_with_base(self):
        """Test extraction from class with base class."""
        pyx_content = '''
cdef class DerivedClass(BaseClass):
    """
    A derived class.

    Inherits from BaseClass.
    """
    pass
'''
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pyx', delete=False) as f:
            f.write(pyx_content)
            temp_path = f.name

        try:
            docstrings = extract_docstrings_from_pyx(temp_path)
            class_docs = [(n, d) for n, t, d, _ in docstrings if t == 'class']
            assert len(class_docs) == 1
            name, doc = class_docs[0]
            assert name == 'DerivedClass'
            assert 'derived class' in doc.lower()
        finally:
            os.unlink(temp_path)

    def test_extract_empty_file(self):
        """Test extraction from empty file returns empty list."""
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pyx', delete=False) as f:
            f.write('')
            temp_path = f.name

        try:
            docstrings = extract_docstrings_from_pyx(temp_path)
            assert docstrings == []
        finally:
            os.unlink(temp_path)

    def test_extract_no_docstrings(self):
        """Test extraction from file with no docstrings."""
        pyx_content = '''
cdef class NoDoc:
    def method(self):
        return 42
'''
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pyx', delete=False) as f:
            f.write(pyx_content)
            temp_path = f.name

        try:
            docstrings = extract_docstrings_from_pyx(temp_path)
            # Should find no docstrings (no triple-quoted strings after class/def)
            assert len(docstrings) == 0
        finally:
            os.unlink(temp_path)


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
