# pyTOPP tests

This directory contains unit tests and light integration tests for the `pyTOPP` package (the Python wrappers living under `pyopenms.pytopp`). Tests are written with **pytest** and are also wired into **CTest** via the OpenMS's project's CMake setup.

---

## Layout

- `test_util.py` - unit tests for general helpers in `pyopenms.pytopp.util`
- `test_ctd_support.py` - unit tests for CTD help/printer and parsing helpers in `pyopenms.pytopp.ctdsupport`
- `test_integration_pytopp_pyprophet.py` - light integration tests for the PyProphet CLI wrapper
- `conftest.py` - small bootstrapping for tests (e.g. shims/stubs so you can run tests without a full pyopenms install)
- `data` *(optional)* - small test files for integration tests (e.g., tiny OSW files)

---

## Quick start (pytest)

From the repository root (or anywhere), run:

```bash
pytest -q src/pyTOPP/tests
```

Tips:

* Show verbose output: `pytest -vv src/pyTOPP/tests`
* Run a single file: `pytest -q src/pyTOPP/tests/test_util.py`
* Run a single test: `pytest -q src/pyTOPP/tests/test_util.py::test_as_bool_variants`
* See print statements while debugging: `pytest -s …`
* Filter by keyword: `pytest -k levels src/pyTOPP/tests`
* Run tests in parallel (requires `pytest-xdist`): `pytest -n auto src/pyTOPP/tests` or `pytest -n 4 src/pyTOPP/tests`

> **Dev env note:** Tests are set up so you can run them against the in-tree sources (no install needed). If you run into import issues, ensure `src/pyTOPP/src` is on `PYTHONPATH`, e.g.:
>
> ```bash
> export PYTHONPATH="$(pwd)/src/pyTOPP/src${PYTHONPATH:+:$PYTHONPATH}"
> ```

---

## Running via CTest

If you're building with CMake, the `src/pyTOPP/CMakeLists.txt` registers:

* one aggregate pytest run, and
* per-file pytest runs.

From your build directory:

```bash
ctest -R pytopp          # run only pyTOPP tests
ctest -VV -R pytopp      # verbose
```

---

## Snapshot tests (pytest-regtest)

Some tests use the `regtest` fixture to snapshot command lines for stability. If you deliberately change expected command formatting, you'll need to update the snapshots:

* Re-run the tests and follow the `pytest-regtest` plugin's instructions to refresh snapshots `--regtest-reset` (or delete the corresponding `.regtest` files and re-run).
  See the plugin's docs for the exact update flag supported by your version.

---

## Writing a new test

1. **Create a new file** named `test_*.py` in this folder.
   Example: `test_my_tool.py`.

2. **Import from the package under test**, not via relative filesystem paths, e.g.:

   ```python
   from pyopenms.pytopp import util
   from pyopenms.pytopp.ctdsupport import CTDHelpPrinter, CTDHelpConfig
   ```

3. **Use pytest fixtures**:

   * `tmp_path` for temporary files/directories
   * `monkeypatch` to alter env/argv, etc.
   * `capsys` to capture stdout/stderr
   * `regtest` for snapshotting stable text output

4. **CLI-style tests**: To test a CLI module that reads `sys.argv`, set it in the test:

   ```python
   import sys, importlib
   mod = importlib.import_module("pyopenms.pytopp.tools.pytopp_mytool_cli")
   monkeypatch.setattr(sys, "argv", ["pytopp_mytool_cli.py", "-in", "…", "-out", "…"])
   rc = mod.main()  # or catch SystemExit
   ```

5. **Keep inputs tiny**. If you need real data files (e.g., `.osw`), put them in the data folder with a small size and make tests skip if missing.

---

## Dependencies

Install test dependencies with your preferred tool:

```bash
# Using uv (fast, recommended)
uv pip install pytest pytest-regtest pytest-xdist

# Using pip
pip install pytest pytest-regtest pytest-xdist
```
