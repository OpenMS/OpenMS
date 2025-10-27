# pyTOPP — Python TOPP tools for OpenMS

`pyTOPP` is a small, Python package shipped with OpenMS that
hosts Python TOPP wrappers (e.g., the `OpenSwathProphet` CLI). It's designed to:

- **Bundle** Python CLI tools alongside OpenMS TOPP tools 
- **Wrap** each tool with an OS-appropriate launcher in `bin/`
- **Bootstrap** the `pyopenms` namespace so `pyopenms.pytopp` is importable
- **Make development easy**: run tools directly from source, or via the wrapper

---

## Layout (w.r.t. OpenMS source tree)
```bash
src/pyTOPP/
├─ CMakeLists.txt                 # build + packaging logic
└─ src/pyopenms/pytopp/           # Python package (namespace under pyopenms)
   ├─ tools/
   │  ├─ pytopp_pyprophet_cli.py  # example tool
   │  └─ ...                      # add your tools here
   ├─ common.py
   ├─ ctdsupport.py
   ├─ util.py
   └─ __init__.py
```

At build time, everything under `src/pyopenms/pytopp` is copied into
`${build}/share/pytopp/site/pyopenms/pytopp`, and a tiny `sitecustomize.py` is
generated in `${build}/share/pytopp/site/` to extend the `pyopenms` namespace so
our vendored `pytopp` subtree is always discovered.

Wrappers are generated in:
```bash
${build}/bin/<ToolName>         # Linux/macOS
${build}/bin/<ToolName>.cmd     # Windows
```

Each wrapper prepends the vendored site path to `PYTHONPATH` and then runs the
tool (prefers `uv`, falls back to `python`).

---

## Prerequisites

- A standard OpenMS build (CMake + compiler)
- Optional but recommended for Python CLIs: [uv](https://docs.astral.sh/uv/)
- Python 3.9+ at runtime when not using `uv`

---

## Building with CMake

`src/CMakeLists.txt` contains:
```cmake
option(PYTOPP "Install pyTOPP CLI wrappers (Python tools via uv)" OFF)
add_subdirectory(pyTOPP)
```

If `PYTOPP` is `ON` (`-DPYTOPP=On`), building OpenMS will:

1. Copy `src/pyTOPP/src/pyopenms/pytopp` into the build tree
2. Generate `sitecustomize.py` in `share/pytopp/site`
3. Create per-tool wrappers in `${build}/bin`
4. Install both the wrappers and vendored Python package under `share/pytopp/site`

Reconfigure/rebuild after adding tools or changing the pyTOPP CMake.

---

## Using the tools

From the build tree:
```bash
# Wrapper (preferred)
./bin/OpenSwathProphet --help

# Or run the script directly via uv (shebang-friendly)
./share/pytopp/site/pyopenms/pytopp/tools/pytopp_pyprophet_cli.py --help
```

From source (development):
```bash
# Point Python at the source tree to resolve pyopenms.pytopp
export PYTHONPATH="$(pwd)/src/pyTOPP/src${PYTHONPATH:+:${PYTHONPATH}}"

# Now you can run the CLI straight from the repo:
./src/pyTOPP/src/pyopenms/pytopp/tools/pytopp_pyprophet_cli.py --help
```

> **Tip (Windows):** Use `set PYTHONPATH=%CD%\src\pyTOPP\src;%PYTHONPATH%`

---

## Creating a new pyTOPP tool

### 1) Add your CLI script

Create a new file under:
```
src/pyTOPP/src/pyopenms/pytopp/tools/pytopp_<name>_cli.py
```

Give it a modern uv shebang and dependencies header (so it's self-contained):
```python
#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "pyopenms>=3",            # if you import pyopenms
#   "CTDopts @ git+https://github.com/WorkflowConversion/CTDopts.git",
#   # ... your extra deps
# ]
# ///

from __future__ import annotations
import sys

from CTDopts.CTDopts import CTDModel
from pyopenms.pytopp.ctdsupport import (
    CTDHelpConfig,
    CTDHelpPrinter,
    PyTOPPArgParser,
)

def build_model() -> CTDModel:
    model = CTDModel(
        name="MyNewTool",
        description="Great example tool",
        category="Utilities",
    )
    model.add("in", required=True,  type="input-file", description="Input file")
    model.add("out", required=False, type="output-file", description="Output file")
    # add more parameters...
    return model

def main() -> int:
    model  = build_model()

    cfg = CTDHelpConfig(
        tool_name=f"{model.name} v{model.version} (pyTOPP)",
        binary_name="mytool",
        show_description=True,
        advanced_by_name={
            "section1:param4", "section2:param3", "section2:param4"
        },
    )
    cli_opt_printer = CTDHelpPrinter(model, cfg=cfg)

    argv = sys.argv[1:]
    if "-h" in argv or "--help" in argv:
        cli_opt_printer.print(advanced=False); return 0
    if "--help-advanced" in argv or "--helphelp" in argv:
        cli_opt_printer.print(advanced=True); return 0

    parser = OpenMSArgParser(model=model, mode="pure")   # or mode="pyopenms" for algo tools
    args   = parser.parse(argv)

    # use args[...] (flat colon-keys) and go do the work...
    print("Effective args:", args)
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

**Notes**

* Place your tool under `pyopenms/pytopp/tools/` so it's vendored and wrapped.
* Keep the name predictable, e.g. `pytopp_<name>_cli.py`.
* The `OpenMSArgParser` supports both:
  * `mode="pure"` — for non-pyOpenMS algorithm tools (merges defaults ← ini ← CLI)
  * `mode="openms"` — for classic pyOpenMS algo tools expecting OpenMS `Param` plumbing.

### 2) Register the tool in CMake

Open `src/pyTOPP/CMakeLists.txt` and add one line in the **"Register tools"** section:
```cmake
add_pytopp_tool(
  NAME    OpenSwathProphet            # wrapper name in bin/ alongside OpenMS TOPP tools
  CLI_REL pyopenms/pytopp/tools/pytopp_pyprophet_cli.py
)

# add yours:
add_pytopp_tool(
  NAME    MyNewTool
  CLI_REL pyopenms/pytopp/tools/pytopp_mynewtool_cli.py
)
```

Reconfigure & build:
```bash
cmake -G Ninja -DPYTOPP=On /home/you/source/OpenMS .
ninja -j<num_threads>
```

You'll get:
```
./bin/MyNewTool
./share/pytopp/site/pyopenms/pytopp/tools/pytopp_mynewtool_cli.py
```

### 3) Run it
```bash
./bin/MyNewTool --help
# or (dev)
export PYTHONPATH="$(pwd)/src/pyTOPP/src:${PYTHONPATH}"
./src/pyTOPP/src/pyopenms/pytopp/tools/pytopp_mynewtool_cli.py --help
```

---

## Argument parsing & INI merging

The provided `OpenMSArgParser` (in `pyopenms.pytopp.common`) supports:

* **Unknown flag detection** (fail early with hints)
* **CTD/INI merge semantics**:
```
  effective = defaults  <-  values from -ini  <-  explicit CLI flags
```
* **Forwarding "extra" args** using `+--...` tokens (they are unwrapped to be passed to the executable)
* **Boolean flags** with `-flag true/false`, `-flag=value`

You can choose:

* `mode="pure"` for tools that don't construct OpenMS `Param` objects
* `mode="pyopenms"` for tools that do (classic pyOpenMS algo path)
