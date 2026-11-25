"""
matplotlib is an optional dependency of pyOpenMS (the ``plotting`` extra, OpenMS
issue #6728). ``import pyopenms`` and ``import pyopenms.plotting`` must work
without it; a plot function called without it must name the extra to install.
"""

import subprocess
import sys

import pytest

import pyopenms
from pyopenms.plotting import mirror_plot_spectrum, plot_chromatogram, plot_spectrum


def test_import_does_not_load_matplotlib():
    # A new interpreter: other tests of this session import matplotlib themselves.
    code = (
        "import sys, pyopenms, pyopenms.plotting\n"
        "assert 'matplotlib' not in sys.modules, 'importing pyopenms loaded matplotlib'\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    "plot",
    [
        lambda: plot_chromatogram(pyopenms.MSChromatogram()),
        lambda: plot_spectrum(pyopenms.MSSpectrum()),
        lambda: mirror_plot_spectrum(pyopenms.MSSpectrum(), pyopenms.MSSpectrum()),
    ],
    ids=["plot_chromatogram", "plot_spectrum", "mirror_plot_spectrum"],
)
def test_plot_without_matplotlib_names_the_extra(monkeypatch, plot):
    # None in sys.modules makes the import fail as if matplotlib were not installed.
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", None)
    with pytest.raises(ImportError, match=r"pyopenms\[plotting\]"):
        plot()
