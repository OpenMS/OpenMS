"""Evidence and decisions: the common currency of all detectors.

Every detector (raw data, file metadata, repository records, free text, LLM)
emits Evidence objects. The reconciler turns the evidence for one parameter
into a Decision with a confidence and a list of conflicts. The GUI shows the
Decision and lets the user drill into the Evidence behind it.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any

# Where a piece of evidence comes from. The order is the default trust order
# for parameters that the data can measure (see reconcile.py).
SOURCES = (
    "data",        # measured from the spectra (reporter ions, mass errors, ...)
    "search",      # measured from a quick ProSE pre-search on a subset
    "metadata",    # read from the file (instrument model, activation, windows)
    "repository",  # PRIDE / MassIVE / SDRF records of the dataset
    "user",        # the user's own description of the experiment
    "text",        # methods section or repository protocols (rule-based)
    "llm",         # an LLM reading user description / methods / literature
    "default",     # instrument-class or tool default
)

# Parameters the tool knows about. The keys are stable: they are the contract
# with the GUI (report.json) and with prose_params.py.
PARAMETERS = {
    "acquisition": "Acquisition mode (DDA, DIA, PRM)",
    "dia_windows": "DIA isolation window scheme",
    "instrument": "Instrument model",
    "instrument_class": "Instrument class (PeptDeep naming)",
    "ms1_analyzer": "MS1 mass analyzer",
    "ms2_analyzer": "MS2 mass analyzer",
    "fragmentation": "Fragmentation / activation method",
    "collision_energy": "Collision energy (NCE)",
    "ion_mobility": "Ion mobility separation",
    "precursor_tolerance": "Precursor mass tolerance",
    "precursor_offset": "Systematic precursor mass offset",
    "fragment_tolerance": "Fragment mass tolerance",
    "fragment_offset": "Systematic fragment mass offset",
    "charge_range": "Precursor charge range",
    "label": "Quantification label",
    "cys_alkylation": "Cysteine alkylation (fixed modification)",
    "enrichment": "PTM enrichment",
    "variable_mods": "Variable modifications",
    "enzyme": "Protease",
    "organism": "Organism (database)",
}


@dataclass
class Evidence:
    parameter: str
    value: Any
    confidence: float          # 0..1, the detector's own estimate
    source: str                # one of SOURCES
    method: str                # short name of the detector, e.g. "reporter-ions"
    detail: str = ""           # human readable explanation
    stats: dict = field(default_factory=dict)  # numbers behind it, for the GUI
    quote: str = ""            # verbatim text span for text/llm evidence

    def __post_init__(self):
        if self.parameter not in PARAMETERS:
            raise ValueError(f"unknown parameter {self.parameter!r}")
        if self.source not in SOURCES:
            raise ValueError(f"unknown source {self.source!r}")
        self.confidence = float(max(0.0, min(1.0, self.confidence)))

    def to_dict(self):
        return asdict(self)


@dataclass
class Decision:
    parameter: str
    value: Any
    confidence: float
    source: str                 # source of the winning evidence
    evidence: list = field(default_factory=list)   # all Evidence for the parameter
    conflicts: list = field(default_factory=list)  # human readable conflict notes
    status: str = "auto"        # auto | review | default | unknown

    def to_dict(self):
        d = asdict(self)
        d["label"] = PARAMETERS[self.parameter]
        return d
