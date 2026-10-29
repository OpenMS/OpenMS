"""Map decisions to ProSE parameters (the 'Search:' section of the ProSE TOPP
tool, see src/openms/source/ANALYSIS/ID/ProSEAlgorithm.cpp)."""

from __future__ import annotations

LABEL_MODS = {
    "TMT6plex": ["TMT6plex (K)", "TMT6plex (N-term)"],
    "TMT10/11plex": ["TMT6plex (K)", "TMT6plex (N-term)"],
    "TMT6/10/11plex": ["TMT6plex (K)", "TMT6plex (N-term)"],
    "TMT2plex": ["TMT2plex (K)", "TMT2plex (N-term)"],
    "TMT0 / TMT (126 only)": ["TMT (K)", "TMT (N-term)"],
    "TMTpro (16/18plex)": ["TMTpro (K)", "TMTpro (N-term)"],
    "iTRAQ4plex": ["iTRAQ4plex (K)", "iTRAQ4plex (N-term)"],
    "iTRAQ8plex": ["iTRAQ8plex (K)", "iTRAQ8plex (N-term)"],
}
SILAC_VARIABLE = {
    "SILAC (Lys8/Arg10)": ["Label:13C(6)15N(2) (K)", "Label:13C(6)15N(4) (R)"],
    "SILAC (Arg6)": ["Label:13C(6) (R)"],
    "SILAC Lys4 or dimethyl (d4)": ["Label:2H(4) (K)"],
}
ENRICHMENT_MODS = {
    "Phospho (STY)": ["Phospho (S)", "Phospho (T)", "Phospho (Y)"],
    "Phospho (Y)": ["Phospho (Y)"],
    "Acetyl (K)": ["Acetyl (K)"],
    "GlyGly (K)": ["GG (K)"],
}
ENZYME_NAMES = {"Trypsin": "Trypsin", "Trypsin/P": "Trypsin/P", "Lys-C": "Lys-C", "Arg-C": "Arg-C",
                "Chymotrypsin": "Chymotrypsin", "Glu-C": "glutamyl endopeptidase", "Asp-N": "Asp-N",
                "unspecific cleavage": "unspecific cleavage"}


TOLERANCE_FACTOR = 2.0  # search window / measured error (ppm tolerances)


def _per_residue(mods):
    """'Phospho (STY)' / 'Phospho (S,T,Y)' -> 'Phospho (S)', 'Phospho (T)', 'Phospho (Y)': ProSE takes one
    UniMod site per entry. Terminal sites ('N-term', 'Protein N-term') pass unchanged."""
    out = []
    for m in mods:
        name, _, site = m.rpartition(" (")
        res = site.rstrip(")").replace(",", "").replace("/", "").replace(" ", "")
        if name and res.isalpha() and res.isupper() and len(res) > 1:
            out += [f"{name} ({r})" for r in res if f"{name} ({r})" not in out]
        elif m not in out:
            out.append(m)
    return out


def to_prose(decisions: dict) -> dict:
    """Returns {'params': {name: value}, 'notes': [...]} with 'Search:'-prefixed
    ProSE parameter names."""
    v = lambda p: decisions[p].value if p in decisions else None
    P, notes = {}, []
    pt, ft = v("precursor_tolerance"), v("fragment_tolerance")
    # The decisions hold the measured error (|offset| + 4 SD, or precision + offset). As a search window that is
    # too tight: on 13 PRIDE runs, 2x the measured ppm value plus ProSE's calibration pass found up to 10% more
    # PSMs at the same entrapment FDP. Da windows (ion-trap MS2) are wide already and stay as measured.
    win = lambda t: round(float(t["value"]) * (TOLERANCE_FACTOR if t["unit"] == "ppm" else 1.0), 3)
    if pt:
        P["Search:precursor:mass_tolerance_lower"] = win(pt)  # ProSE declares tolerances as double
        P["Search:precursor:mass_tolerance_upper"] = win(pt)
        P["Search:precursor:mass_tolerance_unit"] = pt["unit"]
    if ft:
        P["Search:fragment:mass_tolerance"] = win(ft)
        P["Search:fragment:mass_tolerance_unit"] = ft["unit"]
    if (pt and pt["unit"] == "ppm") or (ft and ft["unit"] == "ppm"):
        notes.append(f"ppm tolerances are {TOLERANCE_FACTOR:g}x the measured error, with ProSE's calibration pass on")
    cr = v("charge_range")
    if cr:
        P["Search:precursor:min_charge"], P["Search:precursor:max_charge"] = int(cr[0]), int(cr[1])
    fixed, var = [], ["Oxidation (M)"]
    alk = v("cys_alkylation")
    if alk and not str(alk).startswith(("none", "unknown")):
        fixed.append(alk)
    lab = v("label")
    if lab in LABEL_MODS:
        fixed += LABEL_MODS[lab]
    elif lab in SILAC_VARIABLE:
        var += SILAC_VARIABLE[lab]
        notes.append("SILAC heavy labels are searched as variable modifications")
    enr = v("enrichment")
    if enr in ENRICHMENT_MODS:
        var += ENRICHMENT_MODS[enr]
    elif enr and str(enr).startswith("Glyco"):
        notes.append("glycopeptides detected: ProSE has no glycan search; consider a glyco engine")
    for m in v("variable_mods") or []:
        if isinstance(m, str) and m not in var and m not in fixed and "(" in m and decisions["variable_mods"].confidence >= 0.5:
            var.append(m)
    P["Search:modifications:fixed"] = _per_residue(fixed)
    P["Search:modifications:variable"] = _per_residue(var)
    enz = v("enzyme")
    if enz:
        P["Search:enzyme"] = ENZYME_NAMES.get(enz, enz)
    ic = v("instrument_class")
    if ic and v("ms2_analyzer") not in ("IT", "LIT"):   # PeptDeep models high-res HCD-type spectra
        # ProSE accepts only these PeptDeep classes; free text from records ('TIMS quadrupole time-of-flight') is mapped
        ic = "timsTOF" if "tims" in str(ic).lower() else ic
        if ic in ("QE", "Lumos", "timsTOF", "SciexTOF", "ThermoTOF"):
            P["Search:peptdeep:instrument"] = ic
    ce = v("collision_energy")
    if ce and v("fragmentation") in ("HCD", "CID"):
        P["Search:peptdeep:nce"] = -1.0   # ProSE picks it from the spectra and refines it
    if v("fragmentation") and any(x in str(v("fragmentation")) for x in ("ETD", "ECD", "EThcD")):
        P["Search:ions:by_activation"] = "true"
    # calibration pass always: it re-centres and narrows the window on ProSE's own PSMs, which the mini-search's
    # offset (2,000 spectra, no variable mods) cannot replace
    P["Search:calibration:enabled"] = "true"
    acq = str(v("acquisition") or "")
    if acq.startswith("DIA"):
        notes.append("DIA data: search pseudo-spectra (DIA-Umpire style deconvolution; OpenMS has "
                     "MasstraceCorrelator/ClusterMassTraces) or use a DIA engine; ProSE searches DDA-style spectra")
    if acq == "PRM":
        notes.append("PRM data: a targeted (Skyline-style) analysis usually fits better than a database search")
    return {"params": P, "notes": notes}


def to_cli(params: dict) -> list:
    args = []
    for k, val in params.items():
        args.append(f"-{k}")
        if isinstance(val, list):
            args.extend(str(x) for x in val)
        else:
            args.append(str(val))
    return args


def to_ini_snippet(params: dict) -> str:
    """ProSE parameters as a TOPP INI fragment (ITEM/ITEMLIST under ProSE:1)."""
    lines = ['<?xml version="1.0" encoding="ISO-8859-1"?>',
             '<PARAMETERS version="1.7.0">', '  <NODE name="ProSE">', '    <NODE name="1">']
    tree = {}
    for k, val in params.items():
        node = tree
        parts = k.split(":")
        for p in parts[:-1]:
            node = node.setdefault(p, {})
        node[parts[-1]] = val

    def emit(node, ind):
        for k, val in node.items():
            if isinstance(val, dict):
                lines.append(f'{ind}<NODE name="{k}">')
                emit(val, ind + "  ")
                lines.append(f"{ind}</NODE>")
            elif isinstance(val, list):
                lines.append(f'{ind}<ITEMLIST name="{k}" type="string">')
                for x in val:
                    lines.append(f'{ind}  <LISTITEM value="{x}"/>')
                lines.append(f"{ind}</ITEMLIST>")
            else:
                t = "double" if isinstance(val, float) else ("int" if isinstance(val, int) else "string")
                lines.append(f'{ind}<ITEM name="{k}" value="{val}" type="{t}"/>')
    emit(tree, "      ")
    lines += ["    </NODE>", "  </NODE>", "</PARAMETERS>"]
    return "\n".join(lines) + "\n"


def to_sdrf_technical(decisions: dict) -> dict:
    """The SDRF comment[] columns that describe the search (like techSDRF)."""
    v = lambda p: decisions[p].value if p in decisions else None
    out = {}
    if v("instrument"):
        out["comment[instrument]"] = f"NT={v('instrument')}"
    if v("fragmentation"):
        out["comment[dissociation method]"] = f"NT={v('fragmentation')}"
    if v("acquisition"):
        a = str(v("acquisition"))
        out["comment[proteomics data acquisition method]"] = \
            "NT=Data-independent acquisition" if a.startswith("DIA") else "NT=Data-dependent acquisition"
    if v("label"):
        out["comment[label]"] = "label free sample" if v("label") == "label free" else str(v("label"))
    if v("enzyme"):
        out["comment[cleavage agent details]"] = f"NT={v('enzyme')}"
    pt, ft = v("precursor_tolerance"), v("fragment_tolerance")
    if pt:
        out["comment[precursor mass tolerance]"] = f"{pt['value']} {pt['unit']}"
    if ft:
        out["comment[fragment mass tolerance]"] = f"{ft['value']} {ft['unit']}"
    alk = v("cys_alkylation")
    if alk and not str(alk).startswith(("none", "unknown")):
        name = alk.split(" (")[0]
        out["comment[modification parameters]"] = f"NT={name};TA=C;MT=Fixed"
    if v("organism") and v("organism") not in ("unknown", "not in panel"):
        out["characteristics[organism]"] = str(v("organism")).replace("_", " ")
    return out
