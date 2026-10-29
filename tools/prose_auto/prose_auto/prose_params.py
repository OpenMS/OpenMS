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
PRECURSOR_FLOOR_PPM = 10.0
# fragment floor per analyzer, from the benchmark sweep: Orbitrap HCD 20 ppm, Astral 10 ppm (frag 10 ppm was the
# Astral optimum, +9%), TOF / timsTOF 15 ppm (10 ppm lost, 40 ppm gained on tims)
FRAGMENT_FLOOR_PPM = {"orbitrap": 20.0, "astral": 10.0, "tof": 15.0}
ION_TRAP_FRAGMENT_DA = 0.35
CALIBRATION_OFFSET_PPM = 3.0  # calibrate when the measured precursor offset is larger (or unknown)



def _analyzer_class(v):
    """orbitrap | astral | tof, for the fragment-window floor"""
    ins, ms2 = str(v("instrument") or "").lower(), str(v("ms2_analyzer") or "").upper()
    if "astral" in ins:
        return "astral"
    if "TOF" in ms2 or any(k in ins for k in ("tims", "tof", "zenotof", "synapt")):
        return "tof"
    return "orbitrap"


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
    # The decisions hold the measured error (|offset| + 4 SD, or precision + offset); as a search window it is too
    # tight. On the benchmark groups (PR build, entrapment FDP) the 2x window lost 1-4% against 20 ppm on Orbitrap
    # HCD and gained 10% on ion-trap MS2 at 0.35 Da: windows are 2x the measured error with floors of 10 ppm
    # (precursor) and 20 ppm (fragment); ion-trap fragments 0.35 Da.
    lowres2 = v("ms2_analyzer") in ("IT", "LIT") or (ft and ft["unit"] == "Da")
    if pt:
        p = round(max(TOLERANCE_FACTOR * float(pt["value"]), PRECURSOR_FLOOR_PPM), 3) if pt["unit"] == "ppm" \
            else float(pt["value"])
        P["Search:precursor:mass_tolerance_lower"] = p  # ProSE declares tolerances as double
        P["Search:precursor:mass_tolerance_upper"] = p
        P["Search:precursor:mass_tolerance_unit"] = pt["unit"]
    if lowres2:
        P["Search:fragment:mass_tolerance"], P["Search:fragment:mass_tolerance_unit"] = ION_TRAP_FRAGMENT_DA, "Da"
    elif ft:
        P["Search:fragment:mass_tolerance"] = round(max(TOLERANCE_FACTOR * float(ft["value"]),
                                                       FRAGMENT_FLOOR_PPM[_analyzer_class(v)]), 3)
        P["Search:fragment:mass_tolerance_unit"] = ft["unit"]
    notes.append(f"search windows: {TOLERANCE_FACTOR:g}x the measured error, at least {PRECURSOR_FLOOR_PPM:g} ppm "
                 f"precursor and {FRAGMENT_FLOOR_PPM[_analyzer_class(v)]:g} ppm fragment on this analyzer "
                 f"({ION_TRAP_FRAGMENT_DA} Da on ion-trap MS2)")
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
        enz = ENZYME_NAMES.get(enz, enz)
        # a tryptic digest means Trypsin/P (no proline rule); it also found 2-6% more PSMs than Trypsin on 6 of 7
        # benchmark groups. Plain Trypsin is left to tune.py as a candidate.
        if enz == "Trypsin":
            enz = "Trypsin/P"
        P["Search:enzyme"] = enz
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
    # calibration pass when the offset is unknown or large, or the MS2 is low-res (velos: +7%); with a small measured
    # offset it cost 1-3% (hfx, lumos_lfq)
    po = decisions.get("precursor_offset")
    off = abs(float(po.value["value"])) if po and isinstance(po.value, dict) and po.value.get("unit") == "ppm" else None
    P["Search:calibration:enabled"] = "true" if (off is None or off > CALIBRATION_OFFSET_PPM or lowres2) else "false"
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
