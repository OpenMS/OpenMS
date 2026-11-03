#!/usr/bin/env python3
"""Regenerate OpenSwathWorkflow_im_runs_input.{mzML,tsv} from scratch.

A synthetic diaPASEF run of eight iRT kit peptides and their spectral library. The peptides elute
at RT = 2 * iRT + 100 s, so the library also serves as linear iRT file, and their ion mobility in
the data is IM_SHIFT above the library value, so the calibration estimates a non-trivial ion
mobility correction. Each cycle (every 3 s) holds one MS1 frame and one diaPASEF MS2 frame
(m/z 400-800, 1/K0 0.6-1.4); every peak carries its ion mobility.

Requires only the Python standard library (no OpenMS build).
"""
import base64
import math
import struct
import zlib
from pathlib import Path

OUT_DIR = Path(__file__).parent
OUT_MZML = OUT_DIR / "OpenSwathWorkflow_im_runs_input.mzML"
OUT_TSV = OUT_DIR / "OpenSwathWorkflow_im_runs_input.tsv"

IM_SHIFT = 0.03          # data ion mobility = library ion mobility + IM_SHIFT
RT_SLOPE, RT_INTERCEPT = 2.0, 100.0
RT_SIGMA = 6.0           # elution peak width (s)
IM_SIGMA = 0.01          # ion mobility peak width (1/K0)
IM_OFFSETS = (-0.02, -0.01, 0.0, 0.01, 0.02)
CYCLE_TIME = 3.0
RT_START, RT_END = 30.0, 330.0
CHARGE = 2

# (sequence, iRT, library ion mobility)
PEPTIDES = [
    ("LGGNEQVTR", -28.3, 0.80),
    ("GAGSSEPVTGLDAK", 0.2, 0.95),
    ("VEATFGVDESNAK", 12.4, 1.00),
    ("YILAGVENSK", 19.8, 0.85),
    ("TPVISGGPYEYR", 28.7, 1.05),
    ("DGLDAASYYAPVR", 43.0, 1.10),
    ("GTFIIDPGGVIR", 71.3, 0.90),
    ("LFLQFGAQGSPFLK", 98.1, 1.15),
]
FRAGMENT_INTENSITIES = (100.0, 80.0, 60.0, 45.0, 30.0, 20.0)  # y3 .. y8 (as far as the peptide allows)

MONO = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768,
        "L": 113.08406, "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
        "E": 129.04259, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931}
WATER, PROTON, C13 = 18.010565, 1.007276, 1.0033548


def precursor_mz(seq):
    return (sum(MONO[a] for a in seq) + WATER + CHARGE * PROTON) / CHARGE


def y_ions(seq):
    """(series number, m/z) of the singly charged y3 .. y8 ions"""
    ions = []
    for n in range(3, min(8, len(seq) - 1) + 1):
        ions.append((n, sum(MONO[a] for a in seq[-n:]) + WATER + PROTON))
    return ions


def encode(values, fmt):
    return base64.b64encode(zlib.compress(struct.pack("<%d%s" % (len(values), fmt), *values))).decode()


def data_array(values, fmt, name_param):
    width = '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>' if fmt == "d" else \
            '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"/>'
    binary = encode(values, fmt)
    return ('<binaryDataArray encodedLength="%d">%s%s<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression"/>'
            '<binary>%s</binary></binaryDataArray>' % (len(binary), name_param, width, binary))


def spectrum(index, rt, ms_level, peaks):
    peaks = sorted(peaks)
    mz = [p[0] for p in peaks]
    intensity = [p[1] for p in peaks]
    im = [p[2] for p in peaks]
    out = ['<spectrum id="spectrum=%d" index="%d" defaultArrayLength="%d">' % (index, index, len(peaks)),
           '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="%d" />' % ms_level,
           '<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />']
    if ms_level == 2:
        out.append('<userParam name="ion mobility lower limit" type="xsd:double" value="0.6"/>'
                   '<userParam name="ion mobility upper limit" type="xsd:double" value="1.4"/>')
    out.append('<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795" name="no combination" /><scan>'
               '<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="%g" unitAccession="UO:0000010" '
               'unitName="second" unitCvRef="UO" /></scan></scanList>' % rt)
    if ms_level == 2:
        out.append('<precursorList count="1"><precursor><isolationWindow>'
                   '<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="600" />'
                   '<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="200" />'
                   '<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="200" />'
                   '</isolationWindow><activation><cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />'
                   '</activation></precursor></precursorList>')
    out.append('<binaryDataArrayList count="3">')
    out.append(data_array(mz, "d", '<cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>'))
    out.append(data_array(intensity, "f", '<cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>'))
    out.append(data_array(im, "f", '<cvParam cvRef="MS" accession="MS:1000786" name="non-standard data array" value="Ion Mobility"/>'))
    out.append('</binaryDataArrayList></spectrum>')
    return "".join(out)


def peaks_at(rt, ions):
    """ions: (m/z, relative intensity, ion mobility) -> peaks (m/z, intensity, ion mobility) eluting at rt"""
    peaks = []
    for mz, rel, im, rt0 in ions:
        elution = math.exp(-0.5 * ((rt - rt0) / RT_SIGMA) ** 2)
        if elution < 0.01:
            continue
        for offset in IM_OFFSETS:
            intensity = 1e4 * rel * elution * math.exp(-0.5 * (offset / IM_SIGMA) ** 2)
            peaks.append((mz, intensity, im + offset))
    return peaks


def main():
    ms1_ions, ms2_ions = [], []
    rows = []
    header = ["PrecursorMz", "ProductMz", "PeptideSequence", "ModifiedPeptideSequence", "ProteinId", "FragmentType",
              "FragmentSeriesNumber", "PrecursorIonMobility", "TransitionGroupId", "NormalizedRetentionTime",
              "TransitionId", "LibraryIntensity", "PrecursorCharge", "ProductCharge", "Decoy"]
    for seq, irt, im in PEPTIDES:
        rt0 = RT_SLOPE * irt + RT_INTERCEPT
        prec = precursor_mz(seq)
        for iso, rel in enumerate((100.0, 60.0, 25.0)):
            ms1_ions.append((prec + iso * C13 / CHARGE, rel, im + IM_SHIFT, rt0))
        for (n, mz), rel in zip(y_ions(seq), FRAGMENT_INTENSITIES):
            ms2_ions.append((mz, rel, im + IM_SHIFT, rt0))
            rows.append(["%.5f" % prec, "%.5f" % mz, seq, seq, "iRT_protein", "y", str(n), "%g" % im, seq + "_2",
                         "%g" % irt, "%s_y%d" % (seq, n), "%g" % rel, str(CHARGE), "1", "0"])

    with open(OUT_TSV, "w") as tsv:
        tsv.write("\t".join(header) + "\n")
        for row in rows:
            tsv.write("\t".join(row) + "\n")

    spectra = []
    rt = RT_START
    while rt <= RT_END:
        spectra.append(spectrum(len(spectra), rt, 1, peaks_at(rt, ms1_ions)))
        spectra.append(spectrum(len(spectra), rt + 0.1, 2, peaks_at(rt + 0.1, ms2_ions)))
        rt += CYCLE_TIME

    with open(OUT_MZML, "w") as mzml:
        mzml.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n'
                   '<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                   'xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" version="1.1.0">\n'
                   '<cvList count="2">'
                   '<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"/>'
                   '<cv id="UO" fullName="Unit Ontology" URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>'
                   '</cvList>\n'
                   '<fileDescription><fileContent><cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" />'
                   '<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" /></fileContent></fileDescription>\n'
                   '<softwareList count="1"><software id="so_default" version="">'
                   '<cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="OpenSwathWorkflow_im_runs_input_generate.py" />'
                   '</software></softwareList>\n'
                   '<instrumentConfigurationList count="1"><instrumentConfiguration id="ic_0">'
                   '<cvParam cvRef="MS" accession="MS:1000031" name="instrument model" /></instrumentConfiguration></instrumentConfigurationList>\n'
                   '<dataProcessingList count="1"><dataProcessing id="dp_0"><processingMethod order="0" softwareRef="so_default">'
                   '<cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" /></processingMethod></dataProcessing></dataProcessingList>\n'
                   '<run id="synthetic_diaPASEF" defaultInstrumentConfigurationRef="ic_0">\n'
                   '<spectrumList count="%d" defaultDataProcessingRef="dp_0">\n' % len(spectra))
        for s in spectra:
            mzml.write(s + "\n")
        mzml.write('</spectrumList>\n</run>\n</mzML>\n')


if __name__ == "__main__":
    main()
