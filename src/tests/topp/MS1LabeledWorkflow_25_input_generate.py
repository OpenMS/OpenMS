#!/usr/bin/env python3
"""Write MS1LabeledWorkflow_25_input.mzML, MS1LabeledWorkflow_25.idXML, MS1LabeledWorkflow_25_heavy_id.idXML (SVEDLAFGNYQGEDVTALFEK
identified in its heavy channel) and MS1LabeledWorkflow_25.fasta (synthetic triple SILAC,
[][Lys4,Arg6][Lys8,Arg10], targeted quantification).

Centroided MS1 spectra, one per second from 100 to 660 s and at 700 and 701 s, 2+ peptides eluting as Gaussians, light identifications. Only the first two
isotopes of every channel are recorded (so the detection, which needs three, finds no multiplet); an isotope of another channel that
lies closer to a recorded peak than the FWHM at a resolution of 120,000 at m/z 200 adds to it (unresolved, one centroid at the
intensity-weighted m/z), others are not recorded. Peaks below 'FLOOR' are not recorded either.
  SVEDLAFGNYQGEDVTALFEK 2.7:1:1 and TALFEQLSVEDLAFGNYQGDK 11.6:1:1 (2.3 kDa): the medium M0 lies on the light M+4 (Lys4 - 4 x 13C
    = 0.0117 Da, 5 ppm), which makes up 30 % and 65 % of it (the light M+5 10 % and 32 % of the medium M+1).
  LVNELTEFAK 1:2:2, weak: the light M+1 is above the floor in only the two spectra at the apex (its M0 in six, the medium's and the
    heavy's M0 and M+1 in eight or more).
  YVNDLAGHEK 4:1:1, identified twice: in the weak medium channel (M+1 above the floor in two spectra only) 3 s before the apex,
    and in the light channel 2 s after it.
  SAGLDFVEEK 1:0:0, weak like LVNELTEFAK (light M+1 above the floor in two spectra only), nothing at the medium and heavy m/z.
  ADLTVEFGNK 1:0:6, the light weak as above, the heavy strong; identified in the light as ADLTVEFGNK (better score) and in the heavy
    as its isomer DALTVEFGNK, one signal.
  NEFLGSVDYK 1:1:1, eluting so narrowly (sigma 0.5 s, apex 580.3 s) that all its peaks are above the floor in the two spectra at 580
    and 581 s only, amid spectra every second, plus one low peak (at the floor, 4 % of its summed isotopes at 580 s) at the light M0
    at 582 s; identified in the light at 580.5 s.
  GSAEVTLDFK 1:1:1, in the two spectra at 700 and 701 s only (no other spectrum within 'algorithm:rt_typical', 40 s), apex at
    700.3 s, identified in the light at 700.5 s.
The isotope distributions are computed from the elemental formulas, as CoarseIsotopePatternGenerator does.

Usage: python3 MS1LabeledWorkflow_25_input_generate.py (writes the files next to this script)
"""
import base64
import functools
import math
import os
import struct

HERE = os.path.dirname(os.path.abspath(__file__))
PROTON = 1.007276466812
C13 = 1.0033548378
LYS4, LYS8 = 4.025107, 8.014199
FLOOR = 6.0e5
RESIDUES = {"G": "C2H3NO", "A": "C3H5NO", "S": "C3H5NO2", "P": "C5H7NO", "V": "C5H9NO", "T": "C4H7NO2", "L": "C6H11NO",
            "I": "C6H11NO", "N": "C4H6N2O2", "D": "C4H5NO3", "Q": "C5H8N2O2", "K": "C6H12N2O", "E": "C5H7NO3", "M": "C5H9NOS",
            "H": "C6H7N3O", "F": "C9H9NO", "R": "C6H12N4O", "Y": "C9H9NO2", "W": "C11H10N2O"}
MONO = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956, "S": 31.97207100}
ABUNDANCE = {"C": [0.9893, 0.0107], "H": [0.999885, 0.000115], "N": [0.99636, 0.00364], "O": [0.99757, 0.00038, 0.00205],
             "S": [0.9493, 0.0076, 0.0429, 0.0, 0.0002]}
# (sequence, light : medium (or the share of the light M+4 in the medium M0 it is chosen for), apex RT, sigma, medium M0 at the apex)
PEPTIDES = [("SVEDLAFGNYQGEDVTALFEK", "A", 130.0, 5.0, 1.0e8),
            ("TALFEQLSVEDLAFGNYQGDK", "B", 230.0, 5.0, 1.0e8),
            ("LVNELTEFAK", 0.5, 330.5, 3.0, 2.0e6),
            ("YVNDLAGHEK", 4.0, 430.5, 3.0, 1.1e6),
            ("SAGLDFVEEK", 1.0, 530.5, 3.0, None),
            ("ADLTVEFGNK", 1.0, 630.5, 3.0, None),
            ("NEFLGSVDYK", 1.0, 580.3, 0.5, 1.0e7),
            ("GSAEVTLDFK", 1.0, 700.3, 3.0, 1.0e7)]
# medium and heavy amounts where they are not 1 (relative to the light : medium unit)
AMOUNTS = {"SAGLDFVEEK": (0.0, 0.0), "ADLTVEFGNK": (0.0, 6.0)}
# single low peaks at a peptide's light M0, no envelope: (RT, sequence, intensity)
STRAY_PEAKS = [(582.0, "NEFLGSVDYK", FLOOR)]
SPECTRUM_RTS = list(range(100, 661)) + [700, 701]
# identifications: (sequence, RT, labelled channel shift, PEP)
IDS = [("SVEDLAFGNYQGEDVTALFEK", 130.0, 0.0, "1.0e-03"), ("TALFEQLSVEDLAFGNYQGDK", 230.0, 0.0, "1.0e-03"),
       ("LVNELTEFAK", 330.5, 0.0, "1.0e-03"), ("YVNDLAGHEK(Label:2H(4))", 427.5, LYS4, "1.0e-03"), ("YVNDLAGHEK", 432.5, 0.0, "1.0e-03"),
       ("SAGLDFVEEK", 530.5, 0.0, "1.0e-03"), ("ADLTVEFGNK", 630.5, 0.0, "1.0e-04"),
       ("DALTVEFGNK(Label:13C(6)15N(2))", 631.5, LYS8, "1.0e-03"), ("GSAEVTLDFK", 700.5, 0.0, "1.0e-03"), ("NEFLGSVDYK", 580.5, 0.0, "1.0e-03")]
# MS1LabeledWorkflow_25_heavy_id.idXML: the same with SVEDLAFGNYQGEDVTALFEK identified in its heavy (Lys8) channel
HEAVY_IDS = [("SVEDLAFGNYQGEDVTALFEK(Label:13C(6)15N(2))", rt, LYS8, pep) if seq == "SVEDLAFGNYQGEDVTALFEK" else (seq, rt, shift, pep)
             for seq, rt, shift, pep in IDS]
PROTEIN_SEQUENCES = list(dict.fromkeys(seq.split("(")[0] for seq, _, _, _ in IDS))
TARGET_FRACTION = {"A": 0.3, "B": 0.65}  # the light M+4's share of the medium M0


def composition(sequence):
    atoms = {"H": 2, "O": 1}
    for residue in sequence:
        formula, i = RESIDUES[residue], 0
        while i < len(formula):
            element = formula[i]
            j = i + 1
            while j < len(formula) and formula[j].isdigit():
                j += 1
            atoms[element] = atoms.get(element, 0) + (int(formula[i + 1:j]) if j > i + 1 else 1)
            i = j
    return atoms


@functools.lru_cache(maxsize=None)
def isotopes(sequence, n=12):
    dist = [1.0]
    for element, count in composition(sequence).items():
        for _ in range(count):
            new = [0.0] * min(n, len(dist) + len(ABUNDANCE[element]) - 1)
            for a, pa in enumerate(dist):
                for b, pb in enumerate(ABUNDANCE[element]):
                    if a + b < len(new):
                        new[a + b] += pa * pb
            dist = new
    return [d / dist[0] for d in dist] + [0.0] * (n - len(dist))


def mz(sequence):
    mass = sum(MONO[e] * c for e, c in composition(sequence).items())
    return (mass + 2 * PROTON) / 2


def fwhm(value):
    return value / (120000.0 * math.sqrt(200.0 / value))


def light_to_medium(sequence, ratio):
    if isinstance(ratio, float):
        return ratio
    i4 = isotopes(sequence)[4]
    f = TARGET_FRACTION[ratio]
    return round(f / ((1.0 - f) * i4), 1)


def encode(values, fmt):
    return base64.b64encode(struct.pack("<%d%s" % (len(values), fmt), *values)).decode()


def spectrum_xml(index, rt, peaks):
    mzs = [p[0] for p in peaks]
    its = [p[1] for p in peaks]
    return f"""			<spectrum id="scan={index + 1}" index="{index}" defaultArrayLength="{len(peaks)}">
				<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
				<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1" />
				<cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
				<scanList count="1">
					<cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
					<scan>
						<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt:g}" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
					</scan>
				</scanList>
				<binaryDataArrayList count="2">
					<binaryDataArray encodedLength="{len(encode(mzs, 'd'))}">
						<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
						<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
						<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
						<binary>{encode(mzs, 'd')}</binary>
					</binaryDataArray>
					<binaryDataArray encodedLength="{len(encode(its, 'f'))}">
						<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
						<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
						<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS" />
						<binary>{encode(its, 'f')}</binary>
					</binaryDataArray>
				</binaryDataArrayList>
			</spectrum>
"""


def peaks_at(rt):
    peaks = []
    for sequence, ratio, apex, sigma, height in PEPTIDES:
        medium, heavy = AMOUNTS.get(sequence, (1.0, 1.0))
        if height is None:  # light M+1 at 0.93 x the floor at the apex: above it in the two spectra next to the apex only
            height = FLOOR / (light_to_medium(sequence, ratio) * isotopes(sequence)[1] * 0.93)
        scale = height * math.exp(-((rt - apex) ** 2) / (2 * sigma ** 2))
        if scale * min(1.0, light_to_medium(sequence, ratio)) < FLOOR / 10:
            continue
        dist = isotopes(sequence)
        light = mz(sequence)
        r = light_to_medium(sequence, ratio)
        channels = [(light, r), (light + LYS4 / 2, medium), (light + LYS8 / 2, heavy)]
        recorded = [[c_mz + k * C13 / 2, scale * amount * dist[k]] for c_mz, amount in channels for k in range(2)]
        for c_mz, amount in channels:
            for k in range(2, len(dist)):
                position = c_mz + k * C13 / 2
                nearest = min(recorded, key=lambda p: abs(p[0] - position))
                added = scale * amount * dist[k]
                if added > 0.0 and abs(nearest[0] - position) < fwhm(position):
                    nearest[0] = (nearest[0] * nearest[1] + position * added) / (nearest[1] + added)
                    nearest[1] += added
        peaks += [(p[0], round(p[1])) for p in recorded if p[1] >= FLOOR]
    peaks += [(mz(sequence), round(intensity)) for stray_rt, sequence, intensity in STRAY_PEAKS if stray_rt == rt]
    return sorted(peaks)


def write_mzml(path):
    spectra = [spectrum_xml(index, rt, peaks_at(float(rt))) for index, rt in enumerate(SPECTRUM_RTS)]
    with open(path, "w") as out:
        out.write(f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
	<cvList count="2">
		<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"/>
		<cv id="UO" fullName="Unit Ontology" URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>
	</cvList>
	<fileDescription>
		<fileContent>
			<cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
		</fileContent>
	</fileDescription>
	<softwareList count="1">
		<software id="so_default" version="" >
			<cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="" />
		</software>
	</softwareList>
	<instrumentConfigurationList count="1">
		<instrumentConfiguration id="ic_0">
			<cvParam cvRef="MS" accession="MS:1000031" name="instrument model" />
		</instrumentConfiguration>
	</instrumentConfigurationList>
	<dataProcessingList count="1">
		<dataProcessing id="dp_sp_0">
			<processingMethod order="0" softwareRef="so_default">
				<cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" />
			</processingMethod>
		</dataProcessing>
	</dataProcessingList>
	<run id="ru_0" defaultInstrumentConfigurationRef="ic_0">
		<spectrumList count="{len(spectra)}" defaultDataProcessingRef="dp_sp_0">
{''.join(spectra)}		</spectrumList>
	</run>
</mzML>
""")


def write_idxml(path, ids):
    protein = {seq: i for i, seq in enumerate(PROTEIN_SEQUENCES)}
    entries = "".join(f"""		<PeptideIdentification score_type="Posterior Error Probability" higher_score_better="false" significance_threshold="0.0" MZ="{mz(seq.split('(')[0]) + shift / 2:.15f}" RT="{rt}" spectrum_reference="scan={SPECTRUM_RTS.index(int(rt)) + 1}" >
			<PeptideHit score="{pep}" sequence="{seq}" charge="2" protein_refs="PH_{protein[seq.split('(')[0]]}" >
				<UserParam type="string" name="target_decoy" value="target"/>
			</PeptideHit>
		</PeptideIdentification>
""" for seq, rt, shift, pep in ids)
    proteins = "".join(f"""			<ProteinHit id="PH_{i}" accession="PROT_{i}" score="0.0" sequence="" >
			</ProteinHit>
""" for i in range(len(PROTEIN_SEQUENCES)))
    with open(path, "w") as out:
        out.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<IdXML version="1.5" xsi:noNamespaceSchemaLocation="https://www.openms.de/xml-schema/IdXML_1_5.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
	<SearchParameters id="SP_0" db="MS1LabeledWorkflow_25.fasta" db_version="" taxonomy="" mass_type="monoisotopic" charges="" enzyme="trypsin" missed_cleavages="0" precursor_peak_tolerance="0" precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0" peak_mass_tolerance_ppm="false" >
		<VariableModification name="Label:2H(4) (K)" />
		<VariableModification name="Label:13C(6) (R)" />
		<VariableModification name="Label:13C(6)15N(2) (K)" />
		<VariableModification name="Label:13C(6)15N(4) (R)" />
	</SearchParameters>
	<IdentificationRun date="0000-00-00T00:00:00" search_engine="synthetic" search_engine_version="" search_parameters_ref="SP_0" >
		<ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0.0" >
{proteins}			<UserParam type="stringList" name="spectra_data" value="[MS1LabeledWorkflow_25_input.mzML]"/>
		</ProteinIdentification>
{entries}	</IdentificationRun>
</IdXML>
""")


def write_fasta(path):
    with open(path, "w") as out:
        for i, seq in enumerate(PROTEIN_SEQUENCES):
            out.write(f">PROT_{i}\nMK{seq}GGSGGSGGR\n")


if __name__ == "__main__":
    for seq, ratio, _, _, _ in PEPTIDES:
        dist = isotopes(seq)
        r = light_to_medium(seq, ratio)
        print(f"{seq}: m/z {mz(seq):.6f}, I1..I5 {[round(d, 4) for d in dist[1:6]]}, light:medium {r}, "
              f"light M+4 in the medium M0 {r * dist[4] / (r * dist[4] + 1):.3f}, light M+5 in the medium M+1 {r * dist[5] / (r * dist[5] + dist[1]):.3f}")
    write_mzml(os.path.join(HERE, "MS1LabeledWorkflow_25_input.mzML"))
    write_idxml(os.path.join(HERE, "MS1LabeledWorkflow_25.idXML"), IDS)
    write_idxml(os.path.join(HERE, "MS1LabeledWorkflow_25_heavy_id.idXML"), HEAVY_IDS)
    write_fasta(os.path.join(HERE, "MS1LabeledWorkflow_25.fasta"))
