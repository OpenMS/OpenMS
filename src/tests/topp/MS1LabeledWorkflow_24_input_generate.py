#!/usr/bin/env python3
"""Write MS1LabeledWorkflow_24_input.mzML and MS1LabeledWorkflow_24.idXML (synthetic triple SILAC, [][Lys4,Arg6][Lys8,Arg10]).

Centroided MS1 spectra, one per second from 100 to 160 s, 2+ peptides eluting as a Gaussian (sigma 5 s, apex 130 s, 1e6):
  YVNDLAGHEK 1:1:1, an ordinary triplet (four isotopes per channel in every spectrum), identified in the light channel.
  LVNELTEFAK 1:1:1, whose channels the detection finds separately: in every third spectrum one channel shows four isotopes and the
    other two only their first two, in turn, so that no spectrum holds two complete channels. The light and the medium channel are
    identified, the heavy one is not.

Usage: python3 MS1LabeledWorkflow_24_input_generate.py (writes the two files next to this script)
"""
import base64
import math
import os
import struct

HERE = os.path.dirname(os.path.abspath(__file__))
SPACING = 1.0033548 / 2  # 13C-12C mass difference at charge 2
LYS4, LYS8 = 4.025107, 8.014199
LIGHT_PATTERN = [1.0, 0.6077, 0.2166, 0.0568]
HEAVY_PATTERN = [1.0, 0.5617, 0.1898, 0.0476]
PEPTIDES = {"YVNDLAGHEK": 573.282923178970805, "LVNELTEFAK": 582.318974338470981}


def channels(mz):
    return [(mz, LIGHT_PATTERN), (mz + LYS4 / 2, LIGHT_PATTERN), (mz + LYS8 / 2, HEAVY_PATTERN)]


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


def write_mzml(path):
    spectra = []
    for index, rt in enumerate(range(100, 161)):
        scale = 1e6 * math.exp(-((rt - 130) ** 2) / (2 * 5.0 ** 2))
        peaks = []
        for name, mz in PEPTIDES.items():
            for c, (channel_mz, pattern) in enumerate(channels(mz)):
                isotopes = 4 if name == "YVNDLAGHEK" or index % 3 == c else 2
                for k in range(isotopes):
                    intensity = scale * pattern[k]
                    if intensity >= 1.0:
                        peaks.append((channel_mz + k * SPACING, round(intensity)))
        peaks.sort()
        spectra.append(spectrum_xml(index, rt, peaks))
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


def write_idxml(path):
    ids = [("YVNDLAGHEK", PEPTIDES["YVNDLAGHEK"], 129.0, 30, "PH_1", "1.0e-03"),
           ("LVNELTEFAK", PEPTIDES["LVNELTEFAK"], 130.0, 31, "PH_0", "1.0e-03"),
           ("LVNELTEFAK(Label:2H(4))", PEPTIDES["LVNELTEFAK"] + LYS4 / 2, 132.0, 33, "PH_0", "2.0e-03")]
    entries = "".join(f"""		<PeptideIdentification score_type="Posterior Error Probability" higher_score_better="false" significance_threshold="0.0" MZ="{mz:.15f}" RT="{rt}" spectrum_reference="scan={scan}" >
			<PeptideHit score="{score}" sequence="{seq}" charge="2" protein_refs="{prot}" >
				<UserParam type="string" name="target_decoy" value="target"/>
			</PeptideHit>
		</PeptideIdentification>
""" for seq, mz, rt, scan, prot, score in ids)
    with open(path, "w") as out:
        out.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<IdXML version="1.5" xsi:noNamespaceSchemaLocation="https://www.openms.de/xml-schema/IdXML_1_5.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
	<SearchParameters id="SP_0" db="MS1LabeledWorkflow_8.fasta" db_version="" taxonomy="" mass_type="monoisotopic" charges="" enzyme="unknown_enzyme" missed_cleavages="0" precursor_peak_tolerance="0" precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0" peak_mass_tolerance_ppm="false" >
		<VariableModification name="Label:2H(4) (K)" />
		<VariableModification name="Label:13C(6) (R)" />
		<VariableModification name="Label:13C(6)15N(2) (K)" />
		<VariableModification name="Label:13C(6)15N(4) (R)" />
	</SearchParameters>
	<IdentificationRun date="0000-00-00T00:00:00" search_engine="synthetic" search_engine_version="" search_parameters_ref="SP_0" >
		<ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0.0" >
			<ProteinHit id="PH_0" accession="PROT_A" score="0.0" sequence="" >
			</ProteinHit>
			<ProteinHit id="PH_1" accession="PROT_B" score="0.0" sequence="" >
			</ProteinHit>
			<UserParam type="stringList" name="spectra_data" value="[MS1LabeledWorkflow_24_input.mzML]"/>
		</ProteinIdentification>
{entries}	</IdentificationRun>
</IdXML>
""")


if __name__ == "__main__":
    write_mzml(os.path.join(HERE, "MS1LabeledWorkflow_24_input.mzML"))
    write_idxml(os.path.join(HERE, "MS1LabeledWorkflow_24.idXML"))
