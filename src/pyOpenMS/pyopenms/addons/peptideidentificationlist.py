"""Addon methods for PeptideIdentificationList class."""

from __future__ import annotations
import numpy as np
import warnings
from . import addon


@addon("PeptideIdentificationList")
def to_df(self, decode_ontology=True, default_missing_values=None, export_unidentified=True, columns=None):
    """Converts peptide identifications to a pandas DataFrame.

    :param decode_ontology: Decode meta value names using the PSI-MS ontology.
    :param default_missing_values: Default values for missing data by type.
    :param export_unidentified: Export PeptideIdentifications without PeptideHit.
    :param columns: List of column names to include. If None, includes all.
    :return: DataFrame with peptide identification information.
    """
    import pandas as pd
    import pyopenms

    if default_missing_values is None:
        default_missing_values = {bool: False, int: -9999, float: np.nan, str: ''}

    switchDict = {bool: '?', int: 'i', float: 'f', str: 'U100'}

    count = len(self)
    if not export_unidentified:
        count = sum(len(pep.getHits()) > 0 for pep in self)

    # get all possible metavalues
    metavals = []
    types = []
    mainscorename = "score"
    for pep in self:
        hits = pep.getHits()
        if len(hits) != 0:
            mvs = []
            hits[0].getKeys(mvs)
            metavals += mvs
            mainscorename = pep.getScoreType()

    metavals = list(set(metavals))

    # get type of all metavalues
    for k in metavals:
        k_str = k.decode('utf-8') if isinstance(k, bytes) else k
        if k_str == "target_decoy":
            types.append('?')
        else:
            found = False
            for p in self:
                hits = p.getHits()
                if len(hits) != 0:
                    if hits[0].metaValueExists(k_str):
                        mv = hits[0].getMetaValue(k_str)
                        types.append(switchDict[type(mv)])
                        found = True
                        break
            if not found:
                types.append('U100')

    # get default value for each type
    def get_key(val):
        for key, value in switchDict.items():
            if val == value:
                return key
        return str
    dmv = [default_missing_values[get_key(t)] for t in types]

    decodedMVs = [(m.decode("utf-8") if isinstance(m, bytes) else m) for m in metavals]
    if decode_ontology:
        try:
            cv = pyopenms.ControlledVocabulary()
            cv.loadFromOBO("psims", pyopenms.File.getOpenMSDataPath() + "/CV/psi-ms.obo")
            clearMVs = [cv.getTerm(m).name if m.startswith("MS:") else m for m in decodedMVs]
        except Exception:
            clearMVs = decodedMVs
    else:
        clearMVs = decodedMVs

    clearcols = ["id", "rt", "mz", mainscorename, "charge", "protein_accession", "start", "end", "P_ID", "PSM_ID"] + clearMVs
    coltypes = ['U100', 'f', 'f', 'f', 'i', 'U1000', 'U1000', 'U1000', 'i', 'i'] + types
    dt = list(zip(clearcols, coltypes))

    def extract(pep, pep_idx):
        hits = pep.getHits()
        if not hits:
            return (pep.getIdentifier(), pep.getRT(), pep.getMZ(), default_missing_values[float], default_missing_values[int],
                    default_missing_values[str], default_missing_values[str], default_missing_values[str], pep_idx, default_missing_values[int], *dmv)

        besthit = hits[0]
        ret = [pep.getIdentifier(), pep.getRT(), pep.getMZ(), besthit.getScore(), besthit.getCharge()]
        evs = besthit.getPeptideEvidences()
        ret += [','.join(v) if v else default_missing_values[str] for v in ([e.getProteinAccession() for e in evs],
                                                                            [str(e.getStart()) for e in evs],
                                                                            [str(e.getEnd()) for e in evs])]
        ret += [pep_idx, 0]

        for idx, k in enumerate(metavals):
            k_str = k.decode('utf-8') if isinstance(k, bytes) else k
            if besthit.metaValueExists(k_str):
                val = besthit.getMetaValue(k_str)
                if k_str == "target_decoy":
                    if isinstance(val, str) and val and val[0] == 't':
                        ret.append(True)
                    else:
                        ret.append(False)
                else:
                    ret.append(val)
            else:
                ret.append(dmv[idx])
        return tuple(ret)

    if export_unidentified:
        rows = (extract(pep, pep_idx) for pep_idx, pep in enumerate(self))
    else:
        rows = (extract(pep, pep_idx) for pep_idx, pep in enumerate(self) if pep.getHits())

    df = pd.DataFrame(np.fromiter(rows, dtype=dt, count=count))

    # Filter columns if specified
    if columns is not None:
        available = set(df.columns)
        unknown = [c for c in columns if c not in available]
        if unknown:
            warnings.warn(f"Unknown column(s) ignored in to_df(): {unknown}. "
                          f"Use df_columns() to see available columns.")
        requested = [c for c in columns if c in available]
        df = df[requested]

    return df


@addon("PeptideIdentificationList")
def df_columns(self, decode_ontology=True):
    """Returns a list of column names that to_df() would produce.

    :param decode_ontology: Decode meta value names using the PSI-MS ontology.
    :return: List of column name strings.
    """
    import pyopenms

    # Get all possible metavalues and main score name
    metavals = []
    mainscorename = "score"
    for pep in self:
        hits = pep.getHits()
        if hits:
            mvs = []
            hits[0].getKeys(mvs)
            metavals += mvs
            mainscorename = pep.getScoreType()

    metavals = list(set(metavals))

    # Decode metavalue names if requested
    decodedMVs = [(m.decode("utf-8") if isinstance(m, bytes) else m) for m in metavals]
    if decode_ontology:
        try:
            cv = pyopenms.ControlledVocabulary()
            cv.loadFromOBO("psims", pyopenms.File.getOpenMSDataPath() + "/CV/psi-ms.obo")
            clearMVs = [cv.getTerm(m).name if m.startswith("MS:") else m for m in decodedMVs]
        except Exception:
            clearMVs = decodedMVs
    else:
        clearMVs = decodedMVs

    return ["id", "rt", "mz", mainscorename, "charge", "protein_accession", "start", "end", "P_ID", "PSM_ID"] + clearMVs


@addon("PeptideIdentificationList")
def to_arrow(self, decode_ontology=True, default_missing_values=None, export_unidentified=True, columns=None):
    """Returns an Apache Arrow Table with peptide identification information."""
    import pyarrow as pa
    df = self.to_df(decode_ontology=decode_ontology,
                    default_missing_values=default_missing_values,
                    export_unidentified=export_unidentified,
                    columns=columns)

    # For empty tables, define explicit schema to preserve types
    if len(df) == 0:
        # Core column types - must match df_columns() order
        core_schema = [
            pa.field("id", pa.utf8()),
            pa.field("rt", pa.float32()),
            pa.field("mz", pa.float32()),
            pa.field("score", pa.float32()),
            pa.field("charge", pa.int32()),
            pa.field("protein_accession", pa.utf8()),
            pa.field("start", pa.utf8()),
            pa.field("end", pa.utf8()),
            pa.field("P_ID", pa.int32()),
            pa.field("PSM_ID", pa.int32()),
        ]
        # Filter to requested columns if specified
        if columns is not None:
            core_schema = [f for f in core_schema if f.name in columns]
        else:
            # Only include columns that are in the DataFrame
            core_schema = [f for f in core_schema if f.name in df.columns]
        schema = pa.schema(core_schema)
        return pa.Table.from_pandas(df, schema=schema)

    return pa.Table.from_pandas(df)


@addon("PeptideIdentificationList")
def to_psm_df(self, export_all_hits=True, include_modifications=True, include_peak_annotations=False,
              reference_file_name="", columns=None, additional_score_names=None, scan_format="scan"):
    """
    **EXPERIMENTAL**: Export all PSMs (Peptide-Spectrum Matches) as a DataFrame.

    Unlike to_df() which exports only the top hit per spectrum, this method
    exports all hits with proper ranking information, modifications, and
    additional scores.

    :param export_all_hits: If True, export all hits per identification.
    :param include_modifications: Include detailed modification list column.
    :param include_peak_annotations: Include fragment ion annotations.
    :param reference_file_name: Source file name to include in each row.
    :param columns: List of column names to include. If None, includes all.
    :param additional_score_names: List of additional metavalue names to treat as scores.
    :param scan_format: Controls scan column: "scan" or "nativeId".
    :return: DataFrame with PSM data.
    """
    table = self.to_psm_arrow(
        export_all_hits=export_all_hits,
        include_modifications=include_modifications,
        include_peak_annotations=include_peak_annotations,
        reference_file_name=reference_file_name,
        columns=columns,
        additional_score_names=additional_score_names,
        scan_format=scan_format
    )
    return table.to_pandas()


@addon("PeptideIdentificationList")
def to_qpx(self, qpx_version="1.0", creator="pyOpenMS", software_provider="OpenMS",
           scan_format="scan", **kwargs):
    """
    **EXPERIMENTAL**: Export PSMs as QPX format structure with file metadata and PSMs array.

    :param qpx_version: Version of the QPX format. Default "1.0".
    :param creator: Name of the tool or person who created the file.
    :param software_provider: Name of the software provider.
    :param scan_format: Controls the scan column format.
    :return: Dict with 'file_metadata' and 'psms' keys following QPX schema.
    """
    import uuid as _uuid
    from datetime import datetime as _datetime

    df = self.to_psm_df(scan_format=scan_format, **kwargs)

    # Build file_metadata
    file_metadata = {
        "qpx_version": qpx_version,
        "creator": creator,
        "file_type": "psm",
        "creation_date": _datetime.now().isoformat(),
        "uuid": str(_uuid.uuid4()),
        "scan_format": scan_format,
        "software_provider": software_provider,
    }

    # Convert DataFrame to list of dicts for psms array
    psms = df.to_dict(orient="records")

    return {
        "file_metadata": file_metadata,
        "psms": psms
    }


@addon("PeptideIdentificationList")
def to_psm_arrow(self, export_all_hits=True, include_modifications=True,
                 include_peak_annotations=False, reference_file_name="",
                 columns=None, additional_score_names=None, scan_format="scan"):
    """
    **EXPERIMENTAL**: Export all PSMs as Apache Arrow Table.

    This method builds the Arrow table directly without pandas intermediate,
    making it more memory-efficient for large datasets.

    :param export_all_hits: If True, export all hits per identification.
    :param include_modifications: Include detailed modification list column.
    :param include_peak_annotations: Include fragment ion annotations.
    :param reference_file_name: Source file name to include in each row.
    :param columns: List of column names to include. If None, includes all.
    :param additional_score_names: List of additional metavalue names to treat as scores.
    :param scan_format: Controls scan column: "scan" or "nativeId".
    :return: Arrow Table with PSM data.
    """
    import pyarrow as pa
    import pyopenms

    if scan_format not in ("scan", "nativeId"):
        raise ValueError(f"scan_format must be 'scan' or 'nativeId', got '{scan_format}'")

    # Native ID type accessions for scan number extraction
    _native_id_accessions = [
        "MS:1000768", "MS:1000769", "MS:1000771", "MS:1001508",
        "MS:1000774", "MS:1000777", "MS:1001530",
    ]

    # Create SpectrumLookup once for all scan extractions
    _spectrum_lookup = pyopenms.SpectrumLookup()

    def _extract_scan_number(native_id):
        if not native_id:
            return None
        for accession in _native_id_accessions:
            try:
                scan_num = _spectrum_lookup.extractScanNumber(native_id, accession)
                if scan_num >= 0:
                    return scan_num
            except Exception:
                continue
        return None

    # User-specified additional score names (if any)
    _additional_score_names_set = set(additional_score_names) if additional_score_names else set()

    def _is_known_score(name):
        """Check if a name is a known score type."""
        if name in _additional_score_names_set:
            return True
        return pyopenms.Scores.isKnownScoreType(name)

    def _get_value_type(val):
        """Determine the type string for a metavalue."""
        if type(val).__name__ == 'bool':
            return "bool"
        elif isinstance(val, int):
            return "int"
        elif isinstance(val, float):
            return "float"
        else:
            return "str"

    # Excluded metavalues (have dedicated columns)
    _excluded_psm_metavalues = {"target_decoy", "predicted_RT", "predicted_rt"}
    _excluded_spectrum_metavalues = {"spectrum_reference", "ion_mobility", "IM"}

    # Initialize column lists
    all_sequence = []
    all_peptidoform = []
    all_modifications = []
    all_precursor_charge = []
    all_pep = []
    all_is_decoy = []
    all_calculated_mz = []
    all_observed_mz = []
    all_additional_scores = []
    all_protein_accessions = []
    all_predicted_rt = []
    all_reference_file = []
    all_cv_params = []
    all_scan = []
    all_rt = []
    all_ion_mobility = []
    all_spectrum_reference = []
    all_score = []
    all_score_type = []
    all_hit_index = []
    all_p_id = []

    # Peak annotation arrays
    all_number_peaks = []
    all_mz_array = []
    all_intensity_array = []
    all_charge_array = []
    all_ion_type_array = []
    all_ion_mobility_array = []

    # Metavalue arrays
    all_psm_metavalues = []
    all_spectrum_metavalues = []

    # Create IDScoreSwitcherAlgorithm instance
    idsa = pyopenms.IDScoreSwitcherAlgorithm()

    for pep_idx, pep_id in enumerate(self):
        rt = pep_id.getRT()
        observed_mz = pep_id.getMZ()

        # Spectrum reference
        spec_ref = ""
        if pep_id.metaValueExists("spectrum_reference"):
            spec_ref = pep_id.getMetaValue("spectrum_reference")
            if isinstance(spec_ref, bytes):
                spec_ref = spec_ref.decode("utf-8")

        score_type = pep_id.getScoreType()

        # Ion mobility
        ion_mobility = None
        if pep_id.metaValueExists("ion_mobility"):
            ion_mobility = pep_id.getMetaValue("ion_mobility")
        elif pep_id.metaValueExists("IM"):
            ion_mobility = pep_id.getMetaValue("IM")

        # Extract scan number
        scan = _extract_scan_number(spec_ref)

        # Collect spectrum-level metavalues (PeptideIdentification level)
        spectrum_metavalues = []
        pep_id_keys = []
        pep_id.getKeys(pep_id_keys)
        for key in pep_id_keys:
            key_str = key.decode("utf-8") if isinstance(key, bytes) else str(key)
            if key_str not in _excluded_spectrum_metavalues:
                val = pep_id.getMetaValue(key_str)
                val_type = _get_value_type(val)
                if isinstance(val, bytes):
                    val = val.decode("utf-8")
                spectrum_metavalues.append({
                    "name": key_str,
                    "value": val,
                    "value_type": val_type
                })

        hits = pep_id.getHits()
        if not hits:
            continue

        num_hits = len(hits) if export_all_hits else min(1, len(hits))

        # Detect PEP score type using findScoreType
        pep_search_result = idsa.findScoreType(pep_id, pyopenms.IDType.PEP)
        pep_score_name = pep_search_result.score_name

        for hit_idx in range(num_hits):
            hit = hits[hit_idx]
            seq = hit.getSequence()
            charge = hit.getCharge()

            # Extract modifications
            modifications = []
            if include_modifications and seq.isModified():
                mod_dict = {}

                # N-terminal modification
                if seq.hasNTerminalModification():
                    mod = seq.getNTerminalModification()
                    mod_name = seq.getNTerminalModificationName()
                    accession = mod.getUniModRecordId() if mod else None
                    accession_str = f"UNIMOD:{accession}" if accession and accession > 0 else None
                    position_str = "N-term.0"
                    if mod_name not in mod_dict:
                        mod_dict[mod_name] = {"name": mod_name, "accession": accession_str, "positions": []}
                    mod_dict[mod_name]["positions"].append({"position": position_str, "scores": None})

                # Residue modifications
                for pos, residue in enumerate(seq):
                    if residue.isModified():
                        res_mod = residue.getModification()
                        mod_name = residue.getModificationName()
                        accession = res_mod.getUniModRecordId() if res_mod else None
                        accession_str = f"UNIMOD:{accession}" if accession and accession > 0 else None
                        aa_code = residue.getOneLetterCode()
                        position_str = f"{aa_code}.{pos + 1}"
                        if mod_name not in mod_dict:
                            mod_dict[mod_name] = {"name": mod_name, "accession": accession_str, "positions": []}
                        mod_dict[mod_name]["positions"].append({"position": position_str, "scores": None})

                # C-terminal modification
                if seq.hasCTerminalModification():
                    mod = seq.getCTerminalModification()
                    mod_name = seq.getCTerminalModificationName()
                    accession = mod.getUniModRecordId() if mod else None
                    accession_str = f"UNIMOD:{accession}" if accession and accession > 0 else None
                    position_str = f"C-term.{len(seq) + 1}"
                    if mod_name not in mod_dict:
                        mod_dict[mod_name] = {"name": mod_name, "accession": accession_str, "positions": []}
                    mod_dict[mod_name]["positions"].append({"position": position_str, "scores": None})

                modifications = list(mod_dict.values())

            # Determine is_decoy
            is_decoy = None
            if hit.metaValueExists("target_decoy"):
                td = hit.getMetaValue("target_decoy")
                if isinstance(td, bytes):
                    is_decoy = 1 if td.startswith(b"decoy") else 0
                else:
                    is_decoy = 1 if str(td).startswith("decoy") else 0

            # Protein accessions
            evidences = hit.getPeptideEvidences()
            protein_accessions = [ev.getProteinAccession() for ev in evidences]

            # Build additional_scores and psm_metavalues
            additional_scores = []
            psm_metavalues = []
            keys = []
            hit.getKeys(keys)
            for key in keys:
                key_str = key.decode("utf-8") if isinstance(key, bytes) else str(key)
                if key_str not in _excluded_psm_metavalues:
                    val = hit.getMetaValue(key_str)
                    is_known_score = _is_known_score(key_str)

                    if isinstance(val, (int, float)) and is_known_score:
                        higher_better = None
                        try:
                            score_type_enum = pyopenms.IDScoreSwitcherAlgorithm.toScoreTypeEnum(key_str)
                            higher_better = idsa.isScoreTypeHigherBetter(score_type_enum)
                        except Exception:
                            pass
                        additional_scores.append({
                            "score_name": key_str,
                            "score_value": float(val),
                            "higher_better": higher_better
                        })
                    else:
                        val_type = _get_value_type(val)
                        if isinstance(val, bytes):
                            val = val.decode("utf-8")
                        psm_metavalues.append({
                            "name": key_str,
                            "value": val,
                            "value_type": val_type
                        })

            # Calculate m/z
            calculated_mz = None
            if charge > 0:
                try:
                    calculated_mz = seq.getMZ(charge)
                except Exception:
                    pass

            # Predicted RT
            predicted_rt = None
            if hit.metaValueExists("predicted_RT"):
                predicted_rt = hit.getMetaValue("predicted_RT")
            elif hit.metaValueExists("predicted_rt"):
                predicted_rt = hit.getMetaValue("predicted_rt")

            # PEP value
            pep_value = None
            if pep_score_name:
                if pep_search_result.is_main_score_type:
                    pep_value = hit.getScore()
                else:
                    if hit.metaValueExists(pep_score_name):
                        pep_value = hit.getMetaValue(pep_score_name)

            # Append to column lists
            all_sequence.append(seq.toUnmodifiedString())
            all_peptidoform.append(seq.toString())
            all_modifications.append(modifications if include_modifications else None)
            all_precursor_charge.append(charge)
            all_pep.append(pep_value)
            all_is_decoy.append(is_decoy)
            all_calculated_mz.append(calculated_mz)
            all_observed_mz.append(observed_mz)
            all_additional_scores.append(additional_scores)
            all_protein_accessions.append(protein_accessions)
            all_predicted_rt.append(predicted_rt)
            all_reference_file.append(reference_file_name)
            all_cv_params.append(None)
            if scan_format == "nativeId":
                all_scan.append(spec_ref if spec_ref else None)
            else:
                all_scan.append(str(scan) if scan is not None else None)
            all_rt.append(rt)
            all_ion_mobility.append(ion_mobility)
            all_spectrum_reference.append(spec_ref)
            all_score.append(hit.getScore())
            all_score_type.append(score_type)
            all_hit_index.append(hit_idx)
            all_p_id.append(pep_idx)

            # Metavalues
            all_psm_metavalues.append(psm_metavalues)
            all_spectrum_metavalues.append(spectrum_metavalues)

            # Peak annotations
            if include_peak_annotations:
                peak_annotations = hit.getPeakAnnotations()
                if peak_annotations:
                    all_number_peaks.append(len(peak_annotations))
                    all_mz_array.append([ann.mz for ann in peak_annotations])
                    all_intensity_array.append([ann.intensity for ann in peak_annotations])
                    all_charge_array.append([ann.charge for ann in peak_annotations])
                    ion_types = []
                    for ann in peak_annotations:
                        if isinstance(ann.annotation, bytes):
                            ion_types.append(ann.annotation.decode("utf-8"))
                        else:
                            ion_types.append(ann.annotation)
                    all_ion_type_array.append(ion_types)
                    all_ion_mobility_array.append(None)
                else:
                    all_number_peaks.append(0)
                    all_mz_array.append([])
                    all_intensity_array.append([])
                    all_charge_array.append([])
                    all_ion_type_array.append([])
                    all_ion_mobility_array.append(None)
            else:
                all_number_peaks.append(None)
                all_mz_array.append(None)
                all_intensity_array.append(None)
                all_charge_array.append(None)
                all_ion_type_array.append(None)
                all_ion_mobility_array.append(None)

    # Build Arrow table with proper nested types
    columns_set = set(columns) if columns is not None else None

    def should_include(col_name):
        return columns_set is None or col_name in columns_set

    # Define nested types for complex columns
    position_type = pa.struct([
        ("position", pa.utf8()),
        ("scores", pa.float64())
    ])
    modification_type = pa.struct([
        ("name", pa.utf8()),
        ("accession", pa.utf8()),
        ("positions", pa.list_(position_type))
    ])
    score_type_struct = pa.struct([
        ("score_name", pa.utf8()),
        ("score_value", pa.float64()),
        ("higher_better", pa.bool_())
    ])
    metavalue_type = pa.struct([
        ("name", pa.utf8()),
        ("value", pa.utf8()),
        ("value_type", pa.utf8())
    ])

    data_dict = {}

    # Core columns
    if should_include("sequence"):
        data_dict["sequence"] = pa.array(all_sequence, type=pa.utf8())
    if should_include("peptidoform"):
        data_dict["peptidoform"] = pa.array(all_peptidoform, type=pa.utf8())
    if should_include("modifications"):
        data_dict["modifications"] = pa.array(all_modifications, type=pa.list_(modification_type))
    if should_include("precursor_charge"):
        data_dict["precursor_charge"] = pa.array(all_precursor_charge, type=pa.int32())
    if should_include("posterior_error_probability"):
        data_dict["posterior_error_probability"] = pa.array(all_pep, type=pa.float64())
    if should_include("is_decoy"):
        data_dict["is_decoy"] = pa.array(all_is_decoy, type=pa.int32())
    if should_include("calculated_mz"):
        data_dict["calculated_mz"] = pa.array(all_calculated_mz, type=pa.float64())
    if should_include("observed_mz"):
        data_dict["observed_mz"] = pa.array(all_observed_mz, type=pa.float64())
    if should_include("additional_scores"):
        data_dict["additional_scores"] = pa.array(all_additional_scores, type=pa.list_(score_type_struct))
    if should_include("protein_accessions"):
        data_dict["protein_accessions"] = pa.array(all_protein_accessions, type=pa.list_(pa.utf8()))
    if should_include("predicted_rt"):
        data_dict["predicted_rt"] = pa.array(all_predicted_rt, type=pa.float64())
    if should_include("reference_file_name"):
        data_dict["reference_file_name"] = pa.array(all_reference_file, type=pa.utf8())
    if should_include("cv_params"):
        data_dict["cv_params"] = pa.array(all_cv_params, type=pa.utf8())
    if should_include("scan"):
        data_dict["scan"] = pa.array(all_scan, type=pa.utf8())
    if should_include("rt"):
        data_dict["rt"] = pa.array(all_rt, type=pa.float64())
    if should_include("ion_mobility"):
        data_dict["ion_mobility"] = pa.array(all_ion_mobility, type=pa.float64())

    # Peak annotation columns
    if should_include("number_peaks"):
        data_dict["number_peaks"] = pa.array(all_number_peaks, type=pa.int32())
    if should_include("mz_array"):
        data_dict["mz_array"] = pa.array(all_mz_array, type=pa.list_(pa.float64()))
    if should_include("intensity_array"):
        data_dict["intensity_array"] = pa.array(all_intensity_array, type=pa.list_(pa.float32()))
    if should_include("charge_array"):
        data_dict["charge_array"] = pa.array(all_charge_array, type=pa.list_(pa.int32()))
    if should_include("ion_type_array"):
        data_dict["ion_type_array"] = pa.array(all_ion_type_array, type=pa.list_(pa.utf8()))
    if should_include("ion_mobility_array"):
        data_dict["ion_mobility_array"] = pa.array(all_ion_mobility_array, type=pa.list_(pa.float64()))

    # OpenMS-specific columns
    if should_include("spectrum_reference"):
        data_dict["spectrum_reference"] = pa.array(all_spectrum_reference, type=pa.utf8())
    if should_include("score"):
        data_dict["score"] = pa.array(all_score, type=pa.float64())
    if should_include("score_type"):
        data_dict["score_type"] = pa.array(all_score_type, type=pa.utf8())
    if should_include("hit_index"):
        data_dict["hit_index"] = pa.array(all_hit_index, type=pa.int32())
    if should_include("P_ID"):
        data_dict["P_ID"] = pa.array(all_p_id, type=pa.int32())

    # Metavalue columns - stringify values for Arrow compatibility
    if should_include("psm_metavalues"):
        all_psm_metavalues_str = [
            [{"name": mv["name"],
              "value": str(mv["value"]) if mv["value"] is not None else None,
              "value_type": mv["value_type"]}
             for mv in mvs] for mvs in all_psm_metavalues
        ]
        data_dict["psm_metavalues"] = pa.array(all_psm_metavalues_str, type=pa.list_(metavalue_type))
    if should_include("spectrum_metavalues"):
        all_spectrum_metavalues_str = [
            [{"name": mv["name"],
              "value": str(mv["value"]) if mv["value"] is not None else None,
              "value_type": mv["value_type"]}
             for mv in mvs] for mvs in all_spectrum_metavalues
        ]
        data_dict["spectrum_metavalues"] = pa.array(all_spectrum_metavalues_str, type=pa.list_(metavalue_type))

    table = pa.Table.from_pydict(data_dict)

    if columns_set is not None:
        unknown = [c for c in columns if c not in data_dict]
        if unknown:
            warnings.warn(f"Unknown column(s) ignored in to_psm_arrow(): {unknown}. "
                          f"Use psm_columns() to see available columns.")

    # Sort by rt, observed_mz, precursor_charge, hit_index for consistent ordering
    sort_cols = ["rt", "observed_mz", "precursor_charge", "hit_index"]
    available_sort_cols = [(c, "ascending") for c in sort_cols if c in table.column_names]
    if available_sort_cols and table.num_rows > 0:
        import pyarrow.compute as pc
        indices = pc.sort_indices(table, sort_keys=available_sort_cols, null_placement="at_end")
        table = table.take(indices)

    return table


@addon("PeptideIdentificationList")
def to_parquet(self, path, compression='zstd', compression_level=None,
               row_group_size=None, write_statistics=True,
               export_all_hits=True, include_modifications=True,
               include_peak_annotations=False, reference_file_name="",
               columns=None, additional_score_names=None, scan_format="scan"):
    """
    **EXPERIMENTAL**: Export all PSMs to a Parquet file using native Arrow construction.

    :param path: Path to the output Parquet file.
    :param compression: Compression algorithm ('none', 'snappy', 'gzip', 'lz4', 'zstd').
    :param compression_level: Compression level (ZSTD: 1-22, GZIP: 1-9).
    :param row_group_size: Number of rows per row group.
    :param write_statistics: Write column statistics.
    :param export_all_hits: If True, export all hits per identification.
    :param include_modifications: Include modification details.
    :param include_peak_annotations: Include fragment ion annotations.
    :param reference_file_name: Source file name to include.
    :param columns: Column names to include. None includes all.
    :param additional_score_names: Additional metavalue names to treat as scores.
    :param scan_format: Controls scan column format.
    """
    import pyarrow.parquet as pq

    # Build the Arrow table
    table = self.to_psm_arrow(
        export_all_hits=export_all_hits,
        include_modifications=include_modifications,
        include_peak_annotations=include_peak_annotations,
        reference_file_name=reference_file_name,
        columns=columns,
        additional_score_names=additional_score_names,
        scan_format=scan_format
    )

    # Map compression string to pyarrow compression
    compression_map = {
        'none': None,
        'snappy': 'snappy',
        'gzip': 'gzip',
        'lz4': 'lz4',
        'zstd': 'zstd',
    }
    compression_lower = compression.lower() if compression else 'none'
    if compression_lower not in compression_map:
        raise ValueError(f"Unsupported compression: {compression}. "
                         f"Use one of: {list(compression_map.keys())}")

    pq_compression = compression_map[compression_lower]

    # Build write options dict
    write_kwargs = {
        'compression': pq_compression,
        'write_statistics': write_statistics,
    }

    if compression_level is not None and pq_compression in ('zstd', 'gzip'):
        write_kwargs['compression_level'] = compression_level

    if row_group_size is not None:
        write_kwargs['row_group_size'] = row_group_size

    with open(path, 'wb') as f:
        pq.write_table(table, f, **write_kwargs)


@addon("PeptideIdentificationList")
def psm_columns(self):
    """
    **EXPERIMENTAL**: Return list of column names that to_psm_df() and to_qpx() produce.

    :return: List of column name strings.
    """
    return [
        # QPX schema fields (in schema order)
        "sequence",
        "peptidoform",
        "modifications",
        "precursor_charge",
        "posterior_error_probability",
        "is_decoy",
        "calculated_mz",
        "observed_mz",
        "additional_scores",
        "protein_accessions",
        "predicted_rt",
        "reference_file_name",
        "cv_params",
        "scan",
        "rt",
        "ion_mobility",
        "number_peaks",
        "mz_array",
        "intensity_array",
        "charge_array",
        "ion_type_array",
        "ion_mobility_array",
        # OpenMS-specific fields (not in QPX)
        "spectrum_reference",
        "score",
        "score_type",
        "hit_index",
        "P_ID",
        # Metavalues
        "psm_metavalues",
        "spectrum_metavalues",
    ]


@addon("PeptideIdentificationList")
def update_scores_from_df(self, df, main_score_name):
    """Updates scores in PeptideIdentification objects from a DataFrame."""
    for index, row in df.iterrows():
        pid_index = int(row["P_ID"])
        pi = self[pid_index]
        pi.setScoreType(main_score_name)
        hits = pi.getHits()
        if len(hits) > 0:
            best_hit = hits[0]
            best_hit.setScore(float(row[main_score_name]))
            pi.setHits([best_hit] + list(hits[1:]))

    return self


@addon("PeptideIdentificationList")
def get_df(self, *args, **kwargs):
    """Deprecated: Use to_df() instead."""
    warnings.warn(
        "get_df() is deprecated. Use to_df() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.to_df(*args, **kwargs)


@addon("PeptideIdentificationList")
def to_psm_qpx(self, *args, **kwargs):
    """
    **EXPERIMENTAL**: Export PSMs as QPX format structure with file metadata and PSMs array.

    This is an alias for to_qpx().
    """
    return self.to_qpx(*args, **kwargs)


@addon("PeptideIdentificationList")
def get_psm_columns(self, *args, **kwargs):
    """Deprecated: Use psm_columns() instead."""
    warnings.warn(
        "get_psm_columns() is deprecated. Use psm_columns() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.psm_columns(*args, **kwargs)


@addon("PeptideIdentificationList")
def get_psm_df(self, *args, **kwargs):
    """Deprecated: Use to_psm_df() instead."""
    warnings.warn(
        "get_psm_df() is deprecated. Use to_psm_df() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.to_psm_df(*args, **kwargs)
