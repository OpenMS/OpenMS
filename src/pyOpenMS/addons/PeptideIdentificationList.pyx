cimport numpy as np
import numpy as np




    def __len__(self):
        """
        __len__(self: PeptideIdentificationList) -> int

        Return the number of peptide identifications in the list.

        Enables use of Python's built-in len() function.

        :return: The number of peptide identifications in this list.
        :rtype: int
        """
        return self.size()

    def append(self, PeptideIdentification item):
        """
        append(self: PeptideIdentificationList, item: PeptideIdentification) -> None

        Add a single peptide identification to the end of the list.

        This method provides a Pythonic interface equivalent to push_back().

        :param item: A single PeptideIdentification object to append.
        :type item: PeptideIdentification
        """
        self.push_back(item)

    def extend(self, items):
        """
        extend(self: PeptideIdentificationList, items: Iterable[PeptideIdentification]) -> None

        Add multiple peptide identifications to the end of the list.

        :param items: Can be:
                      - A list/iterable of PeptideIdentification objects
                      - Another PeptideIdentificationList object
        :type items: Iterable[PeptideIdentification]
        :raises TypeError: If items is not iterable or another PeptideIdentificationList.
        """
        if hasattr(items, '__iter__') and not hasattr(items, 'inst'):
            # Handle regular iterables (list, tuple, etc.)
            for peptide_identification in items:
                self.push_back(peptide_identification)
        elif hasattr(items, 'inst') and hasattr(items, '__len__'):
            # Handle another PeptideIdentificationList
            for peptide_identification in items:
                self.push_back(peptide_identification)
        else:
            raise TypeError("extend() argument must be iterable or another PeptideIdentificationList")

    def __str__(self):
        """
        __str__(self: PeptideIdentificationList) -> str
        
        Return a string representation of the PeptideIdentificationList object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: PeptideIdentificationList) -> str
        
        Return a string representation of the PeptideIdentificationList object.

        Returns key properties in a readable format:
        PeptideIdentificationList(size=5)
        """
        return f"PeptideIdentificationList(size={self.size()})"

    def to_df(self, decode_ontology=True, default_missing_values=None, export_unidentified=True):
        """
        to_df(self: PeptideIdentificationList, decode_ontology: bool = True, default_missing_values: dict = None, export_unidentified: bool = True) -> pd.DataFrame

        Converts the peptide identifications to a pandas DataFrame.

        :param decode_ontology: Decode meta value names using the PSI-MS ontology.
                                Default True.
        :type decode_ontology: bool

        :param default_missing_values: Default values for missing data by type.
                                       Default: {<bool>: False, <int>: -9999, <float>: nan, <str>: ''}
        :type default_missing_values: dict

        :param export_unidentified: Export PeptideIdentifications without PeptideHit.
                                    Default True.
        :type export_unidentified: bool

        :return: DataFrame with peptide identification information including:
                 id, rt, mz, score (column name from score type), charge,
                 protein_accession (comma-separated), start (comma-separated),
                 end (comma-separated), P_ID (peptide identification index),
                 PSM_ID (PSM index, currently always 0, only first hit exported),
                 and additional meta value columns from PeptideHit.
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed

        Example::

            peps = feature_map.get_assigned_peptide_identifications()
            df = peps.to_df()

            # Skip unidentified entries
            df = peps.to_df(export_unidentified=False)

            # Custom missing values (use Python type objects)
            df = peps.to_df(default_missing_values={bool: False, int: 0, float: 0.0, str: 'NA'})
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for to_df(). "
                "Please install it with: pip install pandas"
            )
        from . import ControlledVocabulary as _ControlledVocabulary
        from . import File as _File

        if default_missing_values is None:
            # Use type() to get Python types since 'bool' conflicts with C bool in Cython
            default_missing_values = {type(True): False, type(1): -9999, type(1.0): np.nan, type(''): ''}

        switchDict = {type(True): '?', type(1): 'i', type(1.0): 'f', type(''): 'U100'}

        # filter out PeptideIdentifications without PeptideHits if export_unidentified == False
        count = self.size()
        if not export_unidentified:
            count = sum(len(pep.getHits()) > 0 for pep in self)

        # get all possible metavalues
        metavals = []
        types = []
        mainscorename = "score"
        for pep in self:
            hits = pep.getHits()
            if not len(hits) == 0:
                mvs = []
                hits[0].getKeys(mvs)
                metavals += mvs
                mainscorename = pep.getScoreType()

        metavals = list(set(metavals))

        # get type of all metavalues
        for k in metavals:
            if k == b"target_decoy":
                types.append('?')
            else:
                for p in self:
                    hits = p.getHits()
                    if not len(hits) == 0:
                        if hits[0].metaValueExists(k):
                            mv = hits[0].getMetaValue(k)
                            types.append(switchDict[type(mv)])
                            break

        # get default value for each type in types to append if there are no hits in a PeptideIdentification
        def get_key(val):
            for key, value in switchDict.items():
                if val == value:
                    return key
        dmv = [default_missing_values[get_key(t)] for t in types]

        decodedMVs = [m.decode("utf-8") for m in metavals]
        if decode_ontology:
            cv = _ControlledVocabulary()
            cv.loadFromOBO("psims", _File.getOpenMSDataPath() + "/CV/psi-ms.obo")
            clearMVs = [cv.getTerm(m).name if m.startswith("MS:") else m for m in decodedMVs]
        else:
            clearMVs = decodedMVs

        clearcols = ["id", "rt", "mz", mainscorename, "charge", "protein_accession", "start", "end", "P_ID", "PSM_ID"] + clearMVs
        coltypes = ['U100', 'f', 'f', 'f', 'i','U1000', 'U1000', 'U1000', 'i', 'i'] + types
        dt = list(zip(clearcols, coltypes))

        def extract(pep, pep_idx):
            hits = pep.getHits()
            if not hits:
                if export_unidentified:
                    return (pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), default_missing_values[float], default_missing_values[int],
                            default_missing_values[str], default_missing_values[str], default_missing_values[str], pep_idx, default_missing_values[int], *dmv)
                else:
                    return

            besthit = hits[0]
            ret = [pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), besthit.getScore(), besthit.getCharge()]
            # add accession, start and end positions of peptide evidences as comma separated str (like in mzTab)
            evs = besthit.getPeptideEvidences()
            ret += [','.join(v) if v else default_missing_values[str] for v in ([e.getProteinAccession() for e in evs],
                                                                                [str(e.getStart()) for e in evs],
                                                                                [str(e.getEnd()) for e in evs])]

            ret += [str(pep_idx), 0]  # we currently only export the first hit

            for idx, k in enumerate(metavals):
                if besthit.metaValueExists(k):
                    val = besthit.getMetaValue(k)
                    if k == b"target_decoy":
                        if val[0] == 't':
                            ret.append(True)
                        else:
                            ret.append(False)
                    else:
                        ret.append(val)
                else:
                    ret.append(dmv[idx])  # Use precomputed default for this metavalue's type
            return tuple(ret)

        return pd.DataFrame(np.fromiter((extract(pep, pep_idx) for pep_idx, pep in enumerate(self)), dtype=dt, count=count))

    def to_arrow(self, decode_ontology=True, default_missing_values=None, export_unidentified=True):
        """
        to_arrow(self: PeptideIdentificationList, decode_ontology: bool = True, default_missing_values: dict = None, export_unidentified: bool = True) -> pa.Table

        Returns an Apache Arrow Table with peptide identification information.

        :param decode_ontology: Decode meta value names using the PSI-MS ontology.
                                Default True.
        :type decode_ontology: bool

        :param default_missing_values: Default values for missing data by type.
                                       Default: {<bool>: False, <int>: -9999, <float>: nan, <str>: ''}
        :type default_missing_values: dict

        :param export_unidentified: Export PeptideIdentifications without PeptideHit.
                                    Default True.
        :type export_unidentified: bool

        :return: Arrow Table with peptide identification data.
        :rtype: pyarrow.Table

        :raises ImportError: If pyarrow is not installed

        Example::

            table = peps.to_arrow()
            df = table.to_pandas()

            # Convert to polars
            import polars as pl
            df = pl.from_arrow(table)
        """
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_arrow(). "
                "Please install it with: pip install pyarrow"
            )
        df = self.to_df(decode_ontology=decode_ontology,
                        default_missing_values=default_missing_values,
                        export_unidentified=export_unidentified)
        return pa.Table.from_pandas(df)

    def update_scores_from_df(self, df, main_score_name):
        """
        update_scores_from_df(self: PeptideIdentificationList, df: pd.DataFrame, main_score_name: str) -> PeptideIdentificationList

        Updates the scores in PeptideIdentification objects using a pandas DataFrame.

        This method is useful when you've computed new scores (e.g., from rescoring)
        and want to update the original PeptideIdentifications.

        :param df: DataFrame obtained by converting this list to a DataFrame.
                   Minimum required columns: P_ID and a column with name
                   matching main_score_name.
        :type df: pd.DataFrame
        :param main_score_name: Name of the column containing the new scores.
                                This will also be set as the score type.
        :type main_score_name: str
        :return: The updated list of peptide identifications
                 (modifies in place and returns self).
        :rtype: PeptideIdentificationList

        Example::

            >>> # Get DataFrame, compute new scores, update
            >>> df = peps.to_df()
            >>> df['new_score'] = compute_new_scores(df)
            >>> peps.update_scores_from_df(df, 'new_score')
        """
        from . import PeptideIdentification as _PeptideIdentification

        for index, row in df.iterrows():
            pid_index = int(row["P_ID"])
            pi = _PeptideIdentification(self[pid_index])
            pi.setScoreType(main_score_name)
            hits = pi.getHits()
            if len(hits) > 0:
                best_hit = hits[0]
                best_hit.setScore(float(row[main_score_name]))
                hits[0] = best_hit
                pi.setHits(hits)

            self[pid_index] = pi

        return self

    # Deprecated alias for backward compatibility
    def get_df(self, *args, **kwargs):
        """Deprecated: Use to_df() instead."""
        import warnings
        warnings.warn(
            "get_df() is deprecated and will be removed in a future version. "
            "Use to_df() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.to_df(*args, **kwargs)

    def get_psm_df(self, export_all_hits=True, include_modifications=True, include_peak_annotations=False, decode_ontology=True, reference_file_name=""):
        """
        get_psm_df(self: PeptideIdentificationList, export_all_hits: bool = True, include_modifications: bool = True, include_peak_annotations: bool = False, decode_ontology: bool = True, reference_file_name: str = "") -> pd.DataFrame

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export PSMs as a DataFrame following the QPX PSM schema.

        Unlike get_df() which exports only the top hit, this method exports
        all hits (PSMs) with proper ranking information.

        :param export_all_hits: If True, export all hits per identification.
                                If False, only export rank 0 (best hit).
                                Default True.
        :type export_all_hits: bool

        :param include_modifications: Include detailed modification list column
                                      following QPX schema with name, accession,
                                      and positions array. Default True.
        :type include_modifications: bool

        :param include_peak_annotations: Include fragment ion annotations
                                         (mz_array, intensity_array, charge_array,
                                         ion_type_array). Default False as these
                                         can be large.
        :type include_peak_annotations: bool

        :param decode_ontology: Decode meta value names using the PSI-MS ontology.
                                Default True. (Currently unused, for future compatibility)
        :type decode_ontology: bool

        :param reference_file_name: Source file name for the QPX schema.
                                    Default empty string.
        :type reference_file_name: str

        :return: DataFrame with QPX PSM schema columns including:
                 sequence, peptidoform, modifications, precursor_charge,
                 posterior_error_probability, is_decoy, calculated_mz, observed_mz,
                 additional_scores, protein_accessions, predicted_rt, reference_file_name,
                 cv_params, scan, rt, ion_mobility, number_peaks, mz_array, intensity_array,
                 charge_array, ion_type_array, ion_mobility_array
                 Plus OpenMS-specific: spectrum_reference, score, score_type, rank, P_ID
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed

        Example::

            peps = feature_map.get_assigned_peptide_identifications()
            df = peps.get_psm_df()

            # Export only top hit per spectrum
            df = peps.get_psm_df(export_all_hits=False)

            # Without modification details
            df = peps.get_psm_df(include_modifications=False)
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for get_psm_df(). "
                "Please install it with: pip install pandas"
            )
        from . import SpectrumLookup as _SpectrumLookup
        from . import IDScoreSwitcherAlgorithm as _IDScoreSwitcherAlgorithm
        from . import ScoreType as _ScoreType

        # Native ID type accessions to try for scan number extraction
        # See SpectrumLookup.cpp for the full list of supported formats
        _native_id_accessions = [
            "MS:1000768",  # Thermo nativeID format (scan=)
            "MS:1000769",  # Waters nativeID format (scan=)
            "MS:1000771",  # Bruker/Agilent nativeID format (scan=)
            "MS:1001508",  # Agilent MassHunter nativeID format (scanId=)
            "MS:1000774",  # index= format
            "MS:1000777",  # spectrum= format
            "MS:1001530",  # plain number format
        ]

        def _extract_scan_number(native_id):
            """Extract scan number using OpenMS SpectrumLookup logic."""
            if not native_id:
                return None
            sl = _SpectrumLookup()
            for accession in _native_id_accessions:
                try:
                    scan_num = sl.extractScanNumber(native_id, accession)
                    if scan_num >= 0:
                        return scan_num
                except Exception:
                    continue
            return None

        rows = []
        # Create IDScoreSwitcherAlgorithm instance once for PEP detection
        idsa = _IDScoreSwitcherAlgorithm()

        # Known PEP score names from IDScoreSwitcherAlgorithm::type_to_str_
        # These are used to search metavalues when main score is not PEP
        _pep_metavalue_names = [
            b"Posterior Error Probability",
            b"pep",
            b"PEP",
            b"posterior_error_probability",
            b"MS:1001493",
            # Also check with _score suffix
            b"Posterior Error Probability_score",
            b"pep_score",
            b"PEP_score",
            b"posterior_error_probability_score",
        ]

        for pep_idx, pep_id in enumerate(self):
            rt = pep_id.getRT()
            observed_mz = pep_id.getMZ()
            # Spectrum reference is stored as metavalue with key "spectrum_reference"
            spec_ref = ""
            if pep_id.metaValueExists(b"spectrum_reference"):
                spec_ref = pep_id.getMetaValue(b"spectrum_reference")
                if isinstance(spec_ref, bytes):
                    spec_ref = spec_ref.decode("utf-8")
            score_type = pep_id.getScoreType()

            # Ion mobility from metavalue
            ion_mobility = None
            if pep_id.metaValueExists(b"ion_mobility"):
                ion_mobility = pep_id.getMetaValue(b"ion_mobility")
            elif pep_id.metaValueExists(b"IM"):
                ion_mobility = pep_id.getMetaValue(b"IM")

            # Extract scan number using OpenMS SpectrumLookup
            scan = _extract_scan_number(spec_ref)

            hits = pep_id.getHits()
            if not hits:
                continue  # Skip empty identifications

            num_hits = len(hits) if export_all_hits else min(1, len(hits))

            # Detect PEP score type once per PeptideIdentification (applies to all hits)
            # Use isScoreType to check if main score is PEP type
            is_main_score_pep = idsa.isScoreType(score_type, _ScoreType.PEP) if score_type else False

            for rank in range(num_hits):
                hit = hits[rank]
                seq = hit.getSequence()
                charge = hit.getCharge()

                # Extract modifications in QPX schema format
                # QPX format: [{"name", "accession", "positions": [{"position": "{AA}.{pos}", "scores"}]}]
                modifications = []
                if include_modifications and seq.isModified():
                    # Track modifications by name to combine same modifications at different positions
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
                        position_str = f"C-term.{seq.size() + 1}"
                        if mod_name not in mod_dict:
                            mod_dict[mod_name] = {"name": mod_name, "accession": accession_str, "positions": []}
                        mod_dict[mod_name]["positions"].append({"position": position_str, "scores": None})

                    modifications = list(mod_dict.values())

                # Determine is_decoy
                is_decoy = None
                if hit.metaValueExists(b"target_decoy"):
                    td = hit.getMetaValue(b"target_decoy")
                    if isinstance(td, bytes):
                        is_decoy = td.startswith(b"decoy")
                    else:
                        is_decoy = str(td).startswith("decoy")

                # Extract protein accessions
                evidences = hit.getPeptideEvidences()
                protein_accessions = [ev.getProteinAccession() for ev in evidences]

                # Build additional scores as array of records (QPX schema)
                # Format: [{"score_name", "score_value", "higher_better"}]
                # Exclude metavalues that are handled as dedicated columns
                _excluded_metavalues = {b"target_decoy", b"predicted_RT", b"predicted_rt"}
                additional_scores = []
                keys = []
                hit.getKeys(keys)
                for key in keys:
                    if key not in _excluded_metavalues:
                        val = hit.getMetaValue(key)
                        if isinstance(val, (int, float)):
                            key_str = key.decode("utf-8") if isinstance(key, bytes) else str(key)
                            # Determine higher_better using IDScoreSwitcherAlgorithm if possible
                            higher_better = None
                            try:
                                score_type_enum = _IDScoreSwitcherAlgorithm.toScoreTypeEnum(key_str)
                                higher_better = idsa.isScoreTypeHigherBetter(score_type_enum)
                            except Exception:
                                pass  # Unknown score type, leave as None
                            additional_scores.append({
                                "score_name": key_str,
                                "score_value": float(val),
                                "higher_better": higher_better
                            })

                # Calculate m/z from sequence
                calculated_mz = None
                if charge > 0:
                    try:
                        calculated_mz = seq.getMZ(charge)
                    except Exception:
                        pass

                # Get predicted RT if available
                predicted_rt = None
                if hit.metaValueExists(b"predicted_RT"):
                    predicted_rt = hit.getMetaValue(b"predicted_RT")
                elif hit.metaValueExists(b"predicted_rt"):
                    predicted_rt = hit.getMetaValue(b"predicted_rt")

                # Build posterior_error_probability using IDScoreSwitcherAlgorithm
                # is_main_score_pep was computed once per pep_id above using isScoreType()
                pep_value = None
                if is_main_score_pep:
                    pep_value = hit.getScore()
                else:
                    # Search metavalues for known PEP score names
                    for pep_name in _pep_metavalue_names:
                        if hit.metaValueExists(pep_name):
                            pep_value = hit.getMetaValue(pep_name)
                            break

                # Build row with QPX schema fields in order
                row = {
                    # QPX schema fields (in schema order)
                    "sequence": seq.toUnmodifiedString(),
                    "peptidoform": seq.toString(),
                    "modifications": modifications if include_modifications else None,
                    "precursor_charge": charge,
                    "posterior_error_probability": pep_value,
                    "is_decoy": 1 if is_decoy else 0 if is_decoy is not None else None,  # QPX uses int32
                    "calculated_mz": calculated_mz,
                    "observed_mz": observed_mz,
                    "additional_scores": additional_scores,
                    "protein_accessions": protein_accessions,
                    "predicted_rt": predicted_rt,
                    "reference_file_name": reference_file_name,  # QPX required field
                    "cv_params": None,  # QPX optional field (not populated from OpenMS)
                    "scan": str(scan) if scan is not None else None,  # QPX uses string
                    "rt": rt,
                    "ion_mobility": ion_mobility,
                    # OpenMS-specific fields (not in QPX)
                    "spectrum_reference": spec_ref,
                    "score": hit.getScore(),
                    "score_type": score_type,
                    "rank": rank,
                    "P_ID": pep_idx,
                }

                # Add fragment ion peak annotations (QPX schema fields)
                # Always include these columns for schema consistency with get_psm_columns()
                if include_peak_annotations:
                    peak_annotations = hit.getPeakAnnotations()
                    if peak_annotations:
                        row["number_peaks"] = len(peak_annotations)
                        row["mz_array"] = [ann.mz for ann in peak_annotations]
                        row["intensity_array"] = [ann.intensity for ann in peak_annotations]
                        row["charge_array"] = [ann.charge for ann in peak_annotations]
                        row["ion_type_array"] = [ann.annotation.decode("utf-8") if isinstance(ann.annotation, bytes) else ann.annotation for ann in peak_annotations]
                        row["ion_mobility_array"] = None  # QPX field, not available from OpenMS peak annotations
                    else:
                        row["number_peaks"] = 0
                        row["mz_array"] = []
                        row["intensity_array"] = []
                        row["charge_array"] = []
                        row["ion_type_array"] = []
                        row["ion_mobility_array"] = None
                else:
                    # Set to None when not requested, but always include for consistent schema
                    row["number_peaks"] = None
                    row["mz_array"] = None
                    row["intensity_array"] = None
                    row["charge_array"] = None
                    row["ion_type_array"] = None
                    row["ion_mobility_array"] = None

                rows.append(row)

        return pd.DataFrame(rows)

    def to_qpx(self, qpx_version="1.0", creator="pyOpenMS", software_provider="OpenMS",
                scan_format="scan", **kwargs):
        """
        to_qpx(self: PeptideIdentificationList, qpx_version: str = "1.0", creator: str = "pyOpenMS", software_provider: str = "OpenMS", scan_format: str = "scan", **kwargs) -> dict

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export PSMs as QPX format structure with file metadata and PSMs array.

        The QPX (Quantitative Proteomics Exchange) schema defines a standard
        format for PSM data exchange. This method exports data as a dict with:
        - file_metadata: File-level metadata (qpx_version, creator, etc.)
        - psms: List of PSM records

        QPX PSM Schema Fields:
            sequence, peptidoform, modifications, precursor_charge,
            posterior_error_probability, is_decoy, calculated_mz, observed_mz,
            additional_scores, protein_accessions, predicted_rt, reference_file_name,
            cv_params, scan, rt, ion_mobility, number_peaks, mz_array,
            intensity_array, charge_array, ion_type_array, ion_mobility_array

        OpenMS-specific Fields (appended to each PSM):
            spectrum_reference, score, score_type, rank, P_ID

        :param qpx_version: Version of the QPX format. Default "1.0".
        :type qpx_version: str

        :param creator: Name of the tool or person who created the file.
                        Default "pyOpenMS".
        :type creator: str

        :param software_provider: Name of the software provider. Default "OpenMS".
        :type software_provider: str

        :param scan_format: Format of scan identifiers: "scan", "index", or "nativeId".
                            Default "scan".
        :type scan_format: str

        Additional kwargs are passed to get_psm_df().

        :return: Dict with 'file_metadata' and 'psms' keys following QPX schema.
        :rtype: dict

        Example::

            qpx_data = peps.to_qpx(reference_file_name="sample.mzML")

            # Access file metadata
            print(qpx_data["file_metadata"]["qpx_version"])

            # Access PSMs as DataFrame
            import pandas as pd
            psms_df = pd.DataFrame(qpx_data["psms"])

            # Write to Parquet with pyarrow
            import pyarrow as pa
            import pyarrow.parquet as pq
            table = pa.Table.from_pylist(qpx_data["psms"])
            pq.write_table(table, "psms.parquet")
        """
        import uuid as _uuid
        from datetime import datetime as _datetime

        df = self.get_psm_df(**kwargs)

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

    def to_qpx_arrow(self, **kwargs):
        """
        to_qpx_arrow(self: PeptideIdentificationList, **kwargs) -> pa.Table

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export PSMs as Apache Arrow Table (PSMs only, without file metadata).

        This is a convenience method for direct Arrow export of PSM records.
        For full QPX format with file metadata, use to_qpx().

        Accepts same parameters as get_psm_df().

        :return: Arrow Table with PSM data.
        :rtype: pyarrow.Table

        :raises ImportError: If pyarrow is not installed

        Example::

            table = peps.to_qpx_arrow()

            # Convert to polars
            import polars as pl
            df = pl.from_arrow(table)

            # Convert to pandas
            df = table.to_pandas()
        """
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_qpx_arrow(). "
                "Please install it with: pip install pyarrow"
            )
        df = self.get_psm_df(**kwargs)
        return pa.Table.from_pandas(df)

    @staticmethod
    def get_psm_columns():
        """
        get_psm_columns() -> list

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Return list of column names that get_psm_df() and to_qpx() produce.

        The columns follow the QPX PSM schema order, with additional
        OpenMS-specific fields appended.

        Useful for discovering available columns before export.

        :return: List of column name strings.
        :rtype: list

        Example::

            >>> PeptideIdentificationList.get_psm_columns()
            ['sequence', 'peptidoform', 'modifications', ...]
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
            "rank",
            "P_ID",
        ]
