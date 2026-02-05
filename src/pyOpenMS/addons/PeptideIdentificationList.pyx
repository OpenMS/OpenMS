cimport numpy as np
import numpy as np
import warnings




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

    def to_df(self, decode_ontology=True, default_missing_values=None, export_unidentified=True, columns=None):
        """
        to_df(self: PeptideIdentificationList, decode_ontology: bool = True, default_missing_values: dict = None, export_unidentified: bool = True, columns: list = None) -> pd.DataFrame

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

        :param columns: List of column names to include. If None, includes all columns.
                        Use df_columns() to discover available columns.
        :type columns: list

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

            # Export only specific columns
            df = peps.to_df(columns=['id', 'rt', 'mz', 'charge'])
        """
        # Delegate to to_arrow() for DRY - single source of truth for data extraction
        table = self.to_arrow(
            decode_ontology=decode_ontology,
            default_missing_values=default_missing_values,
            export_unidentified=export_unidentified,
            columns=columns
        )
        return table.to_pandas()

    def df_columns(self):
        """
        df_columns(self: PeptideIdentificationList) -> list

        Returns a list of column names that to_df() would produce.

        Useful for discovering available columns before export, especially
        for column selection with the `columns` parameter.

        :return: List of column name strings.
        :rtype: list

        Example::

            peps = feature_map.get_assigned_peptide_identifications()
            print(peps.df_columns())

            # Use for column selection
            df = peps.to_df(columns=peps.df_columns()[:5])
        """
        return [
            "id",
            "rt",
            "mz",
            "score",
            "charge",
            "protein_accession",
            "start",
            "end",
            "P_ID",
            "PSM_ID",
            "metavalues",  # List of {"name": str, "value": str, "type": str}
        ]

    def to_arrow(self, decode_ontology=True, default_missing_values=None, export_unidentified=True, columns=None):
        """
        to_arrow(self: PeptideIdentificationList, decode_ontology: bool = True, default_missing_values: dict = None, export_unidentified: bool = True, columns: list = None) -> pa.Table

        Returns an Apache Arrow Table with peptide identification information.

        This method builds the Arrow table directly without pandas intermediate,
        making it more memory-efficient than to_df() for large datasets.

        :param decode_ontology: Decode meta value names using the PSI-MS ontology.
                                Default True.
        :type decode_ontology: bool

        :param default_missing_values: Default values for missing data by type.
                                       Default: {<bool>: False, <int>: -9999, <float>: nan, <str>: ''}
        :type default_missing_values: dict

        :param export_unidentified: Export PeptideIdentifications without PeptideHit.
                                    Default True.
        :type export_unidentified: bool

        :param columns: List of column names to include. If None, includes all columns.
                        Use df_columns() to discover available columns.
        :type columns: list

        :return: Arrow Table with peptide identification data.
        :rtype: pyarrow.Table

        :raises ImportError: If pyarrow is not installed

        Example::

            table = peps.to_arrow()
            df = table.to_pandas()

            # Convert to polars
            import polars as pl
            df = pl.from_arrow(table)

            # Export only specific columns
            table = peps.to_arrow(columns=['id', 'rt', 'mz', 'charge'])
        """
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_arrow(). "
                "Please install it with: pip install pyarrow"
            )
        from . import ControlledVocabulary as _ControlledVocabulary
        from . import File as _File

        if default_missing_values is None:
            default_missing_values = {type(True): False, type(1): -9999, type(1.0): np.nan, type(''): ''}

        # Load ontology for decoding metavalue names if requested
        cv = None
        if decode_ontology:
            cv = _ControlledVocabulary()
            cv.loadFromOBO("psims", _File.getOpenMSDataPath() + "/CV/psi-ms.obo")

        def decode_mv_name(name):
            """Decode metavalue name using ontology if available."""
            if cv is not None and name.startswith("MS:"):
                try:
                    return cv.getTerm(name).name
                except Exception:
                    return name
            return name

        def get_type_str(val):
            """Get type string for a value."""
            if type(val) is type(True):
                return "bool"
            elif type(val) is type(1):
                return "int"
            elif type(val) is type(1.0):
                return "float"
            else:
                return "str"

        # Initialize lists for core columns
        all_id = []
        all_rt = []
        all_mz = []
        all_score = []
        all_charge = []
        all_protein_accession = []
        all_start = []
        all_end = []
        all_p_id = []
        all_psm_id = []
        all_metavalues = []

        for pep_idx, pep in enumerate(self):
            hits = pep.getHits()

            if not hits:
                if export_unidentified:
                    all_id.append(pep.getIdentifier())
                    all_rt.append(pep.getRT())
                    all_mz.append(pep.getMZ())
                    all_score.append(default_missing_values[type(1.0)])
                    all_charge.append(default_missing_values[type(1)])
                    all_protein_accession.append(default_missing_values[type('')])
                    all_start.append(default_missing_values[type('')])
                    all_end.append(default_missing_values[type('')])
                    all_p_id.append(pep_idx)
                    all_psm_id.append(default_missing_values[type(1)])
                    all_metavalues.append([])
                continue

            besthit = hits[0]
            all_id.append(pep.getIdentifier())
            all_rt.append(pep.getRT())
            all_mz.append(pep.getMZ())
            all_score.append(besthit.getScore())
            all_charge.append(besthit.getCharge())

            # Add protein evidences as comma-separated strings
            evs = besthit.getPeptideEvidences()
            accessions = [e.getProteinAccession() for e in evs]
            starts = [str(e.getStart()) for e in evs]
            ends = [str(e.getEnd()) for e in evs]

            all_protein_accession.append(','.join(accessions) if accessions else default_missing_values[type('')])
            all_start.append(','.join(starts) if starts else default_missing_values[type('')])
            all_end.append(','.join(ends) if ends else default_missing_values[type('')])

            all_p_id.append(pep_idx)
            all_psm_id.append(0)  # Only first hit exported

            # Collect metavalues as key-value pairs
            metavalues = []
            mvs = []
            besthit.getKeys(mvs)
            for mv in mvs:
                if besthit.metaValueExists(mv):
                    val = besthit.getMetaValue(mv)
                    decoded_name = decode_mv_name(mv)
                    metavalues.append({
                        "name": decoded_name,
                        "value": str(val),
                        "value_type": get_type_str(val)
                    })
            all_metavalues.append(metavalues)

        # Build Arrow table directly using pa.Table.from_pydict
        data_dict = {}

        # Only include requested columns
        columns_set = set(columns) if columns is not None else None

        def should_include(col_name):
            return columns_set is None or col_name in columns_set

        # Core columns with proper Arrow types
        if should_include("id"):
            data_dict["id"] = pa.array(all_id, type=pa.utf8())
        if should_include("rt"):
            data_dict["rt"] = pa.array(all_rt, type=pa.float32())
        if should_include("mz"):
            data_dict["mz"] = pa.array(all_mz, type=pa.float32())
        if should_include("score"):
            data_dict["score"] = pa.array(all_score, type=pa.float64())
        if should_include("charge"):
            data_dict["charge"] = pa.array(all_charge, type=pa.int32())
        if should_include("protein_accession"):
            data_dict["protein_accession"] = pa.array(all_protein_accession, type=pa.utf8())
        if should_include("start"):
            data_dict["start"] = pa.array(all_start, type=pa.utf8())
        if should_include("end"):
            data_dict["end"] = pa.array(all_end, type=pa.utf8())
        if should_include("P_ID"):
            data_dict["P_ID"] = pa.array(all_p_id, type=pa.int32())
        if should_include("PSM_ID"):
            data_dict["PSM_ID"] = pa.array(all_psm_id, type=pa.int32())

        # Metavalues as list of structs (matches QPX schema)
        if should_include("metavalues"):
            metavalue_type = pa.struct([
                ("name", pa.utf8()),
                ("value", pa.utf8()),
                ("value_type", pa.utf8())
            ])
            data_dict["metavalues"] = pa.array(all_metavalues, type=pa.list_(metavalue_type))

        if columns_set is not None:
            unknown = [c for c in columns if c not in data_dict]
            if unknown:
                warnings.warn(f"Unknown column(s) ignored in to_arrow(): {unknown}. "
                              f"Use df_columns() to see available columns.")

        return pa.Table.from_pydict(data_dict)

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

    # TODO: Consider adding ontology decoding support (decode_ontology parameter) to PSM export methods
    # to decode CV term accessions (e.g., MS:1001493 -> "Posterior Error Probability")
    def to_psm_df(self, export_all_hits=True, include_modifications=True, include_peak_annotations=False, reference_file_name="", columns=None, additional_score_names=None, scan_format="scan"):
        """
        to_psm_df(self: PeptideIdentificationList, export_all_hits: bool = True, include_modifications: bool = True, include_peak_annotations: bool = False, reference_file_name: str = "", columns: list = None, additional_score_names: list = None, scan_format: str = "scan") -> pd.DataFrame

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export all PSMs (Peptide-Spectrum Matches) as a DataFrame.

        Unlike to_df() which exports only the top hit per spectrum, this method
        exports all hits with proper ranking information, modifications, and
        additional scores. Results are sorted by retention time (rt), precursor
        m/z (observed_mz), precursor charge, and rank for consistent ordering.

        :param export_all_hits: If True, export all hits per identification.
                                If False, only export rank 0 (best hit).
                                Default True.
        :type export_all_hits: bool

        :param include_modifications: Include detailed modification list column
                                      with name, accession, and positions array.
                                      Default True.
        :type include_modifications: bool

        :param include_peak_annotations: Include fragment ion annotations
                                         (mz_array, intensity_array, charge_array,
                                         ion_type_array). Default False as these
                                         can be large.
        :type include_peak_annotations: bool

        :param reference_file_name: Source file name to include in each row.
                                    Default empty string.
        :type reference_file_name: str

        :param columns: List of column names to include. If None, includes all columns.
                        Use psm_columns() to discover available columns.
        :type columns: list

        :param additional_score_names: List of additional metavalue names to treat as scores.
                                       By default, only known score types (from Scores.getAllIDScoreNames())
                                       are included in additional_scores. Non-score metavalues go to
                                       psm_metavalues. Use this parameter to specify custom score names
                                       that should be treated as scores and included in additional_scores.
        :type additional_score_names: list

        :param scan_format: Controls what the ``scan`` column contains:
                            ``"scan"`` (default) extracts the scan number from the
                            native ID, ``"nativeId"`` uses the raw native ID string.
                            The column is always string-typed (``pa.utf8()``) regardless
                            of the format, for consistency across both modes.
        :type scan_format: str

        :return: DataFrame with columns: sequence, peptidoform, modifications,
                 precursor_charge, posterior_error_probability, is_decoy,
                 calculated_mz, observed_mz, additional_scores, protein_accessions,
                 predicted_rt, reference_file_name, cv_params, scan, rt, ion_mobility,
                 spectrum_reference, score, score_type, rank, P_ID,
                 psm_metavalues, spectrum_metavalues.
                 When include_peak_annotations=True: number_peaks, mz_array,
                 intensity_array, charge_array, ion_type_array, ion_mobility_array.
                 Results are sorted by rt, observed_mz, precursor_charge, rank.
                 Empty input returns an empty DataFrame with proper column schema.
        :rtype: pd.DataFrame

        :raises ImportError: If pyarrow is not installed (required for efficient data extraction)

        Example::

            peps = feature_map.get_assigned_peptide_identifications()
            df = peps.to_psm_df()

            # Export only top hit per spectrum
            df = peps.to_psm_df(export_all_hits=False)

            # Without modification details
            df = peps.to_psm_df(include_modifications=False)

            # Export only specific columns
            df = peps.to_psm_df(columns=['sequence', 'precursor_charge', 'score'])

            # Include custom scores in additional_scores (instead of psm_metavalues)
            df = peps.to_psm_df(additional_score_names=['my_custom_score', 'rescoring_score'])
        """
        # Delegate to to_psm_arrow() for DRY - single source of truth for PSM extraction
        # The Arrow->pandas conversion is essentially zero-copy for most column types.
        try:
            table = self.to_psm_arrow(
                export_all_hits=export_all_hits,
                include_modifications=include_modifications,
                include_peak_annotations=include_peak_annotations,
                reference_file_name=reference_file_name,
                columns=columns,
                additional_score_names=additional_score_names,
                scan_format=scan_format
            )
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_psm_df(). "
                "Please install it with: pip install pyarrow"
            )
        return table.to_pandas()

    def to_psm_qpx(self, qpx_version="1.0", creator="pyOpenMS", software_provider="OpenMS",
                    scan_format="scan", **kwargs):
        """
        to_psm_qpx(self: PeptideIdentificationList, qpx_version: str = "1.0", creator: str = "pyOpenMS", software_provider: str = "OpenMS", scan_format: str = "scan", **kwargs) -> dict

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

        :param scan_format: Controls the ``scan`` column and is recorded in file_metadata.
                            ``"scan"`` (default) extracts the scan number from the
                            native ID, ``"nativeId"`` uses the raw native ID string.
                            The column is always string-typed for consistency.
        :type scan_format: str

        Additional kwargs are passed to to_psm_df().

        :return: Dict with 'file_metadata' and 'psms' keys following QPX schema.
        :rtype: dict

        Example::

            qpx_data = peps.to_psm_qpx(reference_file_name="sample.mzML")

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

    def to_psm_arrow(self, export_all_hits=True, include_modifications=True,
                      include_peak_annotations=False,
                      reference_file_name="", columns=None, additional_score_names=None,
                      scan_format="scan"):
        """
        to_psm_arrow(self: PeptideIdentificationList, export_all_hits: bool = True, include_modifications: bool = True, include_peak_annotations: bool = False, reference_file_name: str = "", columns: list = None, additional_score_names: list = None, scan_format: str = "scan") -> pa.Table

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export all PSMs as Apache Arrow Table.

        This method builds the Arrow table directly without pandas intermediate,
        making it more memory-efficient than to_psm_df() for large datasets.
        Complex nested types (modifications, additional_scores) are represented
        as proper Arrow struct and list types, enabling efficient Parquet storage.
        Results are sorted by retention time (rt), precursor m/z (observed_mz),
        precursor charge, and rank for consistent ordering.

        :param export_all_hits: If True, export all hits per identification.
                                If False, only export rank 0 (best hit).
                                Default True.
        :type export_all_hits: bool

        :param include_modifications: Include detailed modification list column
                                      with name, accession, and positions array.
                                      Default True.
        :type include_modifications: bool

        :param include_peak_annotations: Include fragment ion annotations
                                         (mz_array, intensity_array, charge_array,
                                         ion_type_array). Default False as these
                                         can be large.
        :type include_peak_annotations: bool

        :param reference_file_name: Source file name to include in each row.
                                    Default empty string.
        :type reference_file_name: str

        :param columns: List of column names to include. If None, includes all columns.
                        Use psm_columns() to discover available columns.
        :type columns: list

        :param additional_score_names: List of additional metavalue names to treat as scores.
                                       By default, only known score types (from Scores.getAllIDScoreNames())
                                       are included in additional_scores. Non-score metavalues go to
                                       psm_metavalues. Use this parameter to specify custom score names.
        :type additional_score_names: list

        :param scan_format: Controls what the ``scan`` column contains:
                            ``"scan"`` (default) extracts the scan number from the
                            native ID, ``"nativeId"`` uses the raw native ID string.
                            The column is always string-typed (``pa.utf8()``) regardless
                            of the format, for consistency across both modes.
        :type scan_format: str

        :return: Arrow Table with PSM data using native Arrow types.
                 Results are sorted by rt, observed_mz, precursor_charge, rank.
                 Empty input returns an empty table with proper schema.
        :rtype: pyarrow.Table

        :raises ImportError: If pyarrow is not installed

        Example::

            table = peps.to_psm_arrow()

            # Convert to polars
            import polars as pl
            df = pl.from_arrow(table)

            # Convert to pandas
            df = table.to_pandas()

            # Write to Parquet
            import pyarrow.parquet as pq
            pq.write_table(table, "psms.parquet")

            # Include custom scores in additional_scores
            table = peps.to_psm_arrow(additional_score_names=['my_custom_score'])
        """
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_psm_arrow(). "
                "Please install it with: pip install pyarrow"
            )

        if scan_format not in ("scan", "nativeId"):
            raise ValueError(f"scan_format must be 'scan' or 'nativeId', got '{scan_format}'")

        from . import SpectrumLookup as _SpectrumLookup
        from . import IDScoreSwitcherAlgorithm as _IDScoreSwitcherAlgorithm
        from . import Scores as _Scores
        from . import ProForma as _ProForma
        _IDType = _Scores.IDType

        # Native ID type accessions for scan number extraction
        _native_id_accessions = [
            "MS:1000768", "MS:1000769", "MS:1000771", "MS:1001508",
            "MS:1000774", "MS:1000777", "MS:1001530",
        ]

        # Create SpectrumLookup once for all scan extractions (not per-call)
        _spectrum_lookup = _SpectrumLookup()

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
            """Check if a name is a known score type (using C++ registry or user-specified)."""
            if name in _additional_score_names_set:
                return True
            # Use the C++ Scores.isKnownScoreType which handles _score suffix normalization
            return _Scores.isKnownScoreType(name)

        def _get_value_type(val):
            """Determine the type string for a metavalue."""
            # Check bool first using type name (bool is subclass of int in Python)
            # Note: isinstance(val, bool) doesn't work in Cython as 'bool' is C bint
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
        all_rank = []
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
        idsa = _IDScoreSwitcherAlgorithm()

        for pep_idx, pep_id in enumerate(self):
            rt = pep_id.getRT()
            observed_mz = pep_id.getMZ()

            # Spectrum reference
            spec_ref = ""
            if pep_id.metaValueExists("spectrum_reference"):
                spec_ref = pep_id.getMetaValue("spectrum_reference")

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
                if key not in _excluded_spectrum_metavalues:
                    val = pep_id.getMetaValue(key)
                    val_type = _get_value_type(val)
                    spectrum_metavalues.append({
                        "name": key,
                        "value": val,
                        "value_type": val_type
                    })

            hits = pep_id.getHits()
            if not hits:
                continue

            num_hits = len(hits) if export_all_hits else min(1, len(hits))

            # Detect PEP score type using findScoreType (searches main score and metavalues)
            pep_search_result = idsa.findScoreType(pep_id, _IDType.PEP)
            pep_score_name = pep_search_result.score_name  # Empty if not found

            for rank in range(num_hits):
                hit = hits[rank]
                seq = hit.getSequence()
                charge = hit.getCharge()

                # Extract modifications
                modifications = []
                if include_modifications and seq.isModified():
                    mod_dict = {}

                    def _add_mod(mod_dict, mod, mod_name, position_str):
                        accession = mod.getUniModRecordId() if mod else None
                        accession_str = f"UNIMOD:{accession}" if accession and accession > 0 else None
                        if mod_name not in mod_dict:
                            mod_dict[mod_name] = {"name": mod_name, "accession": accession_str, "positions": []}
                        mod_dict[mod_name]["positions"].append({"position": position_str, "scores": None})

                    if seq.hasNTerminalModification():
                        _add_mod(mod_dict, seq.getNTerminalModification(),
                                 seq.getNTerminalModificationName(), "N-term.0")

                    for pos, residue in enumerate(seq):
                        if residue.isModified():
                            _add_mod(mod_dict, residue.getModification(),
                                     residue.getModificationName(),
                                     f"{residue.getOneLetterCode()}.{pos + 1}")

                    if seq.hasCTerminalModification():
                        _add_mod(mod_dict, seq.getCTerminalModification(),
                                 seq.getCTerminalModificationName(),
                                 f"C-term.{seq.size() + 1}")

                    modifications = list(mod_dict.values())

                # Determine is_decoy
                is_decoy = None
                if hit.metaValueExists("target_decoy"):
                    td = hit.getMetaValue("target_decoy")
                    is_decoy = 1 if str(td).startswith("decoy") else 0

                # Protein accessions
                evidences = hit.getPeptideEvidences()
                protein_accessions = [ev.getProteinAccession() for ev in evidences]

                # Build additional_scores and psm_metavalues from PeptideHit metavalues
                # additional_scores: numeric values where name is a known score type
                # psm_metavalues: all other metavalues (non-scores)
                additional_scores = []
                psm_metavalues = []
                keys = []
                hit.getKeys(keys)
                for key in keys:
                    if key not in _excluded_psm_metavalues:
                        val = hit.getMetaValue(key)

                        # Check if this is a known score type (using C++ registry)
                        is_known_score = _is_known_score(key)

                        if isinstance(val, (int, float)) and is_known_score:
                            # Known score - add to additional_scores
                            higher_better = None
                            try:
                                score_type_enum = _IDScoreSwitcherAlgorithm.toScoreTypeEnum(key)
                                higher_better = idsa.isScoreTypeHigherBetter(score_type_enum)
                            except Exception:
                                pass  # Unknown score type, leave as None
                            additional_scores.append({
                                "score_name": key,
                                "score_value": float(val),
                                "higher_better": higher_better
                            })
                        else:
                            # Not a known score - add to psm_metavalues
                            val_type = _get_value_type(val)
                            psm_metavalues.append({
                                "name": key,
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

                # PEP value using findScoreType result
                pep_value = None
                if pep_score_name:  # Non-empty means PEP was found
                    if pep_search_result.is_main_score_type:
                        pep_value = hit.getScore()
                    elif hit.metaValueExists(pep_score_name):
                        pep_value = hit.getMetaValue(pep_score_name)

                # Append to column lists
                all_sequence.append(seq.toUnmodifiedString())
                pf = _ProForma.fromAASequence(seq)
                all_peptidoform.append(_ProForma.toString(pf, _ProForma.WriteMode.CANONICAL))
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
                all_rank.append(rank)
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
            ("scores", pa.float64())  # nullable
        ])
        modification_type = pa.struct([
            ("name", pa.utf8()),
            ("accession", pa.utf8()),
            ("positions", pa.list_(position_type))
        ])
        score_type_struct = pa.struct([
            ("score_name", pa.utf8()),
            ("score_value", pa.float64()),
            ("higher_better", pa.bool_())  # nullable
        ])
        # Metavalue type for psm_metavalues and spectrum_metavalues
        # Uses utf8 for value since metavalues can be heterogeneous (int, float, str)
        # value_type field enables consumers to reconstruct typed values
        metavalue_type = pa.struct([
            ("name", pa.utf8()),
            ("value", pa.utf8()),      # stringified values for heterogeneous data
            ("value_type", pa.utf8())  # original type: "int", "float", "str", "bool"
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
        if should_include("rank"):
            data_dict["rank"] = pa.array(all_rank, type=pa.int32())
        if should_include("P_ID"):
            data_dict["P_ID"] = pa.array(all_p_id, type=pa.int32())

        # Metavalue columns - stringify values for Arrow compatibility
        if should_include("psm_metavalues"):
            # Convert values to strings for Arrow struct, preserving type info
            all_psm_metavalues_str = [
                [{"name": mv["name"],
                  "value": str(mv["value"]) if mv["value"] is not None else None,
                  "value_type": mv["value_type"]}
                 for mv in mvs] for mvs in all_psm_metavalues
            ]
            data_dict["psm_metavalues"] = pa.array(all_psm_metavalues_str, type=pa.list_(metavalue_type))
        if should_include("spectrum_metavalues"):
            # Convert values to strings for Arrow struct, preserving type info
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

        # Sort by rt, observed_mz, precursor_charge, rank for consistent ordering
        sort_cols = ["rt", "observed_mz", "precursor_charge", "rank"]
        available_sort_cols = [(c, "ascending") for c in sort_cols if c in table.column_names]
        if available_sort_cols and table.num_rows > 0:
            import pyarrow.compute as pc
            indices = pc.sort_indices(table, sort_keys=available_sort_cols, null_placement="at_end")
            table = table.take(indices)

        return table

    def to_parquet(self, path, compression='zstd', compression_level=None,
                    row_group_size=None, write_statistics=True,
                    export_all_hits=True, include_modifications=True,
                    include_peak_annotations=False,
                    reference_file_name="", columns=None, additional_score_names=None,
                    scan_format="scan"):
        """
        to_parquet(self: PeptideIdentificationList, path: str, compression: str = 'zstd', compression_level: int = None, row_group_size: int = None, write_statistics: bool = True, export_all_hits: bool = True, include_modifications: bool = True, include_peak_annotations: bool = False, reference_file_name: str = "", columns: list = None, additional_score_names: list = None, scan_format: str = "scan") -> None

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export all PSMs to a Parquet file using native Arrow construction.

        This method exports PSM data directly to Parquet format without pandas
        intermediate, making it efficient for large datasets. Results are sorted
        by retention time (rt), precursor m/z (observed_mz), precursor charge,
        and rank for consistent ordering. Parquet provides:
        - Columnar storage optimized for analytical queries
        - Efficient compression (typically 3-5x with ZSTD)
        - Fast partial reads (only requested columns/row groups are loaded)
        - Wide ecosystem support (DuckDB, Polars, pandas, Spark, R arrow)

        :param path: Path to the output Parquet file.
        :type path: str

        :param compression: Compression algorithm. Options:
                            'none', 'snappy', 'gzip', 'lz4', 'zstd' (default).
                            ZSTD offers the best ratio/speed tradeoff for MS data.
        :type compression: str

        :param compression_level: Compression level (interpretation depends on algorithm).
                                  ZSTD: 1-22 (default None uses pyarrow default).
                                  GZIP: 1-9. Ignored for SNAPPY/LZ4/NONE.
        :type compression_level: int

        :param row_group_size: Number of rows per row group. Default None writes
                               all rows to a single row group. Smaller values
                               enable more parallelism for readers.
        :type row_group_size: int

        :param write_statistics: Write column statistics (min/max/null_count) for
                                 each row group. Enables predicate pushdown for
                                 efficient filtering. Default True.
        :type write_statistics: bool

        :param export_all_hits: If True, export all hits per identification.
                                If False, only export rank 0 (best hit). Default True.
        :type export_all_hits: bool

        :param include_modifications: Include modification details. Default True.
        :type include_modifications: bool

        :param include_peak_annotations: Include fragment ion annotations. Default False.
        :type include_peak_annotations: bool

        :param reference_file_name: Source file name to include. Default empty.
        :type reference_file_name: str

        :param columns: Column names to include. None includes all.
        :type columns: list

        :param additional_score_names: List of additional metavalue names to treat as scores.
                                       See to_psm_df() for details.
        :type additional_score_names: list

        :param scan_format: Controls what the ``scan`` column contains:
                            ``"scan"`` (default) extracts the scan number from the
                            native ID, ``"nativeId"`` uses the raw native ID string.
                            The column is always string-typed (``pa.utf8()``) regardless
                            of the format, for consistency across both modes.
        :type scan_format: str

        :raises ImportError: If pyarrow is not installed

        Example::

            # Basic export with default ZSTD compression
            peps.to_parquet("psms.parquet")

            # Maximum compression (ZSTD level 9)
            peps.to_parquet("psms.parquet", compression='zstd', compression_level=9)

            # Fast compression for temporary files
            peps.to_parquet("psms.parquet", compression='lz4')

            # With peak annotations
            peps.to_parquet("psms.parquet", include_peak_annotations=True)

            # Smaller row groups for parallel reads
            peps.to_parquet("psms.parquet", row_group_size=100000)

            # Read back with pandas
            import pandas as pd
            df = pd.read_parquet("psms.parquet")

            # Read back with polars
            import polars as pl
            df = pl.read_parquet("psms.parquet")

            # Read back with DuckDB (SQL queries)
            import duckdb
            result = duckdb.query("SELECT * FROM 'psms.parquet' WHERE score > 0.95")
        """
        try:
            import pyarrow.parquet as pq
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_parquet(). "
                "Please install it with: pip install pyarrow"
            )

        # Build the Arrow table using native construction
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

        # Add compression level if specified (only for zstd and gzip)
        if compression_level is not None and pq_compression in ('zstd', 'gzip'):
            write_kwargs['compression_level'] = compression_level

        # Add row group size if specified
        if row_group_size is not None:
            write_kwargs['row_group_size'] = row_group_size

        # Write Parquet - use file handle to avoid pyarrow LocalFileSystem
        # conflict with OpenMS Arrow C++ on macOS (double 'file' scheme registration)
        with open(path, 'wb') as f:
            pq.write_table(table, f, **write_kwargs)

    def psm_columns(self):
        """
        psm_columns(self: PeptideIdentificationList) -> list

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Return list of column names that to_psm_df() and to_psm_qpx() produce.

        Useful for discovering available columns before export, especially
        for column selection with the `columns` parameter.

        :return: List of column name strings.
        :rtype: list

        Example::

            peps = feature_map.get_assigned_peptide_identifications()
            print(peps.psm_columns())

            # Use for column selection
            df = peps.to_psm_df(columns=peps.psm_columns()[:5])
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
            # Metavalues (non-score) - format: [{"name": str, "value": any}]
            "psm_metavalues",  # PeptideHit-level metavalues (PSM-specific)
            "spectrum_metavalues",  # PeptideIdentification-level metavalues (spectrum-wide)
        ]
