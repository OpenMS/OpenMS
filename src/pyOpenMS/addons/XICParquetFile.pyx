from libcpp.vector cimport vector as libcpp_vector
from libc.stdint cimport int64_t
from cython.operator cimport dereference as deref
from String cimport String
from StringList cimport StringList
from XICParquetFile cimport XICChromatogram as _XICChromatogram
from XICParquetFile cimport XICAnalyte as _XICAnalyte
from XICParquetFile cimport XICRunInfo as _XICRunInfo
from pyopenms._parquet_filter cimport ParquetFilterBuilder as _PyParquetFilter
from ._pyopenms_1 cimport convString, convOutputString
import numpy as np


    def __init__(self, filename):
        """
        __init__(self, filename: Union[str, bytes, String, List[Union[str, bytes, String]]]) -> None

        Initialize from a single filename or a list/tuple of filenames.

        :param filename: Path to a .xic file or a list/tuple of paths.
        :type filename: Union[str, bytes, String, List[Union[str, bytes, String]]]

        Example::

            xic = XICParquetFile("file.xic")
            xic_multi = XICParquetFile(["file1.xic", "file2.xic"])
        """
        filenames = None
        if isinstance(filename, (str, bytes, String)):
            self._init_0(filename)
            return
        if isinstance(filename, (list, tuple)):
            filenames = []
            for entry in filename:
                if isinstance(entry, bytes):
                    filenames.append(entry)
                elif isinstance(entry, str):
                    filenames.append(entry.encode())
                elif isinstance(entry, String):
                    filenames.append(str(entry).encode())
                else:
                    raise TypeError("filenames must be str/bytes/String")
            self._init_1(filenames)
            return
        raise TypeError("filename must be str/bytes/String or a list/tuple of those")

    def __len__(self):
        """
        __len__(self: XICParquetFile) -> int

        Return the number of files associated with this instance.
        """
        try:
            return len(self.getFilenames())
        except Exception:
            return 0

    def __repr__(self):
        """
        __repr__(self: XICParquetFile) -> str

        Return a detailed string representation.
        """
        try:
            files = self.getFilenames()
            files = [f.decode() if isinstance(f, (bytes, bytearray)) else str(f) for f in files]
            return "XICParquetFile(n_files=%d, files=%r)" % (len(files), files)
        except Exception:
            return "XICParquetFile(<unavailable>)"

    def __str__(self):
        """
        __str__(self: XICParquetFile) -> str
        
        Return a concise string representation.
        """
        try:
            return "XICParquetFile(n_files=%d)" % len(self.getFilenames())
        except Exception:
            return "XICParquetFile(n_files=?)"

    def get_data_dict(self, explode=False, precursor_id=-1, transition_id=-1, modified_sequence="", precursor_charge=-1,
                      product_charge=-1, ms_level=-1, run_id=-1, filter=""):
        """
        get_data_dict(self, explode: bool = False, precursor_id: int = -1, transition_id: int = -1,
                      modified_sequence: str = "", precursor_charge: int = -1, product_charge: int = -1,
                      ms_level: int = -1, run_id: int = -1, filter: str = "") -> Dict[str, np.ndarray]

        Return chromatogram data as a dict of numpy arrays.

        If explode=True, returns long format with rt/intensity rows.
        Otherwise, rt and intensity are stored as object arrays of lists.

        :param explode: If True, return long format with one row per RT/intensity.
        :type explode: bool
        :param precursor_id: Optional precursor id (-1 to ignore).
        :type precursor_id: int
        :param transition_id: Optional transition id (-1 to ignore).
        :type transition_id: int
        :param modified_sequence: Optional modified sequence filter (empty to ignore).
        :type modified_sequence: str
        :param precursor_charge: Optional precursor charge filter (-1 to ignore).
        :type precursor_charge: int
        :param product_charge: Optional product charge filter (-1 to ignore).
        :type product_charge: int
        :param ms_level: Optional MS level filter (-1 to ignore).
        :type ms_level: int
        :param run_id: Optional run id filter (-1 to ignore).
        :type run_id: int
        :param filter: Optional filter expression string.
        :type filter: str

        :return: Dict of numpy arrays keyed by column name.
        :rtype: Dict[str, np.ndarray]

        Example::

            data = xic.get_data_dict(precursor_id=1318, explode=True)
        """
        cdef libcpp_vector[_XICChromatogram] chroms
        self.inst.get().getChromatograms(chroms,
                                         <int64_t>precursor_id,
                                         <int64_t>transition_id,
                                         deref((convString(modified_sequence)).get()),
                                         <int64_t>precursor_charge,
                                         <int64_t>product_charge,
                                         <int64_t>ms_level,
                                         <int64_t>run_id,
                                         deref((convString(filter)).get()))

        cdef list run_id_list = []
        cdef list source_file_list = []
        cdef list ms_level_list = []
        cdef list precursor_id_list = []
        cdef list transition_id_list = []
        cdef list modified_sequence_list = []
        cdef list precursor_charge_list = []
        cdef list product_charge_list = []
        cdef list detecting_transition_list = []
        cdef list precursor_decoy_list = []
        cdef list product_decoy_list = []
        cdef list transition_ordinal_list = []
        cdef list transition_type_list = []
        cdef list annotation_list = []
        cdef list rt_list = []
        cdef list intensity_list = []

        cdef size_t i, j
        cdef _XICChromatogram chrom
        for i in range(chroms.size()):
            chrom = chroms[i]
            if explode:
                if chrom.rt.size() == 0:
                    continue
                for j in range(chrom.rt.size()):
                    run_id_list.append(chrom.run_id)
                    source_file_list.append(convOutputString(chrom.source_file))
                    ms_level_list.append(chrom.ms_level)
                    precursor_id_list.append(chrom.precursor_id if chrom.has_precursor_id else None)
                    transition_id_list.append(chrom.transition_id if chrom.has_transition_id else None)
                    modified_sequence_list.append(convOutputString(chrom.modified_sequence))
                    precursor_charge_list.append(chrom.precursor_charge if chrom.has_precursor_charge else None)
                    product_charge_list.append(chrom.product_charge if chrom.has_product_charge else None)
                    detecting_transition_list.append(chrom.detecting_transition if chrom.has_detecting_transition else None)
                    precursor_decoy_list.append(chrom.precursor_decoy if chrom.has_precursor_decoy else None)
                    product_decoy_list.append(chrom.product_decoy if chrom.has_product_decoy else None)
                    transition_ordinal_list.append(chrom.transition_ordinal if chrom.has_transition_ordinal else None)
                    transition_type_list.append(convOutputString(chrom.transition_type))
                    annotation_list.append(convOutputString(chrom.annotation))
                    rt_list.append(chrom.rt[j])
                    intensity_list.append(chrom.intensity[j])
            else:
                run_id_list.append(chrom.run_id)
                source_file_list.append(convOutputString(chrom.source_file))
                ms_level_list.append(chrom.ms_level)
                precursor_id_list.append(chrom.precursor_id if chrom.has_precursor_id else None)
                transition_id_list.append(chrom.transition_id if chrom.has_transition_id else None)
                modified_sequence_list.append(convOutputString(chrom.modified_sequence))
                precursor_charge_list.append(chrom.precursor_charge if chrom.has_precursor_charge else None)
                product_charge_list.append(chrom.product_charge if chrom.has_product_charge else None)
                detecting_transition_list.append(chrom.detecting_transition if chrom.has_detecting_transition else None)
                precursor_decoy_list.append(chrom.precursor_decoy if chrom.has_precursor_decoy else None)
                product_decoy_list.append(chrom.product_decoy if chrom.has_product_decoy else None)
                transition_ordinal_list.append(chrom.transition_ordinal if chrom.has_transition_ordinal else None)
                transition_type_list.append(convOutputString(chrom.transition_type))
                annotation_list.append(convOutputString(chrom.annotation))
                rt_list.append([chrom.rt[k] for k in range(chrom.rt.size())])
                intensity_list.append([chrom.intensity[k] for k in range(chrom.intensity.size())])

        return {
            "run_id": np.array(run_id_list, dtype=np.int64),
            "source_file": np.array(source_file_list, dtype=object),
            "ms_level": np.array(ms_level_list, dtype=np.int64),
            "precursor_id": np.array(precursor_id_list, dtype=object),
            "transition_id": np.array(transition_id_list, dtype=object),
            "modified_sequence": np.array(modified_sequence_list, dtype=object),
            "precursor_charge": np.array(precursor_charge_list, dtype=object),
            "product_charge": np.array(product_charge_list, dtype=object),
            "detecting_transition": np.array(detecting_transition_list, dtype=object),
            "precursor_decoy": np.array(precursor_decoy_list, dtype=object),
            "product_decoy": np.array(product_decoy_list, dtype=object),
            "transition_ordinal": np.array(transition_ordinal_list, dtype=object),
            "transition_type": np.array(transition_type_list, dtype=object),
            "annotation": np.array(annotation_list, dtype=object),
            "rt": np.array(rt_list, dtype=object),
            "intensity": np.array(intensity_list, dtype=object),
        }

    def _get_data_dict_from_parquet_filter(self, parquet_filter, explode=False):
        """
        Internal helper to fetch chromatograms using a typed ParquetFilter.
        """
        if parquet_filter is None:
            return self.get_data_dict(explode=explode)
        cdef libcpp_vector[_XICChromatogram] chroms
        try:
            self.inst.get().getChromatograms(chroms, deref((<_PyParquetFilter>parquet_filter).inst.get()))
        except AttributeError:
            raise TypeError("parquet_filter must be a ParquetFilterBuilder instance")

        cdef list run_id_list = []
        cdef list source_file_list = []
        cdef list ms_level_list = []
        cdef list precursor_id_list = []
        cdef list transition_id_list = []
        cdef list modified_sequence_list = []
        cdef list precursor_charge_list = []
        cdef list product_charge_list = []
        cdef list detecting_transition_list = []
        cdef list precursor_decoy_list = []
        cdef list product_decoy_list = []
        cdef list transition_ordinal_list = []
        cdef list transition_type_list = []
        cdef list annotation_list = []
        cdef list rt_list = []
        cdef list intensity_list = []

        cdef size_t i, j
        cdef _XICChromatogram chrom
        for i in range(chroms.size()):
            chrom = chroms[i]
            if explode:
                if chrom.rt.size() == 0:
                    continue
                for j in range(chrom.rt.size()):
                    run_id_list.append(chrom.run_id)
                    source_file_list.append(convOutputString(chrom.source_file))
                    ms_level_list.append(chrom.ms_level)
                    precursor_id_list.append(chrom.precursor_id if chrom.has_precursor_id else None)
                    transition_id_list.append(chrom.transition_id if chrom.has_transition_id else None)
                    modified_sequence_list.append(convOutputString(chrom.modified_sequence))
                    precursor_charge_list.append(chrom.precursor_charge if chrom.has_precursor_charge else None)
                    product_charge_list.append(chrom.product_charge if chrom.has_product_charge else None)
                    detecting_transition_list.append(chrom.detecting_transition if chrom.has_detecting_transition else None)
                    precursor_decoy_list.append(chrom.precursor_decoy if chrom.has_precursor_decoy else None)
                    product_decoy_list.append(chrom.product_decoy if chrom.has_product_decoy else None)
                    transition_ordinal_list.append(chrom.transition_ordinal if chrom.has_transition_ordinal else None)
                    transition_type_list.append(convOutputString(chrom.transition_type))
                    annotation_list.append(convOutputString(chrom.annotation))
                    rt_list.append(chrom.rt[j])
                    intensity_list.append(chrom.intensity[j])
            else:
                run_id_list.append(chrom.run_id)
                source_file_list.append(convOutputString(chrom.source_file))
                ms_level_list.append(chrom.ms_level)
                precursor_id_list.append(chrom.precursor_id if chrom.has_precursor_id else None)
                transition_id_list.append(chrom.transition_id if chrom.has_transition_id else None)
                modified_sequence_list.append(convOutputString(chrom.modified_sequence))
                precursor_charge_list.append(chrom.precursor_charge if chrom.has_precursor_charge else None)
                product_charge_list.append(chrom.product_charge if chrom.has_product_charge else None)
                detecting_transition_list.append(chrom.detecting_transition if chrom.has_detecting_transition else None)
                precursor_decoy_list.append(chrom.precursor_decoy if chrom.has_precursor_decoy else None)
                product_decoy_list.append(chrom.product_decoy if chrom.has_product_decoy else None)
                transition_ordinal_list.append(chrom.transition_ordinal if chrom.has_transition_ordinal else None)
                transition_type_list.append(convOutputString(chrom.transition_type))
                annotation_list.append(convOutputString(chrom.annotation))
                rt_list.append([chrom.rt[k] for k in range(chrom.rt.size())])
                intensity_list.append([chrom.intensity[k] for k in range(chrom.intensity.size())])

        return {
            "run_id": np.array(run_id_list, dtype=np.int64),
            "source_file": np.array(source_file_list, dtype=object),
            "ms_level": np.array(ms_level_list, dtype=np.int64),
            "precursor_id": np.array(precursor_id_list, dtype=object),
            "transition_id": np.array(transition_id_list, dtype=object),
            "modified_sequence": np.array(modified_sequence_list, dtype=object),
            "precursor_charge": np.array(precursor_charge_list, dtype=object),
            "product_charge": np.array(product_charge_list, dtype=object),
            "detecting_transition": np.array(detecting_transition_list, dtype=object),
            "precursor_decoy": np.array(precursor_decoy_list, dtype=object),
            "product_decoy": np.array(product_decoy_list, dtype=object),
            "transition_ordinal": np.array(transition_ordinal_list, dtype=object),
            "transition_type": np.array(transition_type_list, dtype=object),
            "annotation": np.array(annotation_list, dtype=object),
            "rt": np.array(rt_list, dtype=object),
            "intensity": np.array(intensity_list, dtype=object),
        }

    def to_df(self, explode=False):
        """
        to_df(self, explode: bool = False) -> pandas.DataFrame

        Return chromatogram data as a pandas DataFrame.

        If explode=True, returns long format with rt/intensity rows.

        :param explode: If True, return long format with one row per RT/intensity.
        :type explode: bool

        :return: DataFrame with chromatogram data.
        :rtype: pandas.DataFrame

        :raises ImportError: If pandas is not installed

        Example::

            df = xic.to_df()
            df_long = xic.to_df(explode=True)
        """
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError("pandas is required for this method. Install with `pip install pandas`.") from e
        data = self.get_data_dict(explode=explode)
        if explode:
            return pd.DataFrame(data)
        # Convert to plain Python lists to avoid pandas inferring 2D arrays.
        clean = {}
        for k, v in data.items():
            if hasattr(v, "ndim") and v.ndim > 1:
                clean[k] = [row for row in v]
            else:
                clean[k] = list(v)
        return pd.DataFrame(clean)

    def get_analyte_dict(self, nest_transitions=True, columns=None):
        """
        get_analyte_dict(self, nest_transitions: bool = True,
                         columns: Optional[List[str]] = None) -> Dict[str, np.ndarray]

        Return unique analyte metadata as a dict of numpy arrays.

        If nest_transitions=False, each row represents a unique precursor-transition
        pair with scalar transition fields. If nest_transitions=True, each row
        represents a unique precursor and transition fields are lists.

        :param nest_transitions: Aggregate transition fields per precursor.
        :type nest_transitions: bool
        :param columns: Optional list of column names to return. If None,
                        uses the default analyte columns.
        :type columns: Optional[List[str]]

        :return: Dict of numpy arrays keyed by column name.
        :rtype: Dict[str, np.ndarray]

        Example::

            analytes = xic.get_analyte_dict()
            analytes_nested = xic.get_analyte_dict(nest_transitions=True)
            analytes_subset = xic.get_analyte_dict(columns=["precursor_id", "modified_sequence"])
        """
        cdef libcpp_vector[_XICAnalyte] analytes
        cdef StringList cols_cpp
        cdef list cols_norm = []
        cdef object s
        if columns is not None:
            if not isinstance(columns, (list, tuple)):
                raise TypeError("columns must be a list/tuple of str/bytes/String")
            for entry in columns:
                if isinstance(entry, bytes):
                    s = entry.decode()
                    cols_norm.append(s.lower())
                elif isinstance(entry, str):
                    s = entry
                    cols_norm.append(s.lower())
                elif isinstance(entry, String):
                    s = str(entry)
                    cols_norm.append(s.lower())
                else:
                    raise TypeError("columns must be str/bytes/String")
                cols_cpp.push_back(deref((convString(s)).get()))
        self.inst.get().getAnalytes(analytes, cols_cpp, <bint>nest_transitions)

        cdef list precursor_id_list = []
        cdef list modified_sequence_list = []
        cdef list precursor_charge_list = []
        cdef list precursor_decoy_list = []
        cdef list transition_id_list = []
        cdef list product_charge_list = []
        cdef list transition_ordinal_list = []
        cdef list detecting_transition_list = []
        cdef list product_decoy_list = []
        cdef list transition_type_list = []
        cdef list annotation_list = []

        cdef size_t i, j
        cdef list t_ids = []
        cdef list p_charges = []
        cdef list t_ordinals = []
        cdef list d_transitions = []
        cdef list p_decoys = []
        cdef list t_types = []
        cdef list annots = []
        cdef _XICAnalyte analyte
        for i in range(analytes.size()):
            analyte = analytes[i]
            precursor_id_list.append(analyte.precursor_id if analyte.has_precursor_id else None)
            modified_sequence_list.append(convOutputString(analyte.modified_sequence))
            precursor_charge_list.append(analyte.precursor_charge if analyte.has_precursor_charge else None)
            precursor_decoy_list.append(analyte.precursor_decoy if analyte.has_precursor_decoy else None)

            if not nest_transitions:
                transition_id_list.append(analyte.transition_id if analyte.has_transition_id else None)
                product_charge_list.append(analyte.product_charge if analyte.has_product_charge else None)
                transition_ordinal_list.append(analyte.transition_ordinal if analyte.has_transition_ordinal else None)
                detecting_transition_list.append(analyte.detecting_transition if analyte.has_detecting_transition else None)
                product_decoy_list.append(analyte.product_decoy if analyte.has_product_decoy else None)
                transition_type_list.append(convOutputString(analyte.transition_type))
                annotation_list.append(convOutputString(analyte.annotation))
            else:
                t_ids = []
                p_charges = []
                t_ordinals = []
                d_transitions = []
                p_decoys = []
                t_types = []
                annots = []
                for j in range(analyte.transition_ids.size()):
                    t_ids.append(analyte.transition_ids[j] if analyte.transition_ids[j] >= 0 else None)
                for j in range(analyte.product_charges.size()):
                    p_charges.append(analyte.product_charges[j] if analyte.product_charges[j] >= 0 else None)
                for j in range(analyte.transition_ordinals.size()):
                    t_ordinals.append(analyte.transition_ordinals[j] if analyte.transition_ordinals[j] >= 0 else None)
                for j in range(analyte.detecting_transitions.size()):
                    d_transitions.append(analyte.detecting_transitions[j] if analyte.detecting_transitions[j] >= 0 else None)
                for j in range(analyte.product_decoys.size()):
                    p_decoys.append(analyte.product_decoys[j] if analyte.product_decoys[j] >= 0 else None)
                for j in range(analyte.transition_types.size()):
                    t_types.append(convOutputString(analyte.transition_types[j]))
                for j in range(analyte.annotations.size()):
                    annots.append(convOutputString(analyte.annotations[j]))

                transition_id_list.append(t_ids)
                product_charge_list.append(p_charges)
                transition_ordinal_list.append(t_ordinals)
                detecting_transition_list.append(d_transitions)
                product_decoy_list.append(p_decoys)
                transition_type_list.append(t_types)
                annotation_list.append(annots)

        data = {
            "precursor_id": np.array(precursor_id_list, dtype=object),
            "modified_sequence": np.array(modified_sequence_list, dtype=object),
            "precursor_charge": np.array(precursor_charge_list, dtype=object),
            "precursor_decoy": np.array(precursor_decoy_list, dtype=object),
            "transition_id": np.array(transition_id_list, dtype=object),
            "product_charge": np.array(product_charge_list, dtype=object),
            "transition_ordinal": np.array(transition_ordinal_list, dtype=object),
            "detecting_transition": np.array(detecting_transition_list, dtype=object),
            "product_decoy": np.array(product_decoy_list, dtype=object),
            "transition_type": np.array(transition_type_list, dtype=object),
            "annotation": np.array(annotation_list, dtype=object),
        }
        if columns is not None:
            return {k: v for k, v in data.items() if k in cols_norm}
        return data

    def get_analyte_df(self, nest_transitions=True, columns=None):
        """
        get_analyte_df(self, nest_transitions: bool = True,
                       columns: Optional[List[str]] = None) -> pandas.DataFrame

        Return unique analyte metadata as a pandas DataFrame.

        :param nest_transitions: Aggregate transition fields per precursor.
        :type nest_transitions: bool
        :param columns: Optional list of column names to return.
        :type columns: Optional[List[str]]

        :return: DataFrame with analyte metadata.
        :rtype: pandas.DataFrame

        :raises ImportError: If pandas is not installed

        Example::

            df = xic.get_analyte_df()
            df_nested = xic.get_analyte_df(nest_transitions=True)
            df_subset = xic.get_analyte_df(columns=["precursor_id", "modified_sequence"])
        """
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError("pandas is required for this method. Install with `pip install pandas`.") from e
        data = self.get_analyte_dict(nest_transitions=nest_transitions, columns=columns)
        return pd.DataFrame({k: list(v) for k, v in data.items()})

    def get_run_dict(self):
        """
        get_run_dict(self) -> Dict[str, np.ndarray]

        Return unique run metadata as a dict of numpy arrays.

        :return: Dict with run_id and source_file arrays.
        :rtype: Dict[str, np.ndarray]

        Example::

            runs = xic.get_run_dict()
        """
        cdef libcpp_vector[_XICRunInfo] runs
        self.inst.get().getRuns(runs)

        cdef list run_id_list = []
        cdef list source_file_list = []
        cdef size_t i
        cdef _XICRunInfo run
        for i in range(runs.size()):
            run = runs[i]
            run_id_list.append(run.run_id)
            source_file_list.append(convOutputString(run.source_file))

        return {
            "run_id": np.array(run_id_list, dtype=np.int64),
            "source_file": np.array(source_file_list, dtype=object),
        }

    def get_run_df(self):
        """
        get_run_df(self) -> pandas.DataFrame

        Return unique run metadata as a pandas DataFrame.

        :return: DataFrame with run_id and source_file.
        :rtype: pandas.DataFrame

        :raises ImportError: If pandas is not installed
        """
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError("pandas is required for this method. Install with `pip install pandas`.") from e
        data = self.get_run_dict()
        return pd.DataFrame(data)

    def df_columns(self):
        """
        df_columns(self) -> List[str]

        Return the parquet schema column names.

        :return: List of column names.
        :rtype: List[str]
        """
        cdef StringList cols
        self.inst.get().getColumns(cols)
        return [convOutputString(cols[i]) for i in range(cols.size())]

    def to_arrow(self, explode=False):
        """
        to_arrow(self: XICParquetFile, explode: bool = False) -> pa.Table

        Returns an Apache Arrow Table representation of the chromatograms.

        If explode=True, returns long format with rt/intensity rows.

        :param explode: If True, return long format with one row per RT/intensity.
        :type explode: bool

        :return: Arrow Table with chromatogram data.
        :rtype: pyarrow.Table

        :raises ImportError: If pyarrow is not installed

        Example::

            table = xic.to_arrow()
            table_long = xic.to_arrow(explode=True)
        """
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_arrow(). "
                "Please install it with: pip install pyarrow"
            )
        return pa.Table.from_pydict(self.get_data_dict(explode=explode))

    def query_chromatograms(self):
        """
        query_chromatograms(self) -> XICParquetFile._ChromatogramQuery

        Return a chainable chromatogram query builder.

        Example::

            q = xic.query_chromatograms()
            df = q.filter_precursor_id(1381).filter_annotation("y3^1").to_df(explode=False)
        """
        from ._parquet_query import ChromatogramQuery
        return ChromatogramQuery(self)
