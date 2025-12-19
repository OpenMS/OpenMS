cimport numpy as np
import numpy as np




    def __len__(self):
        """
        __len__(self: PeptideIdentificationList) -> int

        Return the number of peptide identifications in the list.

        Enables use of Python's built-in len() function.

        Returns:
            int: The number of peptide identifications in this list.
        """
        return self.size()

    def append(self, PeptideIdentification item):
        """
        append(self: PeptideIdentificationList, item: PeptideIdentification) -> None

        Add a single peptide identification to the end of the list.

        This method provides a Pythonic interface equivalent to push_back().

        Args:
            item: A single PeptideIdentification object to append.
        """
        self.push_back(item)

    def extend(self, items):
        """
        extend(self: PeptideIdentificationList, items: Iterable[PeptideIdentification]) -> None

        Add multiple peptide identifications to the end of the list.

        Args:
            items: Can be:
                - A list/iterable of PeptideIdentification objects
                - Another PeptideIdentificationList object

        Raises:
            TypeError: If items is not iterable or another PeptideIdentificationList.
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

    def get_df(self, decode_ontology=True, default_missing_values=None, export_unidentified=True):
        """
        get_df(self: PeptideIdentificationList, decode_ontology: bool = True, default_missing_values: dict = None, export_unidentified: bool = True) -> pd.DataFrame

        Converts the peptide identifications to a pandas DataFrame.

        Args:
            decode_ontology (bool): Decode meta value names using the PSI-MS ontology.
                                   Default True.
            default_missing_values (dict): Default values for missing data by type.
                                          Default: {<bool>: False, <int>: -9999, <float>: nan, <str>: ''}
            export_unidentified (bool): Export PeptideIdentifications without PeptideHit.
                                       Default True.

        Returns:
            pd.DataFrame: DataFrame with peptide identification information including:
                - id: identifier
                - rt: retention time
                - mz: mass-to-charge ratio
                - score: peptide score (column name from score type)
                - charge: peptide charge
                - protein_accession: comma-separated protein accessions
                - start: comma-separated start positions
                - end: comma-separated end positions
                - P_ID: peptide identification index
                - PSM_ID: PSM index (currently always 0, only first hit exported)
                - Additional meta value columns from PeptideHit

        Example:
            >>> peps = feature_map.get_assigned_peptide_identifications()
            >>> df = peps.get_df()

            >>> # Skip unidentified entries
            >>> df = peps.get_df(export_unidentified=False)

            >>> # Custom missing values (use Python type objects)
            >>> df = peps.get_df(default_missing_values={bool: False, int: 0, float: 0.0, str: 'NA'})
        """
        import pandas as pd
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

        Parameters:
            decode_ontology (bool): Decode meta value names using the PSI-MS ontology.
                                   Default True.
            default_missing_values (dict): Default values for missing data by type.
                                          Default: {<bool>: False, <int>: -9999, <float>: nan, <str>: ''}
            export_unidentified (bool): Export PeptideIdentifications without PeptideHit.
                                       Default True.

        Returns:
            pyarrow.Table: Arrow Table with peptide identification data.

        Example:
            >>> table = peps.to_arrow()
            >>> df = table.to_pandas()

            >>> # Convert to polars
            >>> import polars as pl
            >>> df = pl.from_arrow(table)
        """
        import pyarrow as pa
        df = self.get_df(decode_ontology=decode_ontology,
                        default_missing_values=default_missing_values,
                        export_unidentified=export_unidentified)
        return pa.Table.from_pandas(df)

    def update_scores_from_df(self, df, main_score_name):
        """
        update_scores_from_df(self: PeptideIdentificationList, df: pd.DataFrame, main_score_name: str) -> PeptideIdentificationList

        Updates the scores in PeptideIdentification objects using a pandas DataFrame.

        This method is useful when you've computed new scores (e.g., from rescoring)
        and want to update the original PeptideIdentifications.

        Args:
            df (pd.DataFrame): DataFrame obtained by converting this list to a DataFrame.
                              Minimum required columns: P_ID and a column with name
                              matching main_score_name.
            main_score_name (str): Name of the column containing the new scores.
                                  This will also be set as the score type.

        Returns:
            PeptideIdentificationList: The updated list of peptide identifications
                                      (modifies in place and returns self).

        Example:
            >>> # Get DataFrame, compute new scores, update
            >>> df = peps.get_df()
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
