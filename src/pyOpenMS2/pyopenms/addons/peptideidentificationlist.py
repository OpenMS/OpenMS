"""Addon methods for PeptideIdentificationList class."""

from __future__ import annotations
import numpy as np
from . import addon


@addon("PeptideIdentificationList")
def to_df(self, decode_ontology=True, default_missing_values=None, export_unidentified=True):
    """Converts peptide identifications to a pandas DataFrame."""
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
        k_bytes = k.encode('utf-8') if isinstance(k, str) else k
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
            if export_unidentified:
                return (pep.getIdentifier(), pep.getRT(), pep.getMZ(), default_missing_values[float], default_missing_values[int],
                        default_missing_values[str], default_missing_values[str], default_missing_values[str], pep_idx, default_missing_values[int], *dmv)
            else:
                return

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

    return pd.DataFrame(np.fromiter((extract(pep, pep_idx) for pep_idx, pep in enumerate(self)), dtype=dt, count=count))


@addon("PeptideIdentificationList")
def to_arrow(self, decode_ontology=True, default_missing_values=None, export_unidentified=True):
    """Returns an Apache Arrow Table with peptide identification information."""
    import pyarrow as pa
    df = self.to_df(decode_ontology=decode_ontology,
                    default_missing_values=default_missing_values,
                    export_unidentified=export_unidentified)
    return pa.Table.from_pandas(df)


@addon("PeptideIdentificationList")
def update_scores_from_df(self, df, main_score_name):
    """Updates scores in PeptideIdentification objects from a DataFrame."""
    import pyopenms

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
    import warnings
    warnings.warn(
        "get_df() is deprecated. Use to_df() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.to_df(*args, **kwargs)
