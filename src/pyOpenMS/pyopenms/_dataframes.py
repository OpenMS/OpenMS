"""DataFrame export utilities for pyOpenMS.

This module provides utility functions for converting OpenMS data structures to pandas DataFrames.
The get_df() methods are now implemented directly in the Cython classes (MSSpectrum, MSChromatogram,
Mobilogram, MSExperiment, ConsensusMap, FeatureMap, MRMTransitionGroupCP).

This module provides:
- peptide_identifications_to_df: Convert PeptideIdentification list to DataFrame
- update_scores_from_df: Update scores in PeptideIdentifications from DataFrame
- common_meta_value_types: Dictionary for numpy type mapping of common meta values
"""

__all__ = [
    'peptide_identifications_to_df',
    'update_scores_from_df',
    'common_meta_value_types',
    '_add_meta_values',
]

from typing import List, Union

from . import ConsensusMap
from . import FeatureMap
from . import MSExperiment
from . import PeakMap
from . import PeptideIdentificationList as _PeptideIdentificationList
from . import PeptideIdentification as _PeptideIdentification
from . import ControlledVocabulary as _ControlledVocabulary
from . import File as _File
from . import MSSpectrum
from . import PeakSpectrum
from . import MSChromatogram
from . import MRMTransitionGroupCP
from . import Mobilogram
from . import DataValue as _DataValue

import pandas as _pd
import numpy as _np


# Common meta value types for numpy type mapping
common_meta_value_types = {
    b'label': 'U50',
    b'spectrum_index': 'i',
    b'score_fit': 'f',
    b'score_correlation': 'f',
    b'FWHM': 'f',
    b'spectrum_native_id': 'U100',
    b'max_height': 'f',
    b'num_of_masstraces': 'i',
    b'masstrace_intensity': 'f',
    b'Group': 'U50',
    b'is_ungrouped_monoisotopic': 'i',
    b'leftWidth': 'f',
    b'rightWidth': 'f',
    b'total_xic': 'f',
    b'PeptideRef': 'U100',
    b'peak_apices_sum': 'f'
}
"""Global dict to define which autoconversion to numpy types is tried for certain metavalues.

This can be changed to your liking but only affects future exports of any OpenMS datastructure to dataframes.
Especially string lengths (i.e., U types) benefit from adaption to save memory. The default type is currently
hardcoded to U50 (i.e., 50 unicode characters)
"""


def peptide_identifications_to_df(peps: _PeptideIdentificationList, decode_ontology: bool = True,
                                  default_missing_values: dict = None,
                                  export_unidentified: bool = True):
    """Converts a list of peptide identifications to a pandas DataFrame.

    Parameters:
    peps (PeptideIdentificationList): list of PeptideIdentification objects
    decode_ontology (bool): decode meta value names
    default_missing_values: default value for missing values for each data type
    export_unidentified: export PeptideIdentifications without PeptideHit

    Returns:
    pandas.DataFrame: peptide identifications in a DataFrame
    """

    if default_missing_values is None:
        default_missing_values = {bool: False, int: -9999, float: _np.nan, str: ''}

    switchDict = {bool: '?', int: 'i', float: 'f', str: 'U100'}

    # filter out PeptideIdentifications without PeptideHits if export_unidentified == False
    count = peps.size()
    if not export_unidentified:
        count = sum(len(pep.getHits()) > 0 for pep in peps)

    # get all possible metavalues
    metavals = []
    types = []
    mainscorename = "score"
    for pep in peps:
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
            for p in peps:
                hits = p.getHits()
                if not len(hits) == 0:
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

        ret += [str(pep_idx), 0] # we currently only export the first hit

        for k in metavals:
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
                ret.append(default_missing_values[type(val)])
        return tuple(ret)

    return _pd.DataFrame(_np.fromiter((extract(pep, pep_idx) for pep_idx, pep in enumerate(peps)), dtype=dt, count=count))


def update_scores_from_df(peps: _PeptideIdentificationList, df: _pd.DataFrame, main_score_name: str):
    """
    Updates the scores in PeptideIdentification objects using a pandas dataframe.

    :param peps: list of PeptideIdentification objects
    :param df: pandas dataframe obtained by converting peps to a dataframe. Minimum required: P_ID column and column with name passed by main_score_name
    :return: the updated list of peptide identifications
    """

    rets = peps

    for index, row in df.iterrows():
        pid_index = int(row["P_ID"])
        pi = _PeptideIdentification(peps[pid_index])
        pi.setScoreType(main_score_name)
        hits = pi.getHits()
        if len(hits) > 0:
            best_hit = hits[0]
            best_hit.setScore(float(row[main_score_name]))
            hits[0] = best_hit
            pi.setHits(hits)

        rets[pid_index] = pi

    return rets


def _add_meta_values(df: _pd.DataFrame, object: any) -> _pd.DataFrame:
    """
    Adds metavalues from given object to given DataFrame.

    Args:
        df (pd.DataFrame): DataFrame to which metavalues will be added.
        object (any): Object from which metavalues will be extracted.

    Returns:
        pd.DataFrame: DataFrame with added meta values.
    """
    mvs = []
    object.getKeys(mvs)

    for k in mvs:
        dv = object.getMetaValue(k)
        col_name = k.decode()

        try:
            # Handle native Python types (returned by autowrap)
            if isinstance(dv, float):
                value = dv
                dtype = "float64"
            elif isinstance(dv, int):
                value = dv
                dtype = "int64"
            elif isinstance(dv, str):
                value = dv
                dtype = f"U{max(1, len(value))}"
            elif isinstance(dv, bytes):
                value = dv.decode()
                dtype = f"U{max(1, len(value))}"
            # Handle DataValue objects (if ever returned)
            elif hasattr(dv, 'valueType'):
                if dv.valueType() == _DataValue.STRING_VALUE:
                    value = dv.toString().decode()
                    dtype = f"U{max(1, len(value))}"
                elif dv.valueType() == _DataValue.INT_VALUE:
                    value = dv.toInt()
                    dtype = "int32"
                elif dv.valueType() == _DataValue.DOUBLE_VALUE:
                    value = dv.toDouble()
                    dtype = "float64"
                elif dv.valueType() == _DataValue.EMPTY_VALUE:
                    continue
                else:
                    value = str(dv)
                    dtype = "object"
            else:
                value = str(dv)
                dtype = "object"

            df[col_name] = _np.full(df.shape[0], value, dtype=dtype)

        except Exception:
            df[col_name] = _np.full(df.shape[0], str(dv), dtype='object')

    return df
