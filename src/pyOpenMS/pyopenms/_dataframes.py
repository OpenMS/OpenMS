"""
DataFrame integration for pyOpenMS

Provides pandas DataFrame export functionality for OpenMS data structures.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Sequence

import numpy as np

if TYPE_CHECKING:
    import pandas as pd


class DataFrameMixin:
    """Mixin class providing DataFrame export capabilities."""

    def to_df(
        self,
        columns: Optional[Sequence[str]] = None,
        export_meta_values: bool = True,
    ) -> "pd.DataFrame":
        """
        Export data to a pandas DataFrame.

        Parameters
        ----------
        columns : sequence of str, optional
            Columns to include. If None, includes all available columns.
        export_meta_values : bool, default True
            Whether to include meta values as columns.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the data.
        """
        import pandas as pd

        data = self.get_data_dict(columns=columns, export_meta_values=export_meta_values)
        return pd.DataFrame(data)


def spectrum_to_dataframe(
    spectrum,
    columns: Optional[Sequence[str]] = None,
    export_meta_values: bool = True,
) -> "pd.DataFrame":
    """
    Convert an MSSpectrum to a pandas DataFrame.

    Parameters
    ----------
    spectrum : MSSpectrum
        The spectrum to convert.
    columns : sequence of str, optional
        Columns to include.
    export_meta_values : bool, default True
        Whether to include meta values.

    Returns
    -------
    pd.DataFrame
        DataFrame with m/z, intensity, and optional metadata.
    """
    import pandas as pd

    data = spectrum.get_data_dict(columns=columns, export_meta_values=export_meta_values)
    return pd.DataFrame(data)


def chromatogram_to_dataframe(
    chromatogram,
    columns: Optional[Sequence[str]] = None,
    export_meta_values: bool = True,
) -> "pd.DataFrame":
    """
    Convert an MSChromatogram to a pandas DataFrame.

    Parameters
    ----------
    chromatogram : MSChromatogram
        The chromatogram to convert.
    columns : sequence of str, optional
        Columns to include.
    export_meta_values : bool, default True
        Whether to include meta values.

    Returns
    -------
    pd.DataFrame
        DataFrame with RT, intensity, and optional metadata.
    """
    import pandas as pd

    data = chromatogram.get_data_dict(
        columns=columns, export_meta_values=export_meta_values
    )
    return pd.DataFrame(data)


def experiment_to_dataframe(
    experiment,
    long_format: bool = True,
    include_precursor: bool = True,
) -> "pd.DataFrame":
    """
    Convert an MSExperiment to a pandas DataFrame.

    Parameters
    ----------
    experiment : MSExperiment
        The experiment to convert.
    long_format : bool, default True
        If True, returns long format with one row per peak.
        If False, returns wide format with one row per spectrum.
    include_precursor : bool, default True
        Whether to include precursor information for MS2 spectra.

    Returns
    -------
    pd.DataFrame
        DataFrame containing all spectra data.
    """
    import pandas as pd

    if long_format:
        records = []
        for i, spec in enumerate(experiment):
            mz, intensity = spec.get_peaks()
            n_peaks = len(mz)
            if n_peaks == 0:
                continue

            record = {
                "spectrum_index": np.full(n_peaks, i, dtype=np.int32),
                "mz": mz,
                "intensity": intensity,
                "rt": np.full(n_peaks, spec.getRT(), dtype=np.float64),
                "ms_level": np.full(n_peaks, spec.getMSLevel(), dtype=np.int32),
            }

            if include_precursor and spec.getMSLevel() > 1:
                precursors = spec.getPrecursors()
                if precursors:
                    precursor = precursors[0]
                    record["precursor_mz"] = np.full(
                        n_peaks, precursor.getMZ(), dtype=np.float64
                    )
                    record["precursor_charge"] = np.full(
                        n_peaks, precursor.getCharge(), dtype=np.int32
                    )

            records.append(record)

        if not records:
            return pd.DataFrame()

        # Concatenate all records
        data = {
            key: np.concatenate([r[key] for r in records if key in r])
            for key in records[0].keys()
        }
        return pd.DataFrame(data)
    else:
        records = []
        for i, spec in enumerate(experiment):
            mz, intensity = spec.get_peaks()
            record = {
                "spectrum_index": i,
                "rt": spec.getRT(),
                "ms_level": spec.getMSLevel(),
                "n_peaks": len(mz),
                "tic": intensity.sum() if len(intensity) > 0 else 0.0,
                "base_peak_mz": mz[intensity.argmax()] if len(intensity) > 0 else 0.0,
                "base_peak_intensity": intensity.max() if len(intensity) > 0 else 0.0,
            }

            if include_precursor and spec.getMSLevel() > 1:
                precursors = spec.getPrecursors()
                if precursors:
                    precursor = precursors[0]
                    record["precursor_mz"] = precursor.getMZ()
                    record["precursor_charge"] = precursor.getCharge()

            records.append(record)

        return pd.DataFrame(records)


def feature_map_to_dataframe(
    feature_map,
    include_subordinates: bool = False,
) -> "pd.DataFrame":
    """
    Convert a FeatureMap to a pandas DataFrame.

    Parameters
    ----------
    feature_map : FeatureMap
        The feature map to convert.
    include_subordinates : bool, default False
        Whether to include subordinate features.

    Returns
    -------
    pd.DataFrame
        DataFrame with feature data.
    """
    import pandas as pd

    records = []
    for feature in feature_map:
        record = {
            "rt": feature.getRT(),
            "mz": feature.getMZ(),
            "intensity": feature.getIntensity(),
            "charge": feature.getCharge(),
            "quality": feature.getOverallQuality(),
            "width": feature.getWidth(),
        }

        # Add convex hull bounds if available
        hulls = feature.getConvexHulls()
        if hulls:
            hull = hulls[0]
            bbox = hull.getBoundingBox()
            record["rt_start"] = bbox.minPosition()[0]
            record["rt_end"] = bbox.maxPosition()[0]
            record["mz_start"] = bbox.minPosition()[1]
            record["mz_end"] = bbox.maxPosition()[1]

        records.append(record)

        if include_subordinates:
            for sub in feature.getSubordinates():
                sub_record = {
                    "rt": sub.getRT(),
                    "mz": sub.getMZ(),
                    "intensity": sub.getIntensity(),
                    "charge": sub.getCharge(),
                    "quality": sub.getOverallQuality(),
                    "width": sub.getWidth(),
                    "parent_rt": feature.getRT(),
                    "parent_mz": feature.getMZ(),
                }
                records.append(sub_record)

    return pd.DataFrame(records)


def consensus_map_to_dataframe(consensus_map) -> "pd.DataFrame":
    """
    Convert a ConsensusMap to a pandas DataFrame.

    Parameters
    ----------
    consensus_map : ConsensusMap
        The consensus map to convert.

    Returns
    -------
    pd.DataFrame
        DataFrame with consensus feature data.
    """
    import pandas as pd

    records = []
    for cf in consensus_map:
        record = {
            "rt": cf.getRT(),
            "mz": cf.getMZ(),
            "intensity": cf.getIntensity(),
            "charge": cf.getCharge(),
            "quality": cf.getQuality(),
            "n_features": cf.size(),
        }
        records.append(record)

    return pd.DataFrame(records)
