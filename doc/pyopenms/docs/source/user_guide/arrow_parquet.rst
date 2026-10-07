Arrow and Parquet
=================

pyOpenMS hands spectra, features and identifications to `Apache Arrow <https://arrow.apache.org/>`_ as tables,
and reads and writes the Parquet versions of the OpenMS result formats. Arrow tables convert to pandas or Polars
data frames, and are what Parquet files hold. The ``to_arrow()`` functions on this page need the ``pyarrow``
package, which the ``arrow`` extra of pyOpenMS installs:

.. code-block:: bash

    pip install "pyopenms[arrow]"

Spectra as Arrow tables
***********************

``MSExperiment.to_arrow()`` returns the peaks of all spectra as one table, with one row per peak:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    url = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms/src/data/"
    urlretrieve(url + "BSA1.mzML", "BSA1.mzML")
    exp = oms.MSExperiment()
    oms.MzMLFile().load("BSA1.mzML", exp)

    peaks = exp.to_arrow()
    print(peaks.num_rows)
    print(peaks.column_names)

.. code-block:: output

    479455
    ['mz', 'intensity', 'rt', 'spectrum_index', 'ms_level', 'native_id', 'precursor_mz', 'precursor_charge', 'precursor_intensity', 'isolation_lower', 'isolation_upper']

The table has an ``ion_mobility`` column only if the spectra carry ion mobility, which those of BSA1.mzML do
not; ``include_ion_mobility=False`` leaves it out in any case. The slower Python export, which pyOpenMS uses for
an empty experiment and, with a warning, when its compiled Arrow export is missing, always adds the column, with
NaN for spectra without ion mobility.

With ``format="semi_wide"``, the table has one row per spectrum instead, with the m/z and intensity values of its
peaks as lists. ``ms_levels``, ``min_rt``, ``max_rt``, ``min_mz`` and ``max_mz`` select what goes into the table,
for example only the :term:`MS1` peaks between m/z 400 and 800:

.. code-block:: python
    :linenos:

    ms1 = exp.to_arrow(ms_levels=[1], min_mz=400, max_mz=800)
    ms1_df = ms1.to_pandas()

``data="chromatograms"`` or ``data="both"`` export chromatograms, and a single :py:class:`~.MSSpectrum` or
:py:class:`~.MSChromatogram` has a ``to_arrow()`` of its own.

Features and identifications
****************************

A :py:class:`~.FeatureMap`, a :py:class:`~.ConsensusMap` and a :py:class:`~.PeptideIdentificationList` convert
the same way. Their ``to_arrow()`` has the columns of their ``to_df()`` (see :doc:`export_pandas_dataframe`); for
features and consensus features, the ID that ``to_df()`` uses as the index is an extra column.

.. code-block:: python
    :linenos:

    urlretrieve(url + "BSA1_F1_idmapped.featureXML", "BSA1_F1_idmapped.featureXML")
    fmap = oms.FeatureMap()
    oms.FeatureXMLFile().load("BSA1_F1_idmapped.featureXML", fmap)

    features = fmap.to_arrow()
    print(features.num_rows, features.column_names[:6])

.. code-block:: output

    256 ['peptide_sequence', 'peptide_score', 'ID_filename', 'ID_native_id', 'charge', 'rt']

For identifications, ``PeptideIdentificationList.to_psm_arrow()`` gives one row per peptide-spectrum match, with
the peptidoform in ProForma notation, its modifications, the scores and the protein accessions.

Parquet bundles
***************

OpenMS stores identifications, :term:`feature maps` and :term:`consensus maps` as Parquet *bundles*: directories named
``.idparquet``, ``.featureparquet`` or ``.consensusparquet`` that hold one Parquet file per table. The `File
formats <https://openms.readthedocs.io/en/latest/getting-started/file-formats.html>`_ page of the OpenMS
documentation lists the files, and the :term:`TOPP tools` that read and write the bundles. :py:class:`~.FileHandler`
writes and reads them like the XML formats:

.. code-block:: python
    :linenos:

    oms.FileHandler().storeFeatures("BSA1_F1.featureparquet", fmap)
    same = oms.FileHandler().loadFeatures("BSA1_F1.featureparquet")
    print(same.size())

.. code-block:: output

    256

``storeIdentifications()`` and ``loadIdentifications()`` do the same for ``.idparquet``, and
``storeConsensusFeatures()`` and ``loadConsensusFeatures()`` for ``.consensusparquet``. These functions do not need
``pyarrow``. ``.idparquet`` bundles are the native identification format (see :doc:`identification_data`):
``storeIdentifications()`` imports peptide and protein identifications into it, keeping their values and order.

The tables of an ``.idparquet`` bundle keep queries (spectra) and candidates apart. For one row per candidate with
the values of its spectrum, as the ``psms.parquet`` table of OpenMS 3.6 had them, use
:py:meth:`~.IdentificationDataFile.psm_table` (requires ``pyarrow``):

.. code-block:: python
    :linenos:

    psms = oms.IdentificationDataFile.psm_table("BSA1.idparquet").to_pandas()
    print(psms[["peptidoform", "precursor_charge", "rt", "score", "is_decoy"]].head(3))

OpenMS 3.7 does not read the four-table ``.idparquet`` bundles of OpenMS 3.6; convert them to idXML with
IDFileConverter of OpenMS 3.6.

Each table of a bundle is an ordinary Parquet file, so other programs can read it without pyOpenMS, for example
pandas:

.. code-block:: python
    :linenos:

    import pandas as pd

    table = pd.read_parquet("BSA1_F1.featureparquet/features.parquet")
    print(table.shape)

.. code-block:: output

    (256, 17)
