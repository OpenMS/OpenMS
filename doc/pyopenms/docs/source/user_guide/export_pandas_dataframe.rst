Export to pandas DataFrame
==========================

**NOTE: This feature is available only if using a version of pyOpenMS >= 3.0, at the time of writing this means using
the nightly builds as described in the**
`Installation Instructions <installation.html#nightly-ci-wheels>`_.

In pyOpenMS some data structures can be converted to a tabular format as a ``pandas.DataFrame``.
This allows convenient access to data and meta values of spectra, features and identifications.

Required imports for the examples:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    import pandas as pd
    from urllib.request import urlretrieve

    url = "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master/src/data/"


MSExperiment
************

**pyopenms.MSExperiment.get_df(** *long=False* **)**
        Generates a pandas DataFrame with all peaks  in the MSExperiment

        **Parameters:**

        **long :** default False
        
        set to True if you want to have a long/expanded/melted dataframe with one row per peak. Faster but
        replicated RT information. If False, returns rows in the style: rt, np.array(mz), np.array(int)
        
        **Returns:**

        **pandas.DataFrame** 
        
        peak map information stored in a DataFrame

**Examples:**

.. code-block:: python
    :linenos:

    urlretrieve(url + "BSA1.mzML", "BSA1.mzML")
    exp = oms.MSExperiment()
    oms.MzMLFile().load("BSA1.mzML", exp)

    df = exp.get_df()  # default: long = False
    df.head(2)


.. csv-table:: exp.get_df()
   :widths: 2 10 50 50
   :header: ,"RT", "mzarray", "intarray"

   "0",	"1501.41394", "[300.0897645621494, 300.18132740129533, 300.20...",	"[3431.0261, 1181.809, 1516.1746, 1719.8547, 11..."
   "1", "1503.03125", "[300.06577092599525, 300.08932376441896, 300.2...",	"[914.79034, 1842.2311, 2395.1025, 851.4738, 16..." 


.. code-block:: python
    :linenos:

    df = exp.get_df(long=True)
    df.head(2)

.. csv-table:: exp.get_df(long=True)
   :widths: 2 20 20 20
   :header: "",	"RT",	"mz", "inty"

   "0",	"1501.41394",	"300.089752",	"3431.026123"
   "1",	"1501.41394",	"300.181335",	"1181.808960"

PeptideIdentification
*********************

**pyopenms.peptide_identifications_to_df( peps**, *decode_ontology=True*, *default_missing_values={bool: False, int: -9999, float: np.nan, str: ''}*, *export_unidentified=True* **)**
        Generates a pandas DataFrame with all peaks  in the MSExperiment

        **Parameters:**

        **peps :** 
        
        list of PeptideIdentification objects

        **decode_ontology :** default True
        
        if meta values contain CV identifer (e.g., from PSI-MS) they will be automatically decoded into the human readable CV term name.

        **default_missing_values :** default {bool: False, int: -9999, float: np.nan, str: ''}
        
        default value for missing values for each data type

        **export_unidentified :** default True
        
        export PeptideIdentifications without PeptideHit
        
        **Returns:**

        **pandas.DataFrame** 
        
        peptide identifications in a DataFrame

**Example:**

.. code-block:: python
    :linenos:

    urlretrieve(url + "small.idXML", "small.idXML")
    prot_ids = []
    pep_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load("small.idXML", prot_ids, pep_ids)

    df = oms.peptide_identifications_to_df(pep_ids)
    df.head(2)
.. csv-table:: peptide_identifications_to_df(pep_ids)
   :widths: 2 20 10 20 20 10 20 20 20 20 20 20 20 20 20 20
   :header: "",	"id",	"RT",	"mz",	"q-value",	"charge",	"protein_accession",	"start",	"end",	"NuXL:z2 mass",	"NuXL:z3 mass",	"...", "isotope_error",	"NuXL:peptide_mass_z0",	"NuXL:XL_U",	"NuXL:sequence_score"

    "0",	"OpenNuXL_2019-12-04T16:39:43_1021782429466859437",	"900.425415",	"414.730865",	"0.368649",	"4",	"DECOY_sp\|Q86UQ0|ZN589_HUMAN",	"255",	"267",	"828.458069",	"552.641113",	"...",	"0",	"1654.901611",	"0",	"0.173912"
    "1",	"OpenNuXL_2019-12-04T16:39:43_7293634134684008928",	"903.565186",	"506.259521",	"0.422779",	"2",	"sp\|P61313|RL15_HUMAN",	"179",	"187",	"0.0",	"0.0",	"...",	"0",	"1010.504639",	"0",	"0.290786"

FeatureMap
**********

**pyopenms.FeatureMap.get_df(** *meta_values = None* **)**
        Generates a pandas DataFrame with information contained in the FeatureMap.

        Optionally the feature meta values and information for the assigned PeptideHit can be exported.

        **Parameters:**

        **meta_values :** default None
        
        meta values to include (None, [custom list of meta value names] or 'all')

        **export_peptide_identifications (bool):** default True
        
        Export sequence and score for best PeptideHit assigned to a feature.
        Additionally the ID_filename (file name of the corresponding ProteinIdentification) and the ID_native_id 
        (spectrum ID of the corresponding Feature) are exported. They are also annotated as meta values when 
        collecting all assigned PeptideIdentifications from a FeatureMap with FeatureMap.get_assigned_peptide_identifications().
        A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides) can be
        merged with the FeatureMap DataFrame with:
        merged_df = pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])
        
        **Returns:**

        **pandas.DataFrame** 
        
        feature information stored in a DataFrame

**Examples:**
   
.. code-block:: python
    :linenos:

    urlretrieve(url + "BSA1_F1_idmapped.featureXML", "BSA1_F1_idmapped.featureXML")
    feature_map = oms.FeatureMap()
    oms.FeatureXMLFile().load("BSA1_F1_idmapped.featureXML", feature_map)

    df = feature_map.get_df()  # default: meta_values = None
    df.head(2)
.. csv-table:: feature_map.get_df()
   :widths: 20 20 20 20 20 5 20 20 20 20 20 20 20 20
   :header: "id",	"peptide_sequence",	"peptide_score",	"ID_filename",	"ID_native_id",	"charge",	"RT",	"mz",	"RTstart",	"RTend",	"mzstart",	"mzend",	"quality",	"intensity"

   "9650885788371886430",	"LVTDLTK",	"0.000000",	"unknown",	"spectrum=1270",	"2",	"1942.600083",	"395.239277",	"1932.484009",	"1950.834351",	"395.239199",	"397.245758",	"0.808494",	"157572000.0"
   "18416216708636999474",	"DDSPDLPK",	"0.034483",	"unknown",	"spectrum=1167",	"2",	"1749.138335",	"443.711224",	"1735.693115",	"1763.343506",	"443.711122",	"445.717531",	"0.893553",	"54069300.0"


.. code-block:: python
    :linenos:

    df = feature_map.get_df(meta_values="all", export_peptide_identifications=False)
    df.head(2)

.. csv-table:: feature_map.get_df(meta_values = 'all', export_peptide_identifications = False)
   :widths: 20 5 20 20 20 20 20 20 20 20 20 20 20 20 20 20
   :header: "id",	"charge",	"RT",	"mz",	"RTstart",	"RTend",	"mzstart",	"mzend",	"quality",	"intensity",	"FWHM",	"spectrum_index",	"spectrum_native_id",	"label",	"score_correlation",	"score_fit"

   "9650885788371886430",	"2",	"1942.600083",	"395.239277",	"1932.484009",	"1950.834351",	"395.239199",	"397.245758",	"0.808494",	"157572000.0",	"10.061090",	"259",	"spectrum=1270",	"168",	"0.989969",	"0.660286"
   "18416216708636999474",	"2",	"1749.138335",	"443.711224",	"1735.693115",	"1763.343506",	"443.71112",	"445.717531",	"0.893553",	"54069300.0", "14.156094",	"156",	"spectrum=1167",	"169",	"0.999002",	"0.799234"

.. code-block:: python
    :linenos:

    df = feature_map.get_df(meta_values=[b"FWHM", b"label"])
    df.head(2)

.. csv-table:: feature_map.get_df(meta_values = [b'FWHM', b'label'])
   :widths: 20 5 20 20 20 20 20 20 20 20 20 20
   :header: "id",	"charge",	"RT",	"mz",	"RTstart",	"RTend",	"mzstart",	"mzend",	"quality",	"intensity", "FWHM",	"label"

   "9650885788371886430",	"2",	"1942.600083",	"395.239277",	"1932.484009",	"1950.834351",	"395.239199",	"397.245758",	"0.808494",	"157572000.0",	"10.061090",	"168"
   "18416216708636999474",	"2",	"1749.138335",	"443.711224",	"1735.693115",	"1763.343506",	"443.71112",	"445.717531",	"0.893553",	"54069300.0",	"14.156094",	"169"

**Extract assigned peptide identifications from a feature map**

Peptide identifications can be mapped to their corresponding features in a ``FeatureMap``. It is possible to extract them using the function
``pyopenms.FeatureMap.get_assigned_peptide_identifications()`` returning a list of ``PeptideIdentification`` objects.


**pyopenms.FeatureMap.get_assigned_peptide_identifications()**
        Generates a list with peptide identifications assigned to a feature.

        Adds 'ID_native_id' (feature spectrum id), 'ID_filename' (primary MS run path of corresponding ProteinIdentification)
        and 'feature_id' (unique ID of corresponding Feature) as meta values to the peptide hits.
        A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides) can be
        merged with the FeatureMap DataFrame with:
        merged_df = pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])
        
        **Returns:**

        **[PeptideIdentification]** 
        
        list of PeptideIdentification objects

A ``DataFrame`` can be created on the resulting list of `PeptideIdentification` objects using
``pyopenms.peptide_identifications_to_df(assigned_peptides)``.
:term:`Feature map<feature map>` and peptide data frames contain columns, on which they can be merged together to contain the complete
information for peptides and features in a single data frame.

The columns for unambiguously merging the data frames:

- ``feature_id``: the unique feature identifier

- ``ID_native_id``: the feature spectrum native identifier

- ``ID_filename``: the filename (primary MS run path) of the corresponding `ProteinIdentification`

**Example:**

.. code-block:: python
    :linenos:

    feature_df = feature_map.get_df()
    assigned_peptides = feature_map.get_assigned_peptide_identifications()
    assigned_peptide_df = oms.peptide_identifications_to_df(assigned_peptides)

    merged_df = pd.merge(
        feature_df,
        assigned_peptide_df,
        on=["feature_id", "ID_native_id", "ID_filename"],
    )
    merged_df.head(2)

.. csv-table:: consensus_map.get_df()
   :widths: 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20
   :header: "feature_id",	"peptide_sequence",	"peptide_score",	"ID_filename",	"ID_native_id",	"charge_x",	"RT_x",	"mz_x",	"RTstart",	"RTend",	"...",	"id",	"RT_y",	"mz_y",	"q-value",	"charge_y",	"protein_accession",	"start",	"end",	"OMSSA_score",	"target_decoy"

   "9650885788371886430",	"LVTDLTK",	"0.000000",	"unknown",	"spectrum=1270",	"2",	"1942.600083",	"395.239277",	"1932.484009",	"1950.834351",	"...",	"OMSSA_2009-11-17T11:11:11_4731105163044641872",	"1933.405151",	"395.239349",	"0.000000",	"2",	"P02769|ALBU_BOVIN",	"-1",	"-1",	"0.001084",	"True"
   "18416216708636999474",	"DDSPDLPK",	"0.034483",	"unknown",	"spectrum=1167",	"2",	"1749.138335",	"443.711224",	"1735.693115",	"1763.343506",	"...",	"OMSSA_2009-11-17T11:11:11_4731105163044641872",	"1738.033447",	"443.711243",	"0.034483",	"2",	"P02769|ALBU_BOVIN",	"-1",	"-1",	"0.003951",	"True"    

ConsensusMap
************

**pyopenms.ConsensusMap.get_df()**
        Generates a pandas DataFrame with both consensus feature meta data and intensities from each sample.

        **Returns:**

        **pandas.DataFrame** 
        
        :term:`consensus map` meta data and intensity stored in pandas DataFrame

**pyopenms.ConsensusMap.get_intensity_df()**
        Generates a pandas DataFrame with feature intensities from each sample in long format (over files).

        For labelled analyses channel intensities will be in one row, therefore resulting in a semi-long/block format.
        Resulting DataFrame can be joined with result from get_metadata_df by their index 'id'.

        **Returns:**

        **pandas.DataFrame** 
        
        intensity DataFrame

**pyopenms.ConsensusMap.get_metadata_df()**
        Generates a pandas DataFrame with feature meta data (sequence, charge, mz, RT, quality).

        Resulting DataFrame can be joined with result from get_intensity_df by their index 'id'.

        **Returns:**

        **pandas.DataFrame** 
        
        DataFrame with metadata for each feature (such as: best identified sequence, charge, centroid RT/mz, fitting quality)

**Examples:**

.. code-block:: python
    :linenos:

    urlretrieve(
        url + "ProteomicsLFQ_1_out.consensusXML", "ProteomicsLFQ_1_out.consensusXML"
    )
    consensus_map = oms.ConsensusMap()
    oms.ConsensusXMLFile().load("ProteomicsLFQ_1_out.consensusXML", consensus_map)

    df = consensus_map.get_df()
    df.head(2)
.. csv-table:: consensus_map.get_df()
   :widths: 2 10 20 20 20 20 30 10 30
   :header: "id",	"sequence",	"charge",	"RT",	"mz",	"quality",	"BSA1_F1.mzML",	"...",	"BSA1_F2.mzML"

   "2935923263525422257",	"DGDIEAEISR",	"3",	"1523.370634",	"368.843773",	"0.000000",	"0.0",	"...",	"0.0"
   "10409195546240342212",	"SHC(Carbamidomethyl)IAEVEK",	"3",	"1552.032973",	"358.174576",	"0.491247",	"1358151.0",	"...",	"0.0"

.. code-block:: python
    :linenos:

    df = consensus_map.get_intensity_df()
    df.head(2)

.. csv-table:: consensus_map.get_intensity_df()
   :widths: 20 30 10 30
   :header: "id",	"BSA1_F1.mzML",	"...",	"BSA1_F2.mzML"

   "2935923263525422257",	"0.0",	"...",	"0.0"
   "10409195546240342212",	"1358151.0",	"...",	"0.0"

.. code-block:: python
    :linenos:

    df = consensus_map.get_metadata_df()
    df.head(2)

.. csv-table:: consensus_map.get_metadata_df()
   :widths: 20 20 20 20 20 20
   :header: "id",	"sequence",	"charge",	"RT",	"mz",	"quality"

   "2935923263525422257",	"DGDIEAEISR",	"3",	"1523.370634",	"368.843773",	"0.000000"
   "10409195546240342212",	"SHC(Carbamidomethyl)IAEVEK",	"3",	"1552.032973",	"358.174576",	"0.491247"
