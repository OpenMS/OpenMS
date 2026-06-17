Identification by Accurate Mass
===============================
Example workflow for the processing of a set of :term:`mzML` files (defined in the ``files`` variable) including centroiding,
features detection, :term:feature: linking and accurate mass search.
The resulting data gets processed in a pandas data frame with features filtering (missing values, quality) and imputation
of remaining missing values.
Compounds detected during accurate mass search will be annotated in the resulting dataframe.

Imports and :term:`mzML` file path
**********************************

.. code-block:: python
    :linenos:

    import os
    import shutil
    import requests
    import pandas as pd
    import pyopenms as oms
    import numpy as np
    from sklearn.impute import KNNImputer
    from sklearn.preprocessing import FunctionTransformer
    from sklearn.pipeline import Pipeline
    import plotly.graph_objects as go
    import plotly.express as px
    import matplotlib.pyplot as plt

    # set path to your mzML files, or leave like this to use the example data
    files = os.path.join(os.getcwd(), "IdByMz_Example")


Download Example Data
*********************
Execute this cell only for the example workflow.

.. code-block:: python
    :linenos:

    if not os.path.isdir(os.path.join(os.getcwd(), "IdByMz_Example")):
        os.mkdir(os.path.join(os.getcwd(), "IdByMz_Example"))

    base = "https://abibuilder.cs.uni-tuebingen.de/archive/openms/Tutorials/Example_Data/Metabolomics/"
    urls = [
        "datasets/2012_02_03_PStd_050_1.mzML",
        "datasets/2012_02_03_PStd_050_2.mzML",
        "datasets/2012_02_03_PStd_050_3.mzML",
        "databases/PositiveAdducts.tsv",
        "databases/NegativeAdducts.tsv",
        "databases/HMDBMappingFile.tsv",
        "databases/HMDB2StructMapping.tsv",
    ]

    for url in urls:
        request = requests.get(base + url, allow_redirects=True)
        open(os.path.join(files, os.path.basename(url)), "wb").write(
            request.content
        )

Centroiding
***********
If files are already centroided this step can bet omitted.

``in``: path to MS data (files)

``out``: path to centroided :term:`mzML` files in a subfolder 'centroid' (files)

.. code-block:: python
    :linenos:

    if os.path.exists(os.path.join(files, "centroid")):
        shutil.rmtree(os.path.join(files, "centroid"))
    os.mkdir(os.path.join(files, "centroid"))

    for file in os.listdir(files):
        if file.endswith(".mzML"):
            exp_raw = oms.MSExperiment()
            oms.MzMLFile().load(os.path.join(files, file), exp_raw)
            exp_centroid = oms.MSExperiment()

            oms.PeakPickerHiRes().pickExperiment(exp_raw, exp_centroid, True)

            oms.MzMLFile().store(os.path.join(files, "centroid", file), exp_centroid)
            del exp_raw

    files = os.path.join(files, "centroid")

Feature Detection
*****************
``in``: path to centroid :term:`mzML` files (files)

``out``: list of :py:class:`~.FeatureMap` (feature_maps)

.. code-block:: python
    :linenos:

    feature_maps = []

    for file in os.listdir(files):
        if file.endswith(".mzML"):
            exp = oms.MSExperiment()
            oms.MzMLFile().load(os.path.join(files, file), exp)

            exp.sortSpectra(True)

            mass_traces = []
            mtd = oms.MassTraceDetection()
            mtd_params = mtd.getDefaults()
            mtd_params.setValue(
                "mass_error_ppm", 5.0
            )  # set according to your instrument mass error
            mtd_params.setValue(
                "noise_threshold_int", 1000.0
            )  # adjust to noise level in your data
            mtd.setParameters(mtd_params)
            mtd.run(exp, mass_traces, 0)

            mass_traces_split = []
            mass_traces_final = []
            epd = oms.ElutionPeakDetection()
            epd_params = epd.getDefaults()
            epd_params.setValue("width_filtering", "fixed")
            epd.setParameters(epd_params)
            epd.detectPeaks(mass_traces, mass_traces_split)

            if epd.getParameters().getValue("width_filtering") == "auto":
                epd.filterByPeakWidth(mass_traces_split, mass_traces_final)
            else:
                mass_traces_final = mass_traces_split

            feature_map = oms.FeatureMap()
            feat_chrom = []
            ffm = oms.FeatureFindingMetabo()
            ffm_params = ffm.getDefaults()
            ffm_params.setValue("isotope_filtering_model", "none")
            ffm_params.setValue(
                "remove_single_traces", "true"
            )  # set false to keep features with only one mass trace
            ffm_params.setValue("mz_scoring_by_elements", "false")
            ffm_params.setValue("report_convex_hulls", "true")
            ffm.setParameters(ffm_params)
            ffm.run(mass_traces_final, feature_map, feat_chrom)

            feature_map.setUniqueIds()
            feature_map.setPrimaryMSRunPath([file[:-5].encode()])

            feature_maps.append(feature_map)

Feature Map Retention Time Alignment
************************************
``in``: unaligned list of :py:class:`~.FeatureMap` (feature_maps)

``out``: list of :py:class:`~.FeatureMap` aligned to the first :term:`feature map` in the list (feature_maps)

.. code-block:: python
    :linenos:

    # get in index of feature map with highest number of features in feature map list
    ref_index = [
        i[0]
        for i in sorted(
            enumerate([fm.size() for fm in feature_maps]), key=lambda x: x[1]
        )
    ][-1]

    aligner = oms.MapAlignmentAlgorithmPoseClustering()

    aligner.setReference(feature_maps[ref_index])

    for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1 :]:
        trafo = oms.TransformationDescription()
        aligner.align(feature_map, trafo)
        transformer = oms.MapAlignmentTransformer()
        transformer.transformRetentionTimes(
            feature_map, trafo, True
        )  # store original RT as meta value

Visualization of RTs before and after Alignment
***********************************************

.. code-block:: python
    :linenos:

    fmaps = (
        [feature_maps[ref_index]]
        + feature_maps[:ref_index]
        + feature_maps[ref_index + 1 :]
    )

    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(1, 2, 1)
    ax.set_title("consensus map before alignment")
    ax.set_ylabel("m/z")
    ax.set_xlabel("RT")

    # use alpha value to display feature intensity
    ax.scatter(
        [f.getRT() for f in fmaps[0]],
        [f.getMZ() for f in fmaps[0]],
        alpha=np.asarray([f.getIntensity() for f in fmaps[0]])
        / max([f.getIntensity() for f in fmaps[0]]),
    )

    for fm in fmaps[1:]:
        ax.scatter(
            [f.getMetaValue("original_RT") for f in fm],
            [f.getMZ() for f in fm],
            alpha=np.asarray([f.getIntensity() for f in fm])
            / max([f.getIntensity() for f in fm]),
        )

    ax = fig.add_subplot(1, 2, 2)
    ax.set_title("consensus map after alignment")
    ax.set_xlabel("RT")

    for fm in fmaps:
        ax.scatter(
            [f.getRT() for f in fm],
            [f.getMZ() for f in fm],
            alpha=np.asarray([f.getIntensity() for f in fm])
            / max([f.getIntensity() for f in fm]),
        )

    fig.tight_layout()
    fig.legend(
        [fmap.getMetaValue("spectra_data")[0].decode() for fmap in fmaps],
        loc="lower center",
    )
    # in some cases get file name elsewhere, e.g. fmap.getDataProcessing()[0].getMetaValue('parameter: out')
    fig.show()

Feature Linking
***************
``in``: list of:py:class:`~.FeatureMap` (feature_maps)

``out``: :py:class:`~.ConsensusMap` (consensus_map)

.. code-block:: python
    :linenos:

    feature_grouper = oms.FeatureGroupingAlgorithmQT()

    consensus_map = oms.ConsensusMap()
    file_descriptions = consensus_map.getColumnHeaders()

    for i, feature_map in enumerate(feature_maps):
        file_description = file_descriptions.get(i, oms.ColumnHeader())
        file_description.filename = feature_map.getMetaValue("spectra_data")[
            0
        ].decode()
        file_description.size = feature_map.size()
        file_description.unique_id = feature_map.getUniqueId()
        file_descriptions[i] = file_description

    consensus_map.setColumnHeaders(file_descriptions)
    feature_grouper.group(feature_maps, consensus_map)

ConsensusMap to Pandas DataFrame
********************************
``in``: :py:class:`~.ConsensusMap` (consensus_map)

``out``: DataFrame with RT, mz and quality from :py:class:`~.ConsensusMap` (cm_df)

.. code-block:: python
    :linenos:

    intensities = consensus_map.get_intensity_df()

    meta_data = consensus_map.get_metadata_df()[["RT", "mz", "quality"]]

    cm_df = pd.concat([meta_data, intensities], axis=1)
    cm_df.reset_index(drop=True, inplace=True)
    cm_df

Accurate Mass Search
********************
``in``: :py:class:`~.ConsensusMap` (consensus_map)

``out``: DataFrame with :py:class:`~.AccurateMassSearchEngine` results (ams_df)

.. code-block:: python
    :linenos:

    if files.endswith("centroid"):
        files = os.path.join(files, "..")

    ams = oms.AccurateMassSearchEngine()

    ams_params = ams.getParameters()
    ams_params.setValue("ionization_mode", "negative")
    ams_params.setValue(
        "positive_adducts", os.path.join(files, "PositiveAdducts.tsv")
    )
    ams_params.setValue(
        "negative_adducts", os.path.join(files, "NegativeAdducts.tsv")
    )
    ams_params.setValue("db:mapping", [os.path.join(files, "HMDBMappingFile.tsv")])
    ams_params.setValue(
        "db:struct", [os.path.join(files, "HMDB2StructMapping.tsv")]
    )
    ams.setParameters(ams_params)

    mztab = oms.MzTab()

    ams.init()

    ams.run(consensus_map, mztab)

    oms.MzTabFile().store(os.path.join(files, "ids.tsv"), mztab)

    with open(os.path.join(files, "ids_smsection.tsv"), "w") as output, open(
        os.path.join(files, "ids.tsv"), "r"
    ) as input:
        for line in input:
            if line.lstrip().startswith("SM"):
                output.write(line[4:])

    ams_df = pd.read_csv(os.path.join(files, "ids_smsection.tsv"), sep="\t")

    os.remove(os.path.join(files, "ids.tsv"))
    os.remove(os.path.join(files, "ids_smsection.tsv"))

    ams_df

Data Filtering and Imputation
*****************************
``in``: unfiltered :py:class:`~.ConsensusMap` DataFrame (cm_df)

``out``: features below minimum quality and with too many missing values removed,
remaining missing values imputed with KNN algorithm (cm_df)

.. code-block:: python
    :linenos:

    allowed_missing_values = 1
    min_feature_quality = 0.8
    n_nearest_neighbours = 2

    # drop features that have more then the allowed number of missing values or are below minimum feature quality
    to_drop = []

    for i, row in cm_df.iterrows():
        if (
            row.isna().sum() > allowed_missing_values
            or row["quality"] < min_feature_quality
        ):
            to_drop.append(i)

    cm_df.drop(index=cm_df.index[to_drop], inplace=True)

    # Data imputation with KNN
    imputer = Pipeline(
        [
            ("imputer", KNNImputer(n_neighbors=2)),
            (
                "pandarizer",
                FunctionTransformer(
                    lambda x: pd.DataFrame(x, columns=cm_df.columns)
                ),
            ),
        ]
    )

    cm_df = imputer.fit_transform(cm_df)
    cm_df

Annotate :term:Features<features>` with Identified Compounds
************************************************************
``in``: :py:class:`~.ConsensusMap` DataFrame without identifications (cm_df) and :py:class:`~.AccurateMassSearch` DataFrame (ams_df)

``out``: :py:class:`~.ConsensusMap` DataFrame with new identifications column (id_df)

.. code-block:: python
    :linenos:

    id_df = cm_df

    id_df["identifications"] = pd.Series(["" for x in range(len(id_df.index))])

    for rt, mz, description in zip(
        ams_df["retention_time"],
        ams_df["exp_mass_to_charge"],
        ams_df["description"],
    ):
        indices = id_df.index[
            np.isclose(id_df["mz"], float(mz), atol=1e-05)
            & np.isclose(id_df["RT"], float(rt), atol=1e-05)
        ].tolist()
        for index in indices:
            if description != "null":
                id_df.loc[index, "identifications"] += str(description) + ";"
    id_df["identifications"] = [
        item[:-1] if ";" in item else "" for item in id_df["identifications"]
    ]
    id_df.to_csv(os.path.join(files, "result.tsv"), sep="\t", index=False)
    id_df

Visualize :term:`Consensus Features<consensus features>` with Identifications
*****************************************************************************

.. code-block:: python
    :linenos:

    fig = px.scatter(id_df, x="RT", y="mz", hover_name="identifications")
    fig.update_layout(title="Consensus features with identifications (hover)")
    fig.show()