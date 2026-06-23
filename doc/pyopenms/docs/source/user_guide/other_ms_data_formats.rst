Other MS Data Formats
=============================

.. _anchor-other-id-data:

Identification Data (idXML, mzIdentML, pepXML, protXML)
-------------------------------------------------------

You can store and load identification data from an `idXML` file as follows:

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve
    import pyopenms as oms

    gh = gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/IdXMLFile_whole.idXML", "test.idXML")
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load("test.idXML", protein_ids, peptide_ids)
    oms.IdXMLFile().store("test.out.idXML", protein_ids, peptide_ids)

You can store and load identification data from an `mzIdentML` file as follows:

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve

    gh = gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/MzIdentML_3runs.mzid", "test.mzid")
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.MzIdentMLFile().load("test.mzid", protein_ids, peptide_ids)
    oms.MzIdentMLFile().store("test.out.mzid", protein_ids, peptide_ids)
..  # alternatively: -- dont do this, doesnt work
    identifications = oms.Identification()
    oms.MzIdentMLFile().load("test.mzid", identifications)

You can store and load identification data from a TPP `pepXML` file as follows:

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve

    gh = gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/PepXMLFile_test.pepxml", "test.pepxml")
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.PepXMLFile().load("test.pepxml", protein_ids, peptide_ids)
    oms.PepXMLFile().store("test.out.pepxml", protein_ids, peptide_ids)


You can load (storing is not supported) identification data from a TPP `protXML` file as follows:

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve

    gh = gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/ProtXMLFile_input_1.protXML", "test.protXML")
    protein_ids = oms.ProteinIdentification()
    peptide_ids = oms.PeptideIdentification()
    oms.ProtXMLFile().load("test.protXML", protein_ids, peptide_ids)
    # storing protein XML file is not yet supported
..    ProtXMLFile().store("test.out.protXML", protein_ids, peptide_ids, "doc_id_42")


Note how each data file produces two vectors of type :py:class:`~.ProteinIdentification`
and :py:class:`~.PeptideIdentification` which also means that conversion between two data
types is trivial: load data from one data file and use the storage function of
the other file.

Quantiative Data (featureXML, consensusXML)
-------------------------------------------------------

OpenMS stores quantitative information in the internal ``featureXML`` and
``consensusXML`` attributes. The ``featureXML`` format is used to store
quantitative data from a single :term:`LC-MS/MS` run while the ``consensusXML`` is used
to store quantitative data from multiple :term:`LC-MS/MS` runs. These can be accessed
as follows:

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve

    gh = gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(
        gh + "/src/data/FeatureFinderCentroided_1_output.featureXML",
        "test.featureXML",
    )
    features = oms.FeatureMap()
    oms.FeatureXMLFile().load("test.featureXML", features)
    oms.FeatureXMLFile().store("test.out.featureXML", features)

and for ``consensusXML``

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve

    gh = gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(
        gh + "/src/data/ConsensusXMLFile_1.consensusXML", "test.consensusXML"
    )
    consensus_features = oms.ConsensusMap()
    oms.ConsensusXMLFile().load("test.consensusXML", consensus_features)
    oms.ConsensusXMLFile().store("test.out.consensusXML", consensus_features)


.. PyOpenMS also also supports mzQuantML, however this format is currently work in
.. progress and should not be considered stable.
.. 
.. .. code-block:: python
.. 
..     msquant = MSQuantifications()
..     msquant.addConsensusMap(consensus_features)
..     MzQuantMLFile().store("file.mzquant", msquant)
..

Transition data (TraML)
-------------------------------------------------------

The TraML data format allows you to store transition information for targeted
experiments (:term:`SRM` / MRM / PRM / DIA).

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
    urlretrieve(gh + "/src/data/ConvertTSVToTraML_output.TraML", "test.TraML")
    targeted_exp = oms.TargetedExperiment()
    oms.TraMLFile().load("test.TraML", targeted_exp)
    oms.TraMLFile().store("test.out.TraML", targeted_exp)


Inspecting MS File Information (FileInfo)
-------------------------------------------------------

:py:class:`~.FileInfo` is the library-level equivalent of the OpenMS
``FileInfo`` command-line tool. It auto-detects the file type and returns a
structured :py:class:`~.FileInfo.Result` with ranges, per-type counts, and
statistics — without requiring a subprocess call.

Basic usage: load a file and inspect common fields:

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve
    import pyopenms as oms

    gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/BSA1.mzML", "BSA1.mzML")

    fi = oms.FileInfo()

    # run_all() enables meta/processing/statistics in one call
    r = fi.run_all("BSA1.mzML")

    print("File type:", r.meta.file_type_name)
    print("Number of spectra:", r.peak.num_spectra)
    print("MS levels:", r.peak.ms_levels)
    print("RT range:", r.ranges.spectra_overall.rt.min, "–", r.ranges.spectra_overall.rt.max)

    # Activation methods as (ms_level, method_name, count) tuples
    for ms_level, name, count in r.peak.activation_methods_flat():
        print(f"  MS{ms_level} {name}: {count}")

For a feature map or consensus map:

.. code-block:: python
    :linenos:

    from urllib.request import urlretrieve
    import pyopenms as oms

    gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/BSA1_F1_idmapped.featureXML", "test.featureXML")
    urlretrieve(gh + "/src/data/ConsensusXMLFile_1.consensusXML", "test.consensusXML")

    r = oms.FileInfo().run_all("test.featureXML")
    print("Features:", r.feature.num_features)
    print("Is consensus:", r.feature.is_consensus)

    r2 = oms.FileInfo().run_all("test.consensusXML")
    for col in r2.feature.map_columns:
        print(col.filename, col.size)

You can opt into specific analyses via :py:class:`~.FileInfo.Options`:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    opt = oms.FileInfo.Options()
    opt.statistics = True   # equivalent to FileInfo CLI -s flag
    opt.validate = True     # equivalent to -v
    r = oms.FileInfo().run("BSA1.mzML", opt)

    # The same human-readable text the FileInfo CLI would write to -out
    print(r.text)
    # TSV output (equivalent to -out_tsv)
    print(r.tsv)
    # Validation result
    print("Valid:", r.validation.valid)
