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


Inspecting File Contents (FileInfo)
-------------------------------------------------------

:py:class:`~.FileInfo` is the library-level equivalent of the ``FileInfo``
TOPP command-line tool.  It loads any OpenMS-readable file, extracts
structured information about its contents, and can reproduce the exact
human-readable or TSV output of the CLI tool — all from Python.

**Quick inspection with runAll()**

:py:meth:`~.FileInfo.run_all` enables all content metrics (metadata,
processing history, statistics) with a single call:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    fi = oms.FileInfo()
    r = fi.run_all("data.mzML")

    # file type is auto-detected
    print(r.meta.file_type_name)          # e.g. "mzML"

    # peak-file specifics (None for non-peak files)
    if r.peak is not None:
        print("spectra:", r.peak.num_spectra)
        print("MS levels:", r.peak.ms_levels)
        for ms_level, name, count in r.peak.activation_methods_flat():
            print(f"  MS{ms_level} {name}: {count}")

    # RT / m/z / intensity ranges
    print("RT range:", r.ranges.combined.rt.min, "–", r.ranges.combined.rt.max)

    # human-readable output identical to the FileInfo CLI
    print(r.text)

**Feature and identification files**

The same API works for ``featureXML``, ``consensusXML``, and ``idXML``:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    r_feat = oms.FileInfo().run_all("data.featureXML")
    if r_feat.feature is not None:
        print("features:", r_feat.feature.num_features)
        print("charges:", r_feat.feature.charges)

    r_id = oms.FileInfo().run_all("data.idXML")
    if r_id.ident is not None:
        print("peptide hits:", r_id.ident.peptide_hits)
        print("search engines:", r_id.ident.search_engines)

**Selective computation with Options**

By default :py:meth:`~.FileInfo.run` skips expensive operations
(validation, corruption check, detailed spectrum listing).  Pass an
:py:class:`~.FileInfo.Options` object to opt in:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    opt = oms.FileInfo.Options()
    opt.meta        = True   # instrument / sample / contact block (-m)
    opt.processing  = True   # data-processing history (-p)
    opt.statistics  = True   # intensity / charge statistics (-s)
    opt.validate    = True   # schema / semantic validation (-v)

    r = oms.FileInfo().run("data.mzML", opt)

    if r.experiment_meta is not None and r.experiment_meta.present:
        print("instrument:", r.experiment_meta.instrument_name)
    print("-- Statistics --" in r.text)   # True

**Rendering to text or TSV**

The cached CLI output is always available as ``r.text`` and ``r.tsv``.
The static helpers :py:meth:`~.FileInfo.to_text` and
:py:meth:`~.FileInfo.to_tsv` are provided for API symmetry:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    r = oms.FileInfo().run_all("data.featureXML")

    # write the human-readable report to a file
    with open("report.txt", "w") as f:
        f.write(oms.FileInfo.to_text(r))

    # write the machine-readable TSV
    with open("report.tsv", "w") as f:
        f.write(oms.FileInfo.to_tsv(r))
