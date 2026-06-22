File Inspection
===============

The :py:class:`~.FileInfo` class is the library-level equivalent of the
``FileInfo`` TOPP tool.  It can inspect any OpenMS-readable file (mzML,
featureXML, idXML, FASTA, …) and return a structured
:py:class:`~.FileInfo.Result` object containing metadata, range information,
and file-type-specific statistics.

Unlike the TOPP command-line tool, :py:class:`~.FileInfo` can be embedded
directly in Python scripts and workflows without spawning a sub-process.


Basic Usage
***********

The simplest entry point is :py:meth:`~.FileInfo.run_all`, which enables the
full set of computation passes (metadata, processing history, statistics):

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/BSA1.mzML", "BSA1.mzML")

    result = oms.FileInfo().run_all("BSA1.mzML")

    # result.meta             – general file metadata (name, detected type, …)
    # result.ranges           – coordinate ranges across the whole file
    # result.peak             – MSExperiment-specific info (None for non-peak files)
    print("File type:", result.meta.file_type_name)
    print("Number of spectra:", result.peak.num_spectra)

For a lightweight pass that skips metadata / statistics computation use
:py:meth:`~.FileInfo.run` with default options:

.. code-block:: python
    :linenos:

    result_fast = oms.FileInfo().run("BSA1.mzML")
    print("MS levels:", result_fast.peak.ms_levels)


Accessing Peak File Information
*********************************

When the file contains mass spectra the ``peak`` field is a
:py:class:`~.FileInfo.PeakInfo` object:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/BSA1.mzML", "BSA1.mzML")

    result = oms.FileInfo().run_all("BSA1.mzML")

    if result.peak is not None:
        peak = result.peak
        print("Total spectra:         ", peak.num_spectra)
        print("MS levels present:     ", peak.ms_levels)            # list of ints, e.g. [1, 2]
        print("Spectra per MS level:  ", peak.spectra_per_ms_level) # dict int -> int
        print("Total peaks:           ", peak.total_peaks)
        print("Chromatograms:         ", peak.num_chromatograms)

        # Activation methods as a flat list of (ms_level, method_name, count) tuples
        for ms_level, method, count in peak.activation_methods_flat():
            print(f"  MS{ms_level} {method}: {count}")

The ``ranges`` field gives coordinate ranges regardless of file type.  Its
``combined`` sub-object covers all spectra and chromatograms together; each
dimension is a :py:class:`~.FileInfo.Range` with a ``present`` flag:

.. code-block:: python
    :linenos:

    ranges = result.ranges
    mz  = ranges.combined.mz
    rt  = ranges.combined.rt
    tic = ranges.combined.intensity
    if mz.present:
        print("m/z range: ", mz.min, "-", mz.max)
    if rt.present:
        print("RT range:  ", rt.min, "-", rt.max)
    if tic.present:
        print("Intensity: ", tic.min, "-", tic.max)


Accessing Feature, Identification and FASTA Information
*********************************************************

For **featureXML** files the ``feature`` field is a
:py:class:`~.FileInfo.FeatureInfo` object:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(
        gh + "/src/data/FeatureFinderCentroided_1_output.featureXML",
        "features.featureXML",
    )

    result = oms.FileInfo().run_all("features.featureXML")

    if result.feature is not None:
        feat = result.feature
        print("Number of features:", feat.num_features)
        print("Charges:           ", feat.charges)   # dict charge -> count

For **idXML / mzIdentML** files the ``ident`` field is a
:py:class:`~.FileInfo.IdentInfo` object:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/IdXMLFile_whole.idXML", "test.idXML")

    result = oms.FileInfo().run_all("test.idXML")

    if result.ident is not None:
        ident = result.ident
        print("Protein hits:           ", ident.protein_hits)
        print("Non-redundant proteins: ", ident.non_redundant_protein_hits)
        print("Peptide hits:           ", ident.peptide_hits)
        print("Non-redundant peptides: ", ident.non_redundant_peptides)

For **FASTA** files the ``fasta`` field is a :py:class:`~.FileInfo.FastaInfo`
object:

.. code-block:: python
    :linenos:

    # result = oms.FileInfo().run_all("sequences.fasta")
    # if result.fasta is not None:
    #     print("Sequences:      ", result.fasta.num_sequences)
    #     print("Total residues: ", result.fasta.total_residues)
    #     print("Median length:  ", result.fasta.length_stats.median)


Getting Human-Readable or TSV Output
***************************************

:py:meth:`~.FileInfo.to_text` converts a result to the same human-readable
string that the ``FileInfo`` TOPP tool prints to the console.  The same string
is also cached directly in ``result.text``:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    urlretrieve(gh + "/src/data/BSA1.mzML", "BSA1.mzML")

    result = oms.FileInfo().run_all("BSA1.mzML")

    # Human-readable summary (identical to the FileInfo TOPP tool -out output)
    print(result.text)                          # or: oms.FileInfo.to_text(result)

:py:meth:`~.FileInfo.to_tsv` / ``result.tsv`` return a tab-separated version
suitable for downstream parsing:

.. code-block:: python
    :linenos:

    with open("file_info.tsv", "w") as fh:
        fh.write(result.tsv)                    # or: oms.FileInfo.to_tsv(result)


Selective Computation with Options
************************************

For large files you may wish to skip expensive computation steps.
:py:meth:`~.FileInfo.run` accepts a :py:class:`~.FileInfo.Options` object that
controls which passes are executed:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    opts = oms.FileInfo.Options()
    opts.meta       = True    # include experiment/instrument/sample metadata (-m)
    opts.processing = True    # include processing history (-p)
    opts.statistics = False   # skip intensity statistics (-s)
    opts.validate   = False   # skip schema validation (-v)

    result = oms.FileInfo().run("BSA1.mzML", opts)

    if result.experiment_meta is not None and result.experiment_meta.present:
        print("Instrument:", result.experiment_meta.instrument_name)

Enabling only the passes you need can significantly reduce run time when a
pipeline requires just a subset of the information (for example, only the file
type and spectrum count for pre-filtering).
