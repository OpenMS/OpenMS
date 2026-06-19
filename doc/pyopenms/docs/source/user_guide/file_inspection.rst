File Inspection
===============

The :py:class:`~.FileInfo` class is the library-level equivalent of the
``FileInfo`` TOPP command-line tool. It auto-detects the file type, loads the
file, and returns a structured :py:class:`~.FileInfo.Result` that exposes all
metadata, range dimensions, statistics, and per-type details (peak, feature,
identification, FASTA, mzTab). The same result can also be rendered as the
exact human-readable text or TSV that the ``FileInfo`` CLI tool emits.

Quick Start
-----------

.. code-block:: python

    import pyopenms as oms

    r = oms.FileInfo().run_all("data.mzML")

    # Structured typed fields
    print(r.meta.file_type_name)          # "mzML"
    print(r.peak.num_spectra)             # e.g. 1234
    print(r.peak.spectra_per_ms_level)    # {1: 400, 2: 834}

    # Same human-readable text the CLI tool writes
    print(r.text)

``run_all`` switches on all content metrics (meta, processing, statistics). Use
``run`` with an :py:class:`~.FileInfo.Options` object to control exactly which
sections are computed.

Working with Different File Types
----------------------------------

**Peak data (mzML)**

.. code-block:: python

    r = oms.FileInfo().run_all("experiment.mzML")

    if r.peak:
        print("Spectra:", r.peak.num_spectra)
        print("Chromatograms:", r.peak.num_chromatograms)
        print("MS levels:", r.peak.ms_levels)
        # Activation methods as flat list of (ms_level, name, count) tuples
        for level, method, count in r.peak.activation_methods_flat():
            print(f"  MS{level} {method}: {count}")

    # Global m/z and RT ranges
    mz = r.ranges.combined.mz
    if mz.present:
        print(f"m/z: {mz.min:.4f} – {mz.max:.4f}")

**Feature data (featureXML / consensusXML)**

.. code-block:: python

    r = oms.FileInfo().run_all("features.featureXML")
    if r.feature:
        print("Features:", r.feature.num_features)
        print("Unassigned IDs:", r.feature.unassigned_ids)

    r2 = oms.FileInfo().run_all("consensus.consensusXML")
    if r2.feature and r2.feature.is_consensus:
        for col in r2.feature.map_columns:
            print(col.filename, col.size)

**Identification data (idXML)**

.. code-block:: python

    r = oms.FileInfo().run_all("search.idXML")
    if r.ident:
        print("Runs:", r.ident.num_runs)
        print("Peptide hits:", r.ident.peptide_hits)
        print("Proteins:", r.ident.protein_hits)

Controlling Which Sections Are Computed
----------------------------------------

Use :py:class:`~.FileInfo.Options` to enable optional (and potentially
expensive) sections:

.. code-block:: python

    opt = oms.FileInfo.Options()
    opt.statistics = True    # summary statistics per intensity / RT / m/z
    opt.meta = True          # instrument / sample / contact metadata
    opt.processing = True    # processing steps recorded in the file

    r = oms.FileInfo().run("experiment.mzML", opt)

    if "-- Statistics --" in r.text:
        print("Statistics were computed")

    for step in r.processing:
        print(step.software_name, step.software_version)

TSV Output
----------

Every result carries a machine-readable TSV mirror of the human-readable report.
This is identical to what the ``FileInfo -out_tsv`` flag emits:

.. code-block:: python

    r = oms.FileInfo().run_all("experiment.mzML")
    tsv_lines = r.tsv.splitlines()
    for line in tsv_lines[:10]:
        print(line)

    # Static renderers are also available
    text = oms.FileInfo.to_text(r)
    tsv  = oms.FileInfo.to_tsv(r)

Forced File Type
----------------

When the file extension is ambiguous or non-standard, pass
:py:attr:`~.FileInfo.Options.forced_type` to select the parse branch:

.. code-block:: python

    opt = oms.FileInfo.Options()
    opt.forced_type = oms.FileTypes.FileType.FEATUREXML

    r = oms.FileInfo().run("data_without_extension", opt)
    print(r.meta.file_type_name)   # "featureXML"

Range Information
-----------------

All file types report :py:class:`~.FileInfo.Ranges` with per-dimension min/max
intervals. MSExperiment files additionally expose per-MS-level and
chromatogram ranges:

.. code-block:: python

    r = oms.FileInfo().run_all("experiment.mzML")
    ranges = r.ranges

    # Per-MS-level ranges
    if ranges.is_experiment:
        for level, rs in ranges.per_ms_level.items():
            mz = rs.mz
            if mz.present:
                print(f"MS{level} m/z: {mz.min:.4f} – {mz.max:.4f}")

        chrom_rt = ranges.chromatograms.rt
        if chrom_rt.present:
            print(f"Chromatogram RT: {chrom_rt.min:.1f} – {chrom_rt.max:.1f} s")
