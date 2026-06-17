Summary
=======

pyOpenMS is an open-source Python library for :term:`mass spectrometry<Mass Spectrometry>`, specifically for the analysis of proteomics
and metabolomics data in Python. pyOpenMS implements a set of Python bindings to the OpenMS library for computational
:term:`mass spectrometry<Mass Spectrometry>` and is available for Windows, Linux and macOS.

pyOpenMS provides functionality that is commonly used in computational mass
spectrometry. The pyOpenMS package contains Python bindings for a large part of the `OpenMS <https://openms.de>`_
library for :term:`mass spectrometry<Mass Spectrometry>` based proteomics. It thus provides easy access to
a feature-rich, open-source algorithm library for mass-spectrometry based proteomics analysis.

pyOpenMS facilitates the execution of common tasks in proteomics (and other fields of mass spectrometry) such as

- File handling (:term:`mzXML`, :term:`mzML`, TraML, mzTab, :term:`FASTA`, pepxml, protxml, mzIdentML among others)
- Chemistry (mass calculation, peptide fragmentation, isotopic abundances)
- Signal processing (smoothing, filtering, de-isotoping, retention time correction and peak-picking)
- Identification analysis (including peptide search, PTM analysis, cross-linked analytes, FDR control, RNA oligonucleotide search and small molecule search tools)
- Quantitative analysis (including label-free, metabolomics, :term:`SILAC`, :term:`iTRAQ` and :term:`SWATH`/DIA analysis tools)
- Chromatogram analysis (chromatographic peak picking, smoothing, elution profiles and peak scoring for :term:`SRM`/MRM/PRM/:term:`SWATH`/DIA data)
- Interaction with common tools in proteomics and metabolomics:

  - Search engines such as Comet, Mascot, MSGF+, MSFragger, SpectraST, Sage
  - Post-processing tools such as Percolator, MSstats, Epiphany
  - Metabolomics tools such as SIRIUS, CSI:FingerId


User Guide
==========

Information about installing and using pyopenms.

.. toctree::
    :maxdepth: 2

    user_guide/index


API documentation
=================

Documentation on modules, functions, classes of this package.

.. toctree::
    :maxdepth: 2
   
    apidocs/index


Community and contribution guide
================================

Information about the community behind this package and how you can contribute.

.. toctree::
    :maxdepth: 2
   
    community/index


.. toctree::
    :hidden:

    C++ library <https://openms.de/documentation/html/index.html>
