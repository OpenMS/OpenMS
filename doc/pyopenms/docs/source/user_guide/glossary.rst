OpenMS Glossary
===============

A glossary of common terms used throughout OpenMS documentation.

.. glossary::
    :sorted:

    C18
        Octadecyl (C18) is an alkyl radical C(18)H(37) derived from an octadecane by removal of one hydrogen atom.

    CID
    collision-induced dissociation
        Collision-induced dissociation is an MS technique to induce fragmentation of selected ions in the gas phase,
        which are subjected to a subsequent measurement (see :term:`MS2`).

    consensus features
    consensus feature
        Features from replicate experiments with similar retention times and m/z values are linked and considered a consensus feature.
        A consensus feature contains information on the common retention time and m/z values as well as intensities for each sample.
        OpenMS represents a consensus feature using the class `ConsensusFeature
        <https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/classOpenMS_1_1ConsensusFeature.html>`_.

    consensus maps
    consensus map
        A consensus map is a collection of :term:`consensus features` identified from mass spectra across replicate experiments,
        usually by combining multiple :term:`feature maps`.
        One consensus map usually contains many consensus features. OpenMS represents a consensus map using 
        the class `ConsensusMap <https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/classOpenMS_1_1ConsensusMap.html>`_.

    de novo
    de novo peptide sequencing
        A peptide's amino acid sequence is inferred directly from the precursor peptide mass and tandem
        mass spectrum (:term:`MS2` or :term:`MS3`) fragment ions, without comparison to a reference proteome.

    electrospray ionization
    ESI
        Electrospray ionization (ESI) is a technique used in MS to produce ions.
        ESI usually gives rise to peptides of charge 2 or higher.
        ESI is usually coupled to Orbitrap and FTICR instruments.

    FASTA
        A text-based file format for representing nucleotide or amino acid sequences.

    feature
    features
        A feature, in the OpenMS terminology, subsumes all m/z signals originating from a single compound at a 
        certain charge state. This includes the isotope pattern and usually spans multiple spectra in retention time (the elution profile).
  
    feature finder
    FeatureFinder
        A FeatureFinder in pyOpenMS creates a :term:`feature map` using the MS1 spectra from a :term:`peak map` (e.g. from an mzML file) as input.

    feature maps
    feature map
        A feature map is a collection of :term:`feature`\ s identified from a single experiment.
        One feature map usually contains many features. OpenMS represents a feature map using the class
        `FeatureMap <https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/classOpenMS_1_1FeatureMap.html>`_.

    high performance liquid chromatography
    HPLC
        In high performance liquid chromatography (HPLC), analytes are dissolved in a pressurized solvent (mobile phase)
        and pumped through a solid adsorbent material (stationary phase) packed into a
        capillary column. Physicochemical properties of the analyte determine how strongly it
        interacts with the stationary phase.

    iTRAQ
        Isobaric tags for relative and absolute quantitation (iTRAQ) is a MS based multiplexing technique designed to
        identify and quantify proteins from different samples in one single measurement.

    KNIME
        An advanced workflow editor which OpenMS provides a plugin for.

    LC-MS
    LCMS
        :term:`Liquid chromatography<liquid chromatography>`-coupled mass spectrometry.

    LC-MS/MS
    liquid chromatography coupled :term:`tandem mass spectrometry`
        See :term:`LC-MS` and :term:`MS2`.

    liquid chromatography
    LC
        An analytical technique used to separate molecules of interest.

    LuciphorAdapter
        Adapter for the LuciPHOr2: a site localisation tool of generic post-translational modifications from tandem mass
        spectrometry data. More information is available in the `OpenMS API reference documentation
        <https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_LuciphorAdapter.html>`__.

    MALDI
        matrix-assisted laser desorption/ionization (MALDI) is an ionization technique that uses a laser energy-absorbing matrix to create ions.
        ESI usually gives rise to peptides of charge 1.
        MALDI is usually employed in combination with TOF instruments. 
    
    Mass Spectrometry
    MS
        An analytical technique to measure the mass over charge (m/z) ratio of ions along with their abundance. 
        This gives rise to a mass spectrum (with m/z on the x-axis and abundance on the y-axis).

    mass spectra
    mass spectrum
        A visual or numerical representation of a measurement from an MS instrument.
        A spectrum contains (usually many) pairs of mass-over-charge(m/z)+intensity values.

    MascotAdapter
        Used to identify peptides in :term:`MS2` spectra. Read more about this adapter in the `OpenMS API reference documentation
        <https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MascotAdapter.html>`__.

    MSGFPlusAdapter
        Adapter for the MS-GF+ protein identification (database search) engine. More information is available in the
        `OpenMS API reference documentation <https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_MSGFPlusAdapter.html>`__.

    MS1
        Mass spectra of a sample where only precursor ions (i.e. no fragment ions) can be observed.
        Usually MS1 spectra are recorded to select targets for MS2 fragmentation.

    MS2
    MS/MS
        Tandem MS is a technique where two or more mass analyzers are coupled together using an additional, usually destructive,
        reaction step to generate fragment ions which increases their abilities to analyse chemical samples.

    MS3
        Multi-stage MS.

    mzData
    mzdata
        mzData was the first attempt by the Proteomics Standards Initiative (PSI) from the Human Proteome Organization (HUPO)
        to create a standardized format for MS data. This format is now deprecated, and replaced by mzML.

    mzML
    mzml
        The mzML format is an open, XML-based format for mass spectrometer output files, developed by the Proteomics Standard Initiative (PSI)
        with the full participation of vendors and researchers in order to create a single open format that would be supported by all software.

    mzXML
    mzxml
        mzXML is an open data format for storage and exchange of mass spectroscopy data, developed at the SPC/Institute for
        Systems Biology. This format is now deprecated, and replaced by mzML.

    nightly snapshot
        Untested installers and containers which are created regularly between official releases and reflect the current development state.

    octadecyl
        See :term:`C18`.

    OpenMS API
        A C++ interface that allows developers to use OpenMS core library classes and methods.

    Orbitrap
    orbitrap
        In MS, an ion trap mass analyzer consisting of an outer barrel-like electrode and a coaxial inner
        spindle-like electrode that traps ions in an orbital motion around the spindle.
        An ultra-high resolution MS analyzer, capable of resolving fine-isotope structure.

    peak maps
    peak map
        A collection of mass spectra (and/or chromatograms), usually sorted by retention time. Can contain spectra of one or more MS levels (usually level 1 and 2).
  
    peptide-spectrum match
    PSM
    PSMs
        A peptide-spectrum match associates a peptide sequence (possibly including modifications) to an MS/MS spectrum.
        This usually involves comparing the mass spectra of peptide fragments generated from a digested protein sample with a database of predicted
        spectra. Alternatively, this can be done using :term:`de novo` techniques (without a database, just using observed spectra).

    PepNovo
        PepNovo is a de :term:`de novo peptide sequencing` algorithm for :term:`MS2` spectra.

    ProteoWizard
        ProteoWizard is a set of open-source, cross-platform tools and libraries for proteomics data analyses.
        It provides a framework for unified MS data file access and performs standard chemistry and LCMS dataset computations.

    quadrupole
        A low resolution MS analyzer.
        A mass filter allowing one mass channel at a time to reach the detector as the mass range is scanned.

    SILAC
    stable isotope labeling with amino acids in cell culture
        Stands for Stable isotope labeling using amino acids in cell culture.

    SRM
        Selected reaction monitoring (SRM) is a MS technique for targeted small molecule analysis.

    SWATH
        Sequential acquisition of all theoretical fragment ion spectra (SWATH) uses partially overlapping MS2 
        scans with wide isolation windows to capture all fragment ions in a data independent analysis (DIA).

    tandem mass spectrometry
        See :term:`MS2`.

    time-of-flight
    TOF
        Time-of-flight (TOF) is the time taken by an object, particle or wave (be it acoustic, electromagnetic, etc.)
        to travel a distance through a medium.
        TOF analyzers can obtain good, but not ultra-high resolution, such as an :term:`orbitrap`.

    TMT
        Tandem Mass Tag (TMT) is a MS based multiplexing technique designed to identify and 
        quantify proteins from different samples in one single measurement.

    TOPP
       'TOPP - The OpenMS PiPeline' is a pipeline for the analysis of HPLC-MS data. It consists of several small 
       applications that can be chained to create analysis pipelines tailored for a specific problem. See :term:`TOPP tools`.

    TOPPAS
        An assistant for GUI-driven :term:`TOPP` workflow design, build into OpenMS. 
        See `TOPPAS tutorial <https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPPAS_tutorial.html>` for details.

    TOPP tools
        OpenMS provides a number of applications (executable files) that are chainable in a pipeline/script and each process MS data.
        These tools are subdivided into different categories, such as 'File Handling' or 'Peptide Identification'.
        All :term:`TOPP` tools are described in the `OpenMS API reference documentation
        <https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_documentation.html>`__.

    TOPPView
        TOPPView is a viewer for MS and HPLC-MS data and shipped with every OpenMS release.