Glossary
========

A glossary of common terms used throughout OpenMS documentation.

```{glossary}
:sorted:

Liquid chromatography
  An analytical technique used to separate molecules.

mass spectrometry
  An analytical technique used to identify and quantify molecules of interest.

LC-MS
  Liquid Chromatography-Mass Spectrometry.

peptides
  A short chain of amino acids.

FASTA format
  A text-based format for representing nucleotide or amino acid sequences.

octadecyl (C18)
  An alkyl radical C(18)H(37) derived from an octadecane by removal of one hydrogen atom.

Mass
  Mass is a measure of the amount of matter that an object contains. In comparison to often used term weight, which is
  a measure of the force of gravity on that object.

ion
  Any {term}`atom` or group of atoms that bears one or more positive or negative electrical charges. Positively charged are cations, negavtively charged anions.

electrospray ionization (ESI)
  A technique used in mass spectrometry to produce ions.

MS^2
  See {term}`MS/MS`.

atom
  An atom is the smallest unit of ordinary matter that forms a chemical element.

electrospray ionization
  A technique used in mass spectrometry to produce ions using an electrospray in which a high voltage is applied to a liquid to create an {term}`aerosol`.

aerosol
  An aerosol is a suspension of fine solid particles or liquid droplets in air or another gas.

time-of-flight (TOF)
  A measurement of the time taken by an object, particle of wave (be it acoustic, electromagnetic, e.t.c) to travel a distance through a medium.

quadrupole mass filters
  A mass filter allowing one mass channel at a time to reach the detector as the mass range is scanned.

orbitrap analyzers
  In mass spectrometry, an ion trap mass analyzer consisting of an outer barrel-like electrode and a coaxial inner
  spindle-like electrode that traps ions in an orbital motion around the spindle.
  A high resoltion mass spectrometry analyzer.

MS(1)
  First stage to get a spectra. A sample is injected into the mass spectrometer, ionized, accelerated and analyzed by mass spectrometry.

MS(2)
  Ions from MS1 spectra are then selectively fragmented and analyzed by a second stage of mass spectrometry (MS2) to
  generate the spectra for the ion fragments.

MS/MS
  Tandem mass spectrometry, MS^2^, a technique where two or more mass analyzers are coupled together using an additional reaction step to increase their abilities to analyse chemical samples.

collision-induced dissociation (CID)
  A mass spectrometry technique to induce fragmentation of selected ions in the gas phase. Also known as Collision
  induced dissociation.

TOPP
  The OpenMS Pipeline, see TOPP Tools.

MSGFPlusAdapter
  Adapter for the MS-GF+ protein identification (database search) engine. More information is available in the [OpenMS API reference documentation](https://openms.de/current_doxygen/html/TOPP_MSGFPlusAdapter.html).

LuciphorAdapter
  Adapter for the LuciPHOr2: a site localisation tool of generic post-translational modifications from tandem mass
  spectrometry data. More information is available in the [OpenMS API reference documentation](https://openms.de/current_doxygen/html/TOPP_LuciphorAdapter.html).

pyOpenMS
  pyOpenMS is an open-source Python library for mass spectrometry, specifically for the analysis of proteomics and
  metabolomics data in Python. For pyOpenMS documentaion visit [this](https://pyopenms.readthedocs.io/en/latest/) link.

TOPP Tools
  OpenMS provides a number of programs, called TOPP tools, that process mass spectrometry data. More information on TOPP tools can be found in the [OpenMS API reference documentation](https://openms.de/current_doxygen/html/TOPP_documentation.html).

TOPP tool
  see {term}`TOPP Tools`

TOPPView
  TOPPView is a viewer for MS and HPLC-MS data which ships with OpenMS. More information is available in [TOPPView documentation](/getting-started/visualize-with-openms.md).

[Nightly Snapshot](https://openms.de/current_doxygen/html/index.html)
  Untested installers and containers are known as the nightly snapshot.

proteomics
  Proteomics is the large-scale study of proteins.

proteins
  Proteins are vital parts of living organisms, with many functions, for example composing the structural fibers of
  muscle to the enzymes that catalyze the digestion of food to synthesizing and replicating DNA.

Mascot
  A so-called search engine: It identifies peptide sequences from MS/MS spectra.

HPLC-MS
  Data produced by High performance liquid chromatography (HPLC) separates components of a mixture, whereas mass
  spectrometry (MS) offers the detection tools to identify them.

mzML
  The mzML format is an open, XML-based format for mass spectrometer output files, developed with the full participation
  of vendors and researchers in order to create a single open format that would be supported by all software.

mzData
  mzData was the first attempt by the Proteomics Standards Initiative (PSI) from the Human Proteome Organization (HUPO)
  to create a standardized format for Mass Spectrometry data.[7] This format is now deprecated, and replaced by mzML.

mzXML
  mzXML is an open data format for storage and exchange of mass spectroscopy data, developed at the SPC/Institute for
  Systems Biology.

spectra
  Plural of spectrum.

mass spectrum
  A mass spectrum is a plot of the ion signal as a function of the mass-to-charge ratio. A mass spectrum is produced by a single mass spectrometry run. These spectra are used to determine the elemental or isotopic signature of a sample, the masses of particles and of molecules, and to elucidate the chemical identity or structure of molecules and other chemical compounds. OpenMS represents a one dimensional mass spectrum using the class [MSSpectrum](https://openms.de/current_doxygen/html/classOpenMS_1_1MSSpectrum.html). 

m/z
  mass to charge ratio.

retention time
  retention time (RT) in liquid chromatography, is the time it takes for a separated analyte to move through the stationary phase.

ProteoWizard
  ProteoWizard is a set of open-source, cross-platform tools and libraries for proteomics data analyses. It provides a framework for unified mass spectrometry data file access and performs standard chemistry and LCMS dataset computations.

PepNovo
  PepNovo is a de novo sequencing algorithm for {term}`MS/MS` {term}`spectra`.

de novo peptide sequencing
  A peptide’s amino acid sequence is inferred directly from the precursor peptide mass and tandem mass spectrum ({term}`MS/MS` or {term}`MS^3`) fragment ions, without comparison to a reference proteome.

TOPPAS
  An graphical user interface (GUI), which is shipped with OpenMS, to create and execute worflows using {term}`TOPP tools`; see [Workflow Editor](/getting-started/workflows.rst).

chromatogram
  A two-dimensional plot that describes the amount of analyte eluted from a chromatography versus the analyte's retention time. OpenMS represents a chromatogram using the class [MSChromatogram](https://openms.de/current_doxygen/html/structOpenMS_1_1Interfaces_1_1Chromatogram.html)

KNIME
  An advanced workflow editor which OpenMS provides a plugin for; see [Workflow Editor](/getting-started/workflows.rst).

Nextflow
  Script/DSL-based workflow language, executor and utilities; see [Workflow Editor](/getting-started/workflows.rst).

Galaxy
  Server and browser-based interactive workflow editor and runner; see [Workflow Editor](/getting-started/workflows.rst).

SILAC
  Stands for 'Stable isotope labeling using amino acids in cell culture'.

iTRAQ
  Stands for 'Isobaric tags for relative and absolute quantitation'.

TMT
  Tandem Mass Tag (TMT) is a mass spectrometry based system designed to identify and quantify proteins in different samples.

SRM
  Selected reation monitoring is a mass spectrometry technique for small molecule analysis.

SWATH
  Stands for 'Sequential acquisition of all theoretical fragment ion spectra'.

OpenMS API
  An interface that allows developers to use OpenMS core library classes and methods. 

RT
  Retention time.

MS
  Mass Spectrometry

MS^3
  Multi-stage Mass Spectrometry

feature
  An LC-MS feature represents the combined isotopic mass traces of a detected chemical compound. The chromatographic peak shape of a feature is defined by the interaction of the analyte with the LC column. Each feature contains information on retention time, mass-to-charge ratio, intensity and overall quality. OpenMS represents a feature using the class [Feature](https://openms.de/current_doxygen/html/classOpenMS_1_1Feature.html).

feature map
  A feature map is a collection of features identified in a mass spectrum from a single experiment. One feature map can contain many features. OpenMS represents a feature map using the class [FeatureMap](https://openms.de/current_doxygen/html/classOpenMS_1_1FeatureMap.html).

consensus feature
  Features from replicate experiments with similar retention times and m/z values are linked and considered a consensus feature. A consensus feature contains information on the common retention time and m/z values as well as intensities for each sample. OpenMS represents a consensus feature using the class [ConsensusFeature](https://openms.de/current_doxygen/html/classOpenMS_1_1ConsensusFeature.html).

consensus map
  A consensus map is a collection of {term}`consensus features <consensus feature>` identified from mass spectra across replicate experiments. One consensus map can contain many consensus features. OpenMS represents a consensus map using the class [ConsensusMap](https://openms.de/current_doxygen/html/classOpenMS_1_1ConsensusMap.html).

peak
  A single raw data point in a chromatogram or a mass spectrum. OpenMS represents a peak in a chromatogram using the class [ChromatogramPeak](https://openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramPeak.html). OpenMS represents a single, one-dimensional peak in a mass spectrum using the class [PeakID](https://openms.de/current_doxygen/html/classOpenMS_1_1Peak1D.html)

MSExperiment
  An OpenMS class used to represent a single mass spectrometry run. [Read the documentation for further information](https://openms.de/current_doxygen/html/classOpenMS_1_1MSExperiment.html).
```
