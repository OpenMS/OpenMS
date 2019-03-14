OpenMS
====== 

[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg?style=flat-square)](http://opensource.org/licenses/BSD-3-Clause)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/93e71bad214f46d2a534ec92dbc2efc9)](https://www.codacy.com/app/OpenMS/OpenMS?utm_source=github.com&utm_medium=referral&utm_content=OpenMS/OpenMS&utm_campaign=badger)
[![Build Status](https://travis-ci.org/OpenMS/OpenMS.svg?branch=develop)](https://travis-ci.org/OpenMS/OpenMS)
[![Project Stats](https://www.openhub.net/p/open-ms/widgets/project_thin_badge.gif)](https://www.openhub.net/p/open-ms)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/OpenMS/OpenMS?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/openms/README.html)
[![documentation](https://codedocs.xyz/doxygen/doxygen.svg)](https://ftp.mi.fu-berlin.de/pub/OpenMS/develop-documentation/html/index.html)


<a href="http://www.openms.org/" target="_blank">OpenMS</a> 
is an open-source software C++ library for LC-MS data management and
analyses. It offers an infrastructure for rapid development of mass
spectrometry related software. OpenMS is free software available under the
three clause BSD license and runs under Windows, macOS, and Linux.

It comes with a vast variety of pre-built and ready-to-use tools for proteomics
and metabolomics data analysis (TOPPTools) as well as powerful 1D, 2D and 3D
visualization (TOPPView).

OpenMS offers analyses for various quantitation protocols, including label-free
quantitation, SILAC, iTRAQ, TMT, SRM, SWATH, etc.

It provides built-in algorithms for de-novo identification and database search,
as well as adapters to other state-of-the art tools like X!Tandem, Mascot,
OMSSA, etc. It supports easy integration of OpenMS built tools into workflow
engines like KNIME, Galaxy, WS-Pgrade, and TOPPAS via the TOPPtools concept and
a unified parameter handling via a 'common tool description' (CTD) scheme.

With pyOpenMS, OpenMS offers Python bindings to a large part of the OpenMS API
to enable rapid algorithm development. OpenMS supports the Proteomics Standard
Initiative (PSI) formats for MS data. The main contributors of OpenMS are
currently the Eberhard-Karls-Universität in Tübingen, the Freie Universität
Berlin, and the ETH Zürich.

Features
--------
- Core C++ library under three-clause BSD licence using modern C++11
- Python bindings to the C++ API through pyOpenMS
- Major community file formats supported (mzML, mzXML, mzIdentXML, pepXML etc)
- Over 150+ individual analysis tools (TOPP Tools), covering most MS and LC-MS data processing and mining tasks
- Powerful 1D, 2D and 3D visualization tools (TOPPView)
- Support for most MS identification and quantification workflows (targeted, DIA, label-free, isobaric and stable isotope)
- Support for all major platforms (Windows [XP, 7, 8, 10], macOS and Linux)

Documentation
-------------
Users and developers should start by reading the [OpenMS wiki](https://github.com/OpenMS/OpenMS/wiki) and consult the [current documentation](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/index.html) (or use the [latest release documentation](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/index.html)).
Documentation for the Python bindings pyOpenMS can be found on the [pyOpenMS online documentation](https://pyopenms.readthedocs.io).

Citation
--------
Please cite:

Rost HL, Sachsenberg T, Aiche S, Bielow C et al. OpenMS: a flexible open-source software platform for mass spectrometry data analysis. Nat Meth. 2016; 13, 9: 741-748. doi:10.1038/nmeth.3959.

The file [AUTHORS](AUTHORS) contains a list of all authors who worked on OpenMS.

Licence
-------
OpenMS is released under the [three clause BSD licence](LICENSE).
