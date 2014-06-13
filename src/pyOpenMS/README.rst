
.. contents:: **Table of Contents**

------------
Introduction
------------

This package contains Python bindings for a large part of the OpenMS library
(http://www.open-ms.de) for mass spectrometry based proteomics.  It thus
provides providing facile access to a feature-rich, open-source algorithm
library for mass-spectrometry based proteomics analysis. These Python bindings
allow raw access to the data-structures and algorithms implemented in OpenMS,
specifically those for file access (mzXML, mzML, TraML, mzIdentML among
others), basic signal processing (smoothing, filtering, de-isotoping and
peak-picking) and complex data analysis (including label-free, SILAC, iTRAQ and
SWATH analysis tools).

The pyOpenMS package runs - like OpenMS - on Windows, Linux and OSX.

------------
Installation
------------

We provide binary packages for Python 2.6 and 2.7 on Windows (64 bit) 
and Linux 64 bit which makes the installation very straightforward with pip.
For other platforms, please refer to the compilation instructions.

Binary installation
===================

The current binaries require numpy **1.7.x**.
As we distribute the package as binary eggs, you have to use *easy_install*,
installing with *pip* does not work::

    $ easy_install pyopenms


Source installation
===================

Download the latest OpenMS source from SVN (following `the OpenMS documentation`_), configure and build.

Install Qt and then start with the dependencies of OpenMS itself::

    $ svn co https://open-ms.svn.sourceforge.net/svnroot/open-ms/contrib
    $ cmake .

Now you have to install the dependencies of pyOpenMS:

- Install Python (2.6 or 2.7)
- Install numpy (On OSX, numpy should already be installed. On GNU/Linux there
  should be packages for numpy (e.g. python-numpy for Ubuntu/Debian). On
  Windows, you can install it from `Christoph Gohlkes webpage`_).
- Install setuptools, see the `setuptools PyPI page`_ .
- Use setuptools to install pip, autowrap and nose::

   $ easy_install pip
   $ pip install autowrap
   $ pip install nose

- Configure and build pyOpenMS::

    $ svn co https://open-ms.svn.sourceforge.net/svnroot/open-ms/OpenMS
    $ cmake -DPYOPENMS=ON .
    $ make pyopenms_bdist_egg

This should build a file like *pyopenms-1.10.1-py2.7-linux-x86_64.egg* the
folder *./pyOpenMS/dist* of your build directory which you can distribute
or install it from there::

    $ cd pyOpenMS/dist
    $ easy_install pyopenms-1.10.1-py-2.7-linux-x86_64.egg

------------
Testing
------------

pyOpenMS provides unittests, they are found under ./pyOpenMS/tests/ and can be
executed using nosetests::

    $ python run_nose.py

------------
License
------------

pyOpenMS is published under the 3-clause BSD licence, see ./pyOpenMS/License.txt

-------------
Documentation
-------------

pyOpenMS follows the `OpenMS
documentation <http://www-bs2.informatik.uni-tuebingen.de/services/OpenMS/OpenMS-release/html/classes.html>`_ very closely. Additionally, there is also a `pyOpenMS
Manual <http://proteomics.ethz.ch/pyOpenMS_Manual.pdf>`_ available. The online
manual contains a complete record of every wrapped class and function while the
documentation of the corresponding class or function can be inferred from the
OpenMS online documentation.



.. _the OpenMS documentation: http://www-bs2.informatik.uni-tuebingen.de/services/OpenMS/OpenMS-release/html/index.html
.. _Christoph Gohlkes webpage: http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
.. _setuptools PyPI page: https://pypi.python.org/pypi/setuptools

