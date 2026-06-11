Common Workflow Language
========================

The Common Workflow Language is a set of free and open standards for describing
command-line tool and workflows made from those tool descriptions.

OpenMS tools can self-generate a CWL tool description using the `-write_cwl` option,
provided the build was compiled with TDL support enabled.

.. note::

   As of OpenMS 3.6, TDL (Tool Description Library) support is **disabled by default**
   to reduce build dependencies. CWL generation requires rebuilding OpenMS with
   ``-DENABLE_TDL=ON``:

   .. code-block:: bash

      cmake -DENABLE_TDL=ON ...

   Official binary releases include TDL support.

For example: `PeakPickerHiRes -write_cwl PeakPickerHiRes.cwl`
