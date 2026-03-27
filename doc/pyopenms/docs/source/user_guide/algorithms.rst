Algorithms 
==========

Many signal processing algorithms follow a similar pattern in OpenMS.

.. code-block:: output

  algorithm = NameOfTheAlgorithmClass()
  exp = MSExperiment()
  
  # populate exp, for example load from file
  # ...
  
  # run the algorithm on data
  algorithm.filterExperiment(exp)

In many cases, the processing algorithms have a set of parameters that can be
adjusted. These are accessible through ``Algorithm.getParameters()`` and yield a
:py:class:`~.Param` object (see `Parameter handling <parameter_handling.html>`_) which can
be manipulated. After changing parameters, one can use ``Algorithm.setParameters()`` to
propagate the new parameters to the algorithm:

.. code-block:: output

  algorithm = NameOfTheAlgorithmClass()
  param = algorithm.getParameters()
  param.setValue("algo_parameter", "new_value")
  algorithm.setParameters(param)

  exp = MSExperiment()

  # populate exp, for example load from file
  # ...

  algorithm.filterExperiment(exp)

Since they work on a single :py:class:`~.MSExperiment` object, little input is needed to
execute a filter directly on the data. Examples of filters that follow this
pattern are :py:class:`~.GaussFilter`, :py:class:`~.SavitzkyGolayFilter` as well as the spectral filters
:py:class:`~.NLargest`, :py:class:`~.Normalizer`, :py:class:`~.SpectraMerger`, :py:class:`~.ThresholdMower`, :py:class:`~.WindowMower`.

Using the same example file as before, we can execute a :py:class:`~.GaussFilter` on our test data as follows:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master"
    urlretrieve(gh + "/src/data/tiny.mzML", "test.mzML")

    gf = oms.GaussFilter()
    exp = oms.MSExperiment()
    oms.MzMLFile().load("test.mzML", exp)
    gf.filterExperiment(exp)
    # oms.MzMLFile().store("test.filtered.mzML", exp)

