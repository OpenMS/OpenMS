Query MSExperiment with MassQL
==============================

MassQL is a powerful, SQL-like query language for mass spectrometry data.
For further information visit the `MassQL documentation
<https://mwang87.github.io/MassQueryLanguage_Documentation/>`_.

MS data from a :py:class:`~.MSExperiment` can be exported to :term:`MS1` and :term:`MS2` dataframes, which can
be queried directly with the ``massql`` module.

**pyopenms.MSExperiment.get_massql_df()**
        Exports data from :py:class:`~.MSExperiment` to pandas DataFrames to be used with MassQL.
        
        Both dataframes contain the columns:
        'i': intensity of a peak
        'i_norm': intensity normalized by the maximun intensity in the spectrum
        'i_tic_norm': intensity normalized by the sum of intensities (TIC) in the spectrum
        'mz': mass to charge of a peak
        'scan': number of the spectrum
        'rt': retention time of the spectrum
        'polarity': ion mode of the spectrum as integer value (positive: 1, negative: 2)
        
        The :term:`MS2` dataframe contains additional columns:
        'precmz': mass to charge of the precursor ion
        'ms1scan': number of the corresponding :term:`MS1` spectrum
        'charge': charge of the precursor ion
        
        **Returns:**

        ms1_df : **pandas.DataFrame** 
        
        peak data of :term:`MS1` spectra

        ms2_df : **pandas.DataFrame** 
        
        peak data of :term:`MS2` spectra with precursor information

**Example:**

Load an example file into a :py:class:`~.MSExperiment` and get the :term:`MS1` and :term:`MS2` data frames for a MassQL query.

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from massql import msql_engine

    from urllib.request import urlretrieve

    url = "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master/src/data/"

    urlretrieve(url + "small.mzML", "small.mzML")

    # load MSExperiment
    exp = oms.MSExperiment()
    oms.MzMLFile().load("small.mzML", exp)

    # get MS1 and MS2 dataframes
    ms1_df, ms2_df = exp.get_massql_df()

    ms1_df.head()


.. csv-table:: ms1_df.head()
   :widths: 2 20 20 20 20 20 20 20
   :header: , i,  i_norm,   i_tic_norm,   mz,   scan, rt,   polarity

   0,  2105.75,  0.00455405,   0.000325626,  360.696,       1,  15.0015,           1
   1,  1172.47,  0.00253567,   0.000181306,  361.2,         1,  15.0015,           1
   2,  2287.57,  0.00494729,   0.000353743,  361.208,       1,  15.0015,           1
   3,  1547.15,  0.00334599,   0.000239246,  361.621,       1,  15.0015,           1
   4,  1842.32,  0.00398435,   0.00028489,   362.698,       1,  15.0015,           1

Run a query on ``ms1_df`` and ``ms2_df``. If you don't pass the data frames ``massql_engine.process_query``
will read data from the given file name.

.. code-block:: python
    :linenos:

    # Executing Query
    results_df = msql_engine.process_query(
        "QUERY scaninfo(MS1DATA) WHERE RTMIN=16",
        "small.mzML",
        ms1_df=ms1_df,
        ms2_df=ms2_df,
    )

    results_df.head()


.. csv-table:: results_df.head()
   :widths: 2 20 20 20 20 20
   :header: ,    scan,       rt,    mslevel,            i,    i_norm

   0,     139,  16.001,           1,  6.77786e+06,         1
   1,     140,  16.0095,          1,  9.65984e+06,         1
   2,     141,  16.0185,          1,  7.0933e+06,          1
   3,     143,  16.0268,          1,  7.51255e+06,         1
   4,     144,  16.0354,          1,  1.01007e+07,         1

In the resulting data frame each row represents a scan with the peak intensities summed up.
