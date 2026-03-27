Feature Map Annotation with Peptide Identifications
===================================================

A :term:`feature map` is usually obtained from running a :term:`feature finder`, e.g. :py:class:`~.FeatureFinderAlgorithmPicked` (see `Feature Detection <feature_detection.html>`_).
A logical next step is to compare features across runs using `Map Alignment <map_alignment.html>`_ and `Feature Linking <feature_linking.html>`_.
However, most map aligners in pyOpenMS require features which are annotated with PSMs (see `Identication Data <identification_data.html>`_).
To link features to their respective PSMs (as obtained from a search engine, such as Comet), we can use the :py:class:`~.IDMapper`.


Step 0: Download Example data
-----------------------------

.. code-block:: python
    :linenos:
    
    import pyopenms as oms
    from urllib.request import urlretrieve

    base_url = (
        "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master/src/data/"
    )

    feature_file = "BSA1_F1.featureXML"
    urlretrieve(base_url + feature_file, feature_file)
    
    idxml_file = "BSA1_F1.idXML"
    urlretrieve(base_url + idxml_file, idxml_file)
        
Step 1: Load the Feature Map
----------------------------

First, load the FeatureMap from a `.featureXML` file:

.. code-block:: python

   import pyopenms as oms

   feature_map = oms.FeatureMap()
   oms.FeatureXMLFile().load(feature_file, feature_map)

Step 2: Load Peptide Identifications
------------------------------------

Next, load the PeptideIdentifications from an `.idXML` file:

.. code-block:: python

   peptide_ids = oms.PeptideIdentificationList()
   protein_ids = []
   oms.IdXMLFile().load(idxml_file, protein_ids, peptide_ids)

Step 3: Initialize and Configure `IDMapper`
-------------------------------------------

Now, configure `IDMapper` to apply **retention time (RT) and m/z tolerance settings**:

.. code-block:: python

   id_mapper = oms.IDMapper()
   params = id_mapper.getParameters()
   params.setValue("rt_tolerance", 5.0)  # RT tolerance in seconds
   params.setValue("mz_tolerance", 10.0)  # m/z tolerance in ppm
   id_mapper.setParameters(params)

Step 4: Annotate the FeatureMap
-------------------------------

Use the configured `IDMapper` to link peptide IDs to the FeatureMap:

.. code-block:: python

   # annotate() can optionally use the underlying raw MS data (spectra) to annotate unidentified MS/MS scans to features in the FeatureMap
   # We don't need this here, so we provide an empty default.
   spectra = oms.MSExperiment()
   id_mapper.annotate(feature_map, peptide_ids, protein_ids, True, True, spectra)

Step 5: Save the Annotated FeatureMap
--------------------------------------

Finally, store the modified FeatureMap back to a file:

.. code-block:: python

   oms.FeatureXMLFile().store("BSA1_F1_annotated.featureXML", feature_map)

.. tip::
   You can visualize the annotated FeatureMap using OpenMS visualization tools like `TOPPView`.


You have successfully **annotated a FeatureMap** with PeptideIdentifications using `IDMapper`. This allows further downstream analysis in (py)OpenMS workflows.

