Quality Control 
===============

With :py:class:`~.MzQCFile` pyOpenMS provides a tool to export quality control data from an experiment to a mzQC file following the
`HUPO-PSI specifications
<https://github.com/HUPO-PSI/mzQC>`_. MzQC is a reporting and exchange format for mass spectrometry
quality control data in JSON format.

:py:meth:`~.MzQCFile.store()` requires, besides the ``input_file`` and ``output_file`` path,
at least a :py:class:`~.MSExperiment`. Optionally, a :py:class:`~.FeatureMap` and/or lists of
:py:class:`~.ProteinIdentification` and :py:class:`~.PeptideIdentification` can be provided to calculate additional
proteomics and metabolomics quality metrics.

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master"

    # minimal required informations are input_file, output_file and a MSExperiment
    input_file = "QCCalculator_input.mzML"
    output_file = "result.mzQC"

    # load mzML file in MSExperiment
    urlretrieve(gh + "/src/data/QCCalculator_input.mzML", "input.mzML")
    exp = oms.MSExperiment()
    oms.MzMLFile().load("input.mzML", exp)

    # prepare metadata for the mzQC file (optional but recommended)
    contact_name = "name of the person creating the mzQC file"
    contact_address = "contact address (mail/e-mail or phone) of the person creating the mzQC file"
    description = "description and comments about the mzQC file contents"
    label = "unique and informative label for the run"

    feature_map = oms.FeatureMap()
    # OPTIONAL: load featureXML file in FeatureMap()
    urlretrieve(
        gh + "/src/data/FeatureFinderMetaboIdent_1_output.featureXML",
        "features.featureXML",
    )
    oms.FeatureXMLFile().load("features.featureXML", feature_map)

    prot_ids = []  # list of ProteinIdentification()
    pep_ids = oms.PeptideIdentificationList()  # list of PeptideIdentification()
    # OPTIONAL: get protein and peptide identifications from idXML file
    urlretrieve(gh + "/src/data/OpenPepXL_output.idXML", "ids.idXML")
    oms.IdXMLFile().load("ids.idXML", prot_ids, pep_ids)

    # create MzQCFile object and use store() to calculate and write the mzQC file
    mzqc = oms.MzQCFile()
    mzqc.store(
        input_file,
        output_file,
        exp,  # minimal required information
        contact_name,
        contact_address,
        description,
        label,  # optional, but recommended, can be empty
        feature_map,
        prot_ids,
        pep_ids,
    )  # optional, can be empty
