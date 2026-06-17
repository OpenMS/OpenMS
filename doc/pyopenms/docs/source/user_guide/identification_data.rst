Identification Data
====================

In OpenMS, identifications of peptides, proteins and small molecules are stored
in dedicated data structures. These data structures are typically stored to disc
as idXML or mzIdentML files. The highest-level structure is
:py:class:`~.ProteinIdentification`. It stores all identified proteins of an identification
run (usually all IDs from a single HPLC-MS run) as :py:class:`~.ProteinHit` objects plus additional metadata (search parameters, etc.). Each
:py:class:`~.ProteinHit` represents a potential protein which may be present in the sample. The `ProteinHit` contains the actual protein identifier (also known as `accession`), an associated score, and
the protein sequence. The latter may be omitted to reduce memory consumption.

A :py:class:`~.PeptideIdentification` object stores the
data corresponding to a single identified spectrum or feature. It has members
for the retention time, m/z, and a vector of :py:class:`~.PeptideHit` objects. Each :py:class:`~.PeptideHit`
stores the information of a specific :term:`peptide-spectrum match` (:term:`PSM`), e.g., the score
and the peptide sequence. Each :py:class:`~.PeptideHit` also contains a vector of
:py:class:`~.PeptideEvidence` objects which store the reference to one or more (in the case the
peptide maps to multiple proteins) proteins and the position therein.

.. NOTE::
  Proteins and their corresponding peptides are linked by a common identifier (e.g., a unique string of time and date of the search).
  The Identifier can be set using the :py:meth:`~.ProteinIdentification.setIdentifier` method in
  :py:class:`~.ProteinIdentification` and :py:class:`~.PeptideIdentification`.
  Similarly :py:meth:`~.ProteinIdentification.getIdentifier` can be used to check the identifier.
  Using this link one can retrieve search meta data (which is stored at the protein level) for individual peptides.

  .. code-block:: python
    :linenos:
    
    import pyopenms as oms

    protein_id = oms.ProteinIdentification()
    peptide_id = oms.PeptideIdentification()

    # Sets the Identifier
    protein_id.setIdentifier("IdentificationRun1")
    peptide_id.setIdentifier("IdentificationRun1")

    # Prints the Identifier
    print("Protein Identifier -", protein_id.getIdentifier())
    print("Peptide Identifier -", peptide_id.getIdentifier())

  .. code-block:: output

    Protein Identifier - IdentificationRun1
    Peptide Identifier - IdentificationRun1

Protein Identification
***********************

We can create an object of type :py:class:`~.ProteinIdentification`  and populate it with
:py:class:`~.ProteinHit` objects as follows:

.. see doc/code_examples/Tutorial_IdentificationClasses.cpp

.. code-block:: python
  :linenos:

  import pyopenms as oms

  # Create new protein identification object corresponding to a single search
  protein_id = oms.ProteinIdentification()
  protein_id.setIdentifier("IdentificationRun1")

  # Each ProteinIdentification object stores a vector of protein hits
  protein_hit = oms.ProteinHit()
  protein_hit.setAccession("sp|MyAccession")
  protein_hit.setSequence("PEPTIDERDLQMTQSPSSLSVSVGDRPEPTIDE")
  protein_hit.setScore(1.0)
  protein_hit.setMetaValue("target_decoy", b"target")  # its a target protein

  protein_id.setHits([protein_hit])

We have now added a single :py:class:`~.ProteinHit` with the accession ``sp|MyAccession`` to
the :py:class:`~.ProteinIdentification` object (note how on line 14 we directly added a
list of size 1).  We can continue to add meta-data for the whole identification
run (such as search parameters):

.. code-block:: python
  :linenos:

  now = oms.DateTime.now()
  date_string = now.getDate()
  protein_id.setDateTime(now)

  # Example of possible search parameters
  search_parameters = (
      oms.SearchParameters()
  )  # ProteinIdentification::SearchParameters
  search_parameters.db = "database"
  search_parameters.charges = "+2"
  protein_id.setSearchParameters(search_parameters)

  # Some search engine meta data
  protein_id.setSearchEngineVersion("v1.0.0")
  protein_id.setSearchEngine("SearchEngine")
  protein_id.setScoreType("HyperScore")

  # Iterate over all protein hits
  for hit in protein_id.getHits():
      print("Protein hit accession:", hit.getAccession())
      print("Protein hit sequence:", hit.getSequence())
      print("Protein hit score:", hit.getScore())


PeptideIdentification
**********************

Next, we can also create a :py:class:`~.PeptideIdentification` object and add
corresponding :py:class:`~.PeptideHit` objects:

.. code-block:: python
  :linenos:

  peptide_id = oms.PeptideIdentification()

  peptide_id.setRT(1243.56)
  peptide_id.setMZ(440.0)
  peptide_id.setScoreType("ScoreType")
  peptide_id.setHigherScoreBetter(False)
  peptide_id.setIdentifier("IdentificationRun1")

  # define additional meta value for the peptide identification
  peptide_id.setMetaValue("AdditionalMetaValue", "Value")

  # create a new PeptideHit (best PSM, best score)
  peptide_hit = oms.PeptideHit()
  peptide_hit.setScore(1.0)
  peptide_hit.setRank(1)
  peptide_hit.setCharge(2)
  peptide_hit.setSequence(oms.AASequence.fromString("DLQM(Oxidation)TQSPSSLSVSVGDR"))

  ev = oms.PeptideEvidence()
  ev.setProteinAccession("sp|MyAccession")
  ev.setAABefore(b"R")
  ev.setAAAfter(b"P")
  ev.setStart(123)  # start and end position in the protein
  ev.setEnd(141)
  peptide_hit.setPeptideEvidences([ev])

  # create a new PeptideHit (second best PSM, lower score)
  peptide_hit2 = oms.PeptideHit()
  peptide_hit2.setScore(0.5)
  peptide_hit2.setRank(2)
  peptide_hit2.setCharge(2)
  peptide_hit2.setSequence(oms.AASequence.fromString("QDLMTQSPSSLSVSVGDR"))
  peptide_hit2.setPeptideEvidences([ev])

  # add PeptideHit to PeptideIdentification
  peptide_id.setHits([peptide_hit, peptide_hit2])


This allows us to represent single spectra (:py:class:`~.PeptideIdentification` at m/z
:math:`440.0` and *rt* :math:`1234.56`) with possible identifications that are ranked by score.
In this case, apparently two possible peptides match the spectrum which have
the first three amino acids in a different order "DLQ" vs "QDL").

We can now display the peptides we just stored:

.. code-block:: python
  :linenos:

  # Iterate over PeptideIdentification
  peptide_ids = oms.PeptideIdentificationList()
  peptide_ids.push_back(peptide_id)
  for peptide_id in peptide_ids:
      # Peptide identification values
      print("Peptide ID m/z:", peptide_id.getMZ())
      print("Peptide ID rt:", peptide_id.getRT())
      print("Peptide ID score type:", peptide_id.getScoreType())
      # PeptideHits
      for hit in peptide_id.getHits():
          print(" - Peptide hit rank:", hit.getRank())
          print(" - Peptide hit sequence:", hit.getSequence())
          print(" - Peptide hit score:", hit.getScore())
          print(
              " - Mapping to proteins:",
              [ev.getProteinAccession() for ev in hit.getPeptideEvidences()],
          )



Storage on Disk
***************

Finally, we can store the peptide and protein identification data in a
:py:class:`~.IdXMLFile` (a OpenMS internal file format which we have previously
discussed :ref:`anchor-other-id-data`) which we would do as follows:

.. code-block:: python
  :linenos:

  # Store the identification data in an idXML file
  oms.IdXMLFile().store("out.idXML", [protein_id], peptide_ids)
  # and load it back into memory
  prot_ids = []
  pep_ids = oms.PeptideIdentificationList()
  oms.IdXMLFile().load("out.idXML", prot_ids, pep_ids)

  # Iterate over all protein hits
  for protein_id in prot_ids:
      for hit in protein_id.getHits():
          print("Protein hit accession:", hit.getAccession())
          print("Protein hit sequence:", hit.getSequence())
          print("Protein hit score:", hit.getScore())
          print("Protein hit target/decoy:", hit.getMetaValue("target_decoy"))

  # Iterate over PeptideIdentification
  for peptide_id in pep_ids:
      # Peptide identification values
      print("Peptide ID m/z:", peptide_id.getMZ())
      print("Peptide ID rt:", peptide_id.getRT())
      print("Peptide ID score type:", peptide_id.getScoreType())
      # PeptideHits
      for hit in peptide_id.getHits():
          print(" - Peptide hit rank:", hit.getRank())
          print(" - Peptide hit sequence:", hit.getSequence())
          print(" - Peptide hit score:", hit.getScore())
          print(
              " - Mapping to proteins:",
              [ev.getProteinAccession() for ev in hit.getPeptideEvidences()],
          )

You can inspect the ``out.idXML`` XML file produced here, and you will find a :py:class:`~.ProteinHit` entry for
the protein that we stored and two :py:class:`~.PeptideHit` entries for the two peptides stored on disk.
