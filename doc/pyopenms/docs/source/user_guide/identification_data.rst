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
  protein_hit.setMetaValue("target_decoy", "target")  # its a target protein

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
  ev.setAABefore("R")
  ev.setAAAfter("P")
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


Owning identification datasets (experimental)
*********************************************

``IdentificationData`` groups owned candidates into analysis runs with one shared
score contract. An analysis run can refer to several physical MS files. Peptides,
oligonucleotides and compounds have a string representation and explicit encoding.
The existing peptide/protein classes above remain available for legacy workflows.

Records are plain values whose constructors take the field names as keywords::

    ID = oms.IdentificationData
    data = ID()
    run = data.addRun("comet_1")  # MoleculeKind.PEPTIDE
    score = run.addScore(ID.ScoreDefinition(name="expect", higher_better=False, software="Comet"))
    run.setPrimaryScore(score)
    source = run.addSource(ID.SourceFile(path="BSA1.mzML"))
    query = run.addIdentification(source, ID.Observation(data_id="scan=1234", rt=1234.5, mz=582.32))
    albumin = ID.QualifiedAccession(database="uniprot.fasta", accession="P02769")
    evidence = ID.ParentEvidence(parent=albumin, start=65, end=74, before="K", after="T")
    run.addMatch(query, ID.MatchData(representation="LVNELTEFAK", charge=2, parent_evidence=[evidence]), [0.003])

Oligonucleotide and compound runs (``ID.MoleculeKind.OLIGONUCLEOTIDE``, ``ID.MoleculeKind.COMPOUND``)
use other encodings. Every run of a dataset declares the same score definitions::

    oligo = ID.MatchData(encoding=ID.Encoding.NA_SEQUENCE, representation="AUCGAUCG", charge=-3)
    compound = ID.MatchData(
        encoding=ID.Encoding.SMILES, representation="CC(=O)OC1=CC=CC=C1C(=O)O", formula="C9H8O4",
        charge=1, adduct=oms.AdductInfo.parseAdductString("M+H;1+"),
        identifiers=[ID.QualifiedAccession(database="HMDB", accession="HMDB0001879")])

Records with metadata also take ``metadata={name: value}``, and so does the ``parameters``
field of a score definition. Records compare by value (``==``) and print their fields. Mutable
records are unhashable; IDs (``QueryId``, ``MatchId``, ``ScoreId``, ``SourceId``), references
(``QueryReference``, ``MatchReference``), ``MoleculeIdentity`` and ``QualifiedAccession`` hash by
value and serve as dict keys and set members. A ``Match`` read from a run also compares its ID and
scores, not just its payload::

    observation = ID.Observation(data_id="scan=1", rt=12.5, metadata={"FWHM": 3.5})
    observation          # IdentificationData.Observation(data_id='scan=1', rt=12.5, mz=None, metadata={'FWHM': 3.5})
    accessions = {albumin, ID.QualifiedAccession(database="uniprot.fasta", accession="P02769")}  # one entry

Method names follow the C++ API (``getRun``, ``addMatch``, ``retainBest``). Getters such as
``getRun`` and ``getRuns`` return independent copies (snapshots). To edit a run of a dataset,
use ``run_view``: its methods act on the run inside the dataset, so there is nothing to write
back. ``addRun`` returns such a view as well::

    data = oms.FileHandler().loadIdentificationData("search.idXML")
    run = data.run_view(data.getRuns()[0].getIdentifier())
    run.retainBest(run.getPrimaryScore())
    oms.IdentificationDataFile.store("reduced.idparquet", data)

A view looks the run up by its UUID on every call. It keeps the dataset alive and raises
``KeyError`` if the run is no longer part of it. New query and match IDs are only ever
allocated by the run inside the dataset, which keeps feature links (run UUID plus match ID)
unambiguous. The identification data of a feature or consensus map is edited the same way
through ``identification_data_view()``::

    features.identification_data_view().run_view("search").setScore(match_id, score_id, 0.01)

The native output is a fresh directory containing a manifest and typed Parquet
tables. Existing destinations are rejected. Filtering does not renumber retained
records or automatically rerun protein inference. Choose preservation or removal
of inference explicitly when filtering a dataset or exporting a streaming subset::

    oms.IdentificationDataFile.filter(
        "search.idparquet", "subset.idparquet",
        lambda run_uuid, match: match.scores[0] is not None and match.scores[0] < 0.01,
        oms.IdentificationData.InferencePolicy.DISCARD)

This example assumes column zero is a smaller-is-better probability score in every
selected run. Inspect each run's definitions before applying a cross-run threshold.
``IdentificationDataFile.scan`` supports selecting runs and score columns while
skipping molecular payloads, metadata, evidence and annotations. Streaming avoids
loading every match into Python. Full loading and inference still require memory
proportional to their working data.

The tables of a bundle are plain Parquet files, so pyarrow, Polars or DuckDB can read
them directly. Each row carries its ``run_uuid``, and score columns are named after
their definitions (``score_pep``, ``score_q_value``); ``inspect`` lists them per run::

    print(oms.IdentificationDataFile.inspect("search.idparquet")[0].score_columns)
    # DuckDB: SELECT run_uuid, match_id, representation, score_q_value
    #         FROM 'search.idparquet/matches.parquet' WHERE score_q_value < 0.01

``IdentificationDataAdapter`` provides explicit legacy, feature and consensus
conversion. Strict export rejects information the target cannot express. The
permissive policy reports losses. Removing a match does not remove its measured
feature; live associations must be pruned or rejected.
