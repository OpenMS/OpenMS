Standard Notations: ProForma, USI and mzPAF
===========================================

The Proteomics Standards Initiative (PSI) of the Human Proteome Organization (HUPO)
defines text notations for peptides, spectra and peak annotations. pyOpenMS reads and
writes three of them:

* ProForma describes a peptide with its modifications, a peptidoform, in one string.
* The Universal Spectrum Identifier (USI) names one spectrum in a public data set.
* mzPAF, the Peak Annotation Format, names the ion that explains a peak.

ProForma
********

In `ProForma <https://github.com/HUPO-PSI/ProForma>`__, a modification follows its residue
in square brackets. It can be a name (``[Oxidation]``), an accession (``[UNIMOD:35]``), a
mass shift (``[+15.9949]``) or a formula (``[Formula:O]``). Terminal modifications stand
before or after the sequence, joined by a hyphen: ``[Acetyl]-PEPTIDE-[Amidated]``.

pyOpenMS parses ProForma into a :py:class:`~.Peptidoform`, one peptide chain, or a
:py:class:`~.PeptidoformIon`, one or more chains with an optional charge. Their methods
call static methods of :py:class:`~.ProForma`, for example ``ProForma.parse()``.

Reading and Writing
~~~~~~~~~~~~~~~~~~~

.. code-block:: python
    :linenos:

    import pyopenms as oms

    for text in ["EM[Oxidation]K", "EM[UNIMOD:35]K", "[Acetyl]-EM[+15.99]K"]:
        pf = oms.Peptidoform.fromString(text)
        print(pf.toString(), pf.toString(oms.ProForma.WriteMode.CANONICAL))

.. code-block:: output

    EM[Oxidation]K EM[Oxidation]K
    EM[UNIMOD:35]K EM[UNIMOD:35]K
    [Acetyl]-EM[+15.99]K [Acetyl]-EM[+15.9900]K

``Peptidoform.fromString()`` parses a string and returns a :py:class:`~.Peptidoform`.
``toString()`` writes it as ProForma. The default mode, ``ProForma.WriteMode.LOSSLESS``,
keeps the text of each mass shift. ``ProForma.WriteMode.CANONICAL`` writes each mass
shift with a sign and four decimal places, and each localization score with two. The
modes differ only in these numbers. ``print()`` shows the LOSSLESS form.

Masses and Charged Peptidoforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python
    :linenos:

    pf = oms.Peptidoform.fromString("EM[Oxidation]K")
    print(pf.getMonoWeight(), pf.getMZ(2))

    try:
        oms.Peptidoform.fromString("EM[Oxidation]K/2")
    except RuntimeError as error:
        print(error)

    ion = oms.PeptidoformIon.fromString("EM[Oxidation]K/2")
    print(ion, len(ion.chains))
    print(ion.getMonoWeight(), ion.getMZ())

.. code-block:: output

    422.183522687 212.099037810271
    Unexpected characters after peptidoform in: EM[Oxidation]K/2
    EM[Oxidation]K/2 1
    422.183522687 212.099037810271

``getMonoWeight()`` returns the monoisotopic mass of the neutral peptidoform in Da.
``getMZ(2)`` returns the m/z of the :chem:`[M+2H]2+` ion. In ProForma, the charge follows
a slash. ``Peptidoform.fromString()`` does not accept it: like any string it cannot parse,
it raises a ``RuntimeError``. ``PeptidoformIon.fromString()`` reads the charge. The
``chains`` attribute of the ion lists one :py:class:`~.Peptidoform` per chain;
cross-linked chains are separated by ``//``. ``PeptidoformIon.getMZ()`` takes no argument
and uses the charge of the string.

Conversion to and from AASequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Other parts of pyOpenMS, for example :py:class:`~.PeptideHit`, store a peptide as an
:py:class:`~.AASequence` (see :doc:`peptides_proteins`):

.. code-block:: python
    :linenos:

    seq = pf.toAASequence()
    print(seq, seq.getMonoWeight())
    seq = oms.AASequence.fromString(".(Acetyl)PEPS(Phospho)TIDEK")
    print(oms.Peptidoform.fromAASequence(seq))

    pf = oms.Peptidoform.fromString("[Phospho]?PEPTIDE")
    print(pf.isRepresentableAsAASequence())
    print(pf.getAASequenceConversionIssues()[0].description)
    print(pf.getMonoWeight(), pf.toAASequence().getMonoWeight())

.. code-block:: output

    EM(Oxidation)K 422.183522687
    [UNIMOD:1]-PEPS[UNIMOD:21]TIDEK
    False
    Peptidoform contains unlocalised modifications
    879.3263006906999 799.3599696906999

``EM[Oxidation]K`` becomes the :py:class:`~.AASequence` ``EM(Oxidation)K``, with the same
mass. ``Peptidoform.fromAASequence()`` uses the UniMod accession of each modification that
has one. An :py:class:`~.AASequence` cannot hold everything that ProForma can express:
``[Phospho]?PEPTIDE`` states a phosphorylation at an unknown site, and the mass of the
peptidoform includes it. ``isRepresentableAsAASequence()`` tells whether a conversion
keeps everything, and ``getAASequenceConversionIssues()`` lists what it would lose.
Without an argument, ``toAASequence()`` uses ``ProForma.ConversionPolicy.BEST_EFFORT`` and
leaves out what does not fit, here the phosphorylation. With
``ProForma.ConversionPolicy.FAIL_ON_LOSS``, it raises a ``RuntimeError`` instead.

JSON
~~~~

:py:class:`~.Peptidoform` has no attribute for its residues. ``ProForma.peptidoformToJSON()``
returns the parsed structure as a JSON string, and ``Peptidoform.fromJSON()`` reads it:

.. code-block:: python
    :linenos:

    import json

    text = oms.ProForma.peptidoformToJSON(oms.Peptidoform.fromString("EM[Oxidation]K"))
    print(json.loads(text)["sequence"][1]["value"])
    print(oms.Peptidoform.fromJSON(text))

.. code-block:: output

    {'amino_acid': 'M', 'modifications': [[{'tag': {'type': 'named_mod', 'value': {'name': 'Oxidation'}}}]]}
    EM[Oxidation]K

Universal Spectrum Identifier (USI)
***********************************

A `USI <https://www.psidev.info/usi>`__ identifies one spectrum. It has the form
``mzspec:<collection>:<run>:<index type>:<index>``. The collection is a public data set,
such as the ProteomeXchange data set ``PXD000561``, or a spectral library. The index type
is ``scan``, ``index`` or ``nativeId``. An optional last part, the interpretation, is the
ProForma peptidoform ion that explains the spectrum. With it, the USI describes a
:term:`PSM`.

.. code-block:: python
    :linenos:

    text = "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2"
    print(oms.USI.isValidUSI(text), oms.USI.isValidUSI("mzspec:PXD000561:run:17555"))

    usi = oms.USI(text)
    print(usi.getCollection(), usi.getMSRun())
    print(usi.getIndexType(), usi.getIndex())
    print(usi.getInterpretation())
    print(oms.PeptidoformIon.fromString(usi.getInterpretation()).getMZ())

.. code-block:: output

    True False
    PXD000561 Adult_Frontalcortex_bRP_Elite_85_f09
    IndexType.SCAN 17555
    VLHPLEGAVVIIFK/2
    767.971422428621

``USI.isValidUSI()`` checks the structure of a string: the ``mzspec:`` prefix, a
collection, a run, a known index type and an index. The second string lacks the index
type. The check does not parse the interpretation. The :py:class:`~.USI` constructor
parses a string and raises a ``RuntimeError`` if it is not valid. The getters return the
parts, and ``PeptidoformIon.fromString()`` parses the interpretation.

Building a USI for an Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :py:class:`~.PeptideIdentification` holds the peptide hits of one spectrum.
``buildUSI()`` creates the USI of that spectrum:

.. code-block:: python
    :linenos:

    pep_id = oms.PeptideIdentification()
    pep_id.setSpectrumReference("controllerType=0 controllerNumber=1 scan=17555")
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("VLHPLEGAVVIIFK"))
    hit.setCharge(2)
    pep_id.insertHit(hit)
    pep_id.sort()

    run = "Adult_Frontalcortex_bRP_Elite_85_f09"
    print(pep_id.buildUSI(run, "PXD000561", True).toString())
    print(pep_id.buildUSI(run).toString())

.. code-block:: output

    mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2
    mzspec:local:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555

``buildUSI()`` takes the run name, the collection and whether to add the interpretation.
It writes the run name unchanged, so do not pass a path. It reads the scan number from
the spectrum reference, the native ID of the spectrum. If it finds none, it writes the
whole native ID with the index type ``nativeId``. The interpretation is the sequence and
charge of the first :py:class:`~.PeptideHit`; ``sort()`` puts the best hit first. Without
a collection, ``buildUSI()`` writes ``local``, which OpenMS uses for data that is not in a
public repository. Without a spectrum reference, it returns an invalid USI, whose
``toString()`` is empty.

In place of the run name, ``buildUSI()`` also accepts an
:py:class:`~.IdentifierMSRunMapper`, created from the :py:class:`~.ProteinIdentification`
objects (``oms.IdentifierMSRunMapper(protein_ids)``). It then takes the run from the protein identification that has the same
identifier as the peptide identification, and drops the directories from its path.

mzPAF Peak Annotations
**********************

`mzPAF <https://github.com/HUPO-PSI/mzPAF>`__ describes the ion behind a peak of a fragment
ion spectrum (:term:`MS2`). For example, ``y4-H2O^2/1.2ppm`` is the y4 ion after the loss
of water, with charge 2 and a mass error of 1.2 ppm. Alternative annotations of one peak
are separated by commas.

.. code-block:: python
    :linenos:

    ann = oms.MzPAF.parse("y4-H2O^2/1.2ppm")
    print(ann.ion_series, ann.ordinal, ann.charge)
    print([str(loss.formula) for loss in ann.neutral_losses])
    print(ann.mass_delta.value, ann.mass_delta.unit)
    print(oms.MzPAF.toString(ann))

    anns = oms.MzPAF.parseMultiple("b2-H2O,y4^2")
    print(anns.size(), oms.MzPAF.toStringMultiple(anns))

.. code-block:: output

    MzPAFIonSeries.Y 4 2
    ['H2O1']
    1.2 MzPAFDeltaUnit.PPM
    y4-H2O1^2/1.2ppm
    2 b2-H2O1,y4^2

``MzPAF.parse()`` returns an :py:class:`~.MzPAFAnnotation`. Its attributes hold the parts
of the annotation. A part that the string does not contain is ``None``, for example
``charge`` for ``y4``; ``neutral_losses`` is then an empty list. Each neutral loss is an
:py:class:`~.MzPAFNeutralLoss` with an :py:class:`~.EmpiricalFormula`. ``mass_delta`` has a
value and a unit, ``MzPAFDeltaUnit.PPM`` or ``MzPAFDeltaUnit.DALTON``. ``MzPAF.toString()``
writes the count of every element, so ``-H2O`` becomes ``-H2O1``, the same formula.

``MzPAF.parseMultiple()`` returns an :py:class:`~.MzPAFPeakAnnotations`, whose
``annotations`` attribute lists the annotations; ``MzPAF.toStringMultiple()`` writes them
back. ``MzPAF.parse()`` returns only the first one. ``parse()`` and ``parseMultiple()``
raise a ``RuntimeError`` for a string they cannot parse; ``tryParse()`` and
``tryParseMultiple()`` return ``None`` instead.

A :py:class:`~.PeptideHit` stores the annotated peaks of its spectrum as a list of
:py:class:`~.PeptideHit_PeakAnnotation` objects (``getPeakAnnotations()``). Each has an
annotation string, a charge, an m/z and an intensity. :py:class:`~.MzPAF` converts in
both directions:

.. code-block:: python
    :linenos:

    peak = oms.MzPAF.toPeakAnnotation(oms.MzPAF.parse("y1"), 147.1128, 1000.0)
    print(peak.annotation, peak.charge, peak.mz, peak.intensity)
    print(oms.MzPAF.fromPeakAnnotation(peak).size())

    peak.annotation = "y1+"
    print(oms.MzPAF.fromPeakAnnotation(peak).size())

.. code-block:: output

    y1 1 147.1128 1000.0
    1
    0

``MzPAF.toPeakAnnotation()`` sets the annotation string to the mzPAF string and the charge
to the charge of the annotation, or to 1 if the annotation has none.
``MzPAF.fromPeakAnnotation()`` parses only the annotation string. It returns an
:py:class:`~.MzPAFPeakAnnotations`, which is empty if the string is not mzPAF, like ``y1+``.
