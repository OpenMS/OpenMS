Glycans
=======

Glycans are chains or trees of monosaccharides. On proteins, N-glycans are attached to
asparagine and O-glycans mostly to serine or threonine. pyOpenMS 3.6 added classes that
describe a glycan as a composition or as a tree, and that compute theoretical fragment
ions of free glycans and of peptides with one glycan (glycopeptides).

Glycan Compositions
*******************

A :py:class:`~.GlycanComposition` stores ``components``, a list of ``(symbol, count)``
pairs. The symbols are those of ProForma glycan notation, such as ``HexNAc``
(N-acetylhexosamine), ``Hex`` (hexose), ``Fuc`` (fucose) and ``Neu5Ac`` or ``NeuAc``
(N-acetylneuraminic acid). The examples use G0F, a core-fucosylated N-glycan:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    glycan = oms.GlycanComposition()
    glycan.components = [("HexNAc", 4), ("Hex", 3), ("Fuc", 1)]

    residue_formulas = {"HexNAc": "C8H13NO5", "Hex": "C6H10O5", "Fuc": "C6H10O4"}
    residue_mass = sum(
        count * oms.EmpiricalFormula(residue_formulas[symbol]).getMonoWeight()
        for symbol, count in glycan.components
    )
    print(residue_mass, residue_mass + oms.EmpiricalFormula("H2O").getMonoWeight())

    peptidoform = oms.ProForma.parse("EEQYN[Glycan:HexNAc4Hex3Fuc1]STYR")
    print(oms.ProForma.toString(peptidoform, oms.ProForma.WriteMode.LOSSLESS))
    print(oms.ProForma.canCalculateMass(peptidoform))

.. code-block:: output

    1444.5338839348 1462.5444489986
    EEQYN[Glycan:HexNAc4Hex3Fuc1]STYR
    False

:py:class:`~.GlycanComposition` has no method for its mass. In a glycan, each
monosaccharide is a residue, the free monosaccharide minus one water. The first output
line shows the sum of the residue masses, which a glycan adds to a peptide, and the mass
of the free glycan, which has one more water. :py:class:`~.ModificationsDB` also contains
Unimod glycan compositions, such as ``dHex(1)Hex(3)HexNAc(4)`` for G0F (``dHex``:
deoxyhexose), which :py:class:`~.AASequence` accepts like other modifications.
:py:class:`~.ProForma` parses a glycan composition in ProForma notation and writes it
back, but cannot calculate the mass of a peptidoform that contains one.

Glycan Structures
*****************

A :py:class:`~.GlycanStructure` is a rooted tree.
``add_monosaccharide(symbol, parent, linkage)`` adds a node and returns its index. The
first node is the root, the monosaccharide at the reducing end; in a glycopeptide, it is
attached to the peptide. Every other node needs an existing parent, so a tree has no
cycles. The linkage is an optional label that fragment masses do not depend on.
``get_nodes()`` returns the nodes in the order they were added, ``get_composition()`` a
:py:class:`~.GlycanComposition` with one component per node:

.. code-block:: python
    :linenos:

    tree = oms.GlycanStructure()
    root = tree.add_monosaccharide("HexNAc")
    tree.add_monosaccharide("Fuc", root, "alpha1-6")
    core = tree.add_monosaccharide("HexNAc", root, "beta1-4")
    mannose = tree.add_monosaccharide("Hex", core, "beta1-4")
    for linkage in ("alpha1-3", "alpha1-6"):
        branch = tree.add_monosaccharide("Hex", mannose, linkage)
        tree.add_monosaccharide("HexNAc", branch, "beta1-2")

    for index, node in enumerate(tree.get_nodes()):
        print(index, node.monosaccharide, node.parent, repr(node.linkage))
    print(tree.get_composition().components)

.. code-block:: output

    0 HexNAc None ''
    1 Fuc 0 'alpha1-6'
    2 HexNAc 0 'beta1-4'
    3 Hex 2 'beta1-4'
    4 Hex 3 'alpha1-3'
    5 HexNAc 4 'beta1-2'
    6 Hex 3 'alpha1-6'
    7 HexNAc 6 'beta1-2'
    [('HexNAc', 1), ('Fuc', 1), ('HexNAc', 1), ('Hex', 1), ('Hex', 1), ('HexNAc', 1), ('Hex', 1), ('HexNAc', 1)]

Glycan Fragment Ions
********************

:py:class:`~.TheoreticalGlycanSpectrumGenerator` computes fragment ions of positively
charged glycans. ``get_fragments()`` takes a composition or a tree and returns a list of
``Fragment`` objects sorted by m/z, each with an ``ion_type``, a ``name``, the retained
``composition``, a ``neutral_mass``, a ``charge`` and ``get_mz()``. The settings are
attributes of ``TheoreticalGlycanSpectrumGenerator.Options``. Switching off the B and Y
ions (``add_b_ions`` and ``add_y_ions``) leaves the diagnostic (oxonium) ions:

.. code-block:: python
    :linenos:

    Generator = oms.TheoreticalGlycanSpectrumGenerator
    options = Generator.Options()
    options.add_b_ions = False
    options.add_y_ions = False
    for fragment in Generator(options).get_fragments(glycan):
        print(f"{fragment.get_mz():.4f}", fragment.name)
    print(len(Generator().get_fragments(glycan)))

.. code-block:: output

    126.0550 glycan:diagnostic:HexNAc-C2H6O3
    127.0390 glycan:diagnostic:Hex-H4O2
    138.0550 glycan:diagnostic:HexNAc-CH6O3
    144.0655 glycan:diagnostic:HexNAc-C2H4O2
    145.0495 glycan:diagnostic:Hex-H2O
    147.0652 glycan:diagnostic:Fuc
    163.0601 glycan:diagnostic:Hex
    168.0655 glycan:diagnostic:HexNAc-H4O2
    186.0761 glycan:diagnostic:HexNAc-H2O
    204.0867 glycan:diagnostic:HexNAc
    366.1395 glycan:diagnostic:Hex1HexNAc1
    73

A diagnostic ion name gives its residues and, if any, a neutral loss. With the default
options, the composition gives 73 fragments (last line): the diagnostic ions, Y0, and B
and Y ions for each sub-composition of 1 to 3 residues (``min_composition_size``,
``max_composition_size``) at charge 1 and 2 (``min_charge``, ``max_charge``). B ions
consist of glycan residues. Y ions contain the reducing end: water for a free glycan, the
peptide for a glycopeptide. Composition fragments are possible compositions and say
nothing about the structure. A tree gives only fragments that are connected parts of it:

.. code-block:: python
    :linenos:

    options = Generator.Options()
    options.add_diagnostic_ions = False
    options.max_charge = 1
    options.max_cleavages = 1
    for fragment in Generator(options).get_fragments(tree):
        if fragment.get_mz() < 600:
            print(f"{fragment.get_mz():.4f}", fragment.name)

.. code-block:: output

    19.0178 glycan:Y:tree:cuts=0;retained=0
    147.0652 glycan:B:tree:root=1;cuts=;retained=Fuc1
    204.0867 glycan:B:tree:root=5;cuts=;retained=HexNAc1
    204.0867 glycan:B:tree:root=7;cuts=;retained=HexNAc1
    366.1395 glycan:B:tree:root=4;cuts=;retained=Hex1HexNAc1
    366.1395 glycan:B:tree:root=6;cuts=;retained=Hex1HexNAc1
    368.1551 glycan:Y:tree:root=0;cuts=2;retained=Fuc1HexNAc1
    571.2345 glycan:Y:tree:root=0;cuts=3;retained=Fuc1HexNAc2

A structural name gives the root node of a B ion (``root``), the removed subtrees
(``cuts``, each named by its root node) and the retained composition. Fragments of equal
mass from different branches stay separate. Y0 contains only the reducing end, here
water. ``max_cleavages`` limits the broken bonds per fragment; with the default of 2, Y
ions can lose two subtrees and B ions one (internal fragments,
``add_internal_fragments``).

Glycopeptide Fragment Ions
**************************

``get_glycopeptide_fragments()`` takes the peptide without the glycan, the glycan
(composition or tree), the zero-based index of the glycosylated residue and the
fragmentation method. It returns glycan ions (diagnostic, B and Y ions) and backbone ions
(``IonType.PEPTIDE``). The name of a backbone ion that contains the site gives the glycan
it keeps and the 1-based site:

.. code-block:: python
    :linenos:

    peptide = oms.AASequence.fromString("EEQYNSTYR")
    Method = Generator.FragmentationMethod
    generator = Generator()
    for method in (Method.HCD, Method.ETD, Method.ETHCD):
        fragments = generator.get_glycopeptide_fragments(peptide, glycan, 4, method)
        backbone = [f for f in fragments if f.ion_type == Generator.IonType.PEPTIDE]
        n = len(backbone)
        print(method.name, n, "backbone ions,", len(fragments) - n, "glycan ions")
        for f in backbone:
            if f.charge == 1 and f.name.split(";")[0] in ("peptide:y6", "peptide:z6"):
                print(f"  {f.get_mz():.4f} {f.name}")

.. code-block:: output

    HCD 48 backbone ions, 73 glycan ions
      803.3682 peptide:y6;glycan=0;site=N5
      1006.4476 peptide:y6;glycan=HexNAc1;site=N5
    ETD 32 backbone ions, 0 glycan ions
      2231.8834 peptide:z6;glycan=Fuc1Hex3HexNAc4;site=N5
    ETHCD 64 backbone ions, 73 glycan ions
      2231.8834 peptide:z6;glycan=Fuc1Hex3HexNAc4;site=N5
      2247.9021 peptide:y6;glycan=Fuc1Hex3HexNAc4;site=N5

With HCD, b and y ions keep no glycan (``glycan=0``) or one HexNAc, if the glycan has
one. With ETD, c and z ions keep the complete glycan, and there are no glycan ions. EThcD
gives b, y, c and z ions with the complete glycan, and the glycan ions. The z ions are
z+1 radical ions. Backbone ions without the site, such as y3 (TYR), appear once, without
glycan.

``Options.peptide_retention`` changes what backbone ions keep. It maps an ion series
(``"b"``, ``"y"``, ``"c"`` or ``"z"``) to a ``PeptideRetention`` that replaces the default
for this series: ``intact`` keeps the complete glycan (default ``False``), ``stripped``
no glycan (default ``True``) and ``stubs`` lists compositions that are part of the
glycan, for example HexNAc1Fuc1 for a core-fucosylated stub.

Annotated Spectra
*****************

``to_spectrum()`` converts fragments into an :term:`MS2` :py:class:`~.MSSpectrum` sorted
by m/z. Every peak has intensity 1; intensities are not predicted. The string data array
``IonNames`` holds one annotation per peak, the integer data array ``Charges`` its charge.
mzPAF has no glycan ion series, so each annotation is a named compound,
``_[name]^charge``, in the form that :py:class:`~.MzPAF` reads and writes. Here, the tree
is passed instead of the composition, so the glycan ions are connected parts of the tree:

.. code-block:: python
    :linenos:

    fragments = generator.get_glycopeptide_fragments(peptide, tree, 4, Method.HCD)
    spectrum = Generator.to_spectrum(fragments)
    annotations = spectrum.getStringDataArrays()[0]
    charges = spectrum.getIntegerDataArrays()[0]
    print(spectrum.getMSLevel(), spectrum.size(), annotations.getName(), charges.getName())
    for peak, annotation in zip(spectrum, annotations):
        if 1185 < peak.getMZ() < 1220:
            print(f"{peak.getMZ():.4f}", peak.getIntensity(), annotation)

.. code-block:: output

    2 147 IonNames Charges
    1189.5120 1.0 _[glycan:Y:tree:cuts=0;retained=0;site=N5]^1
    1215.9869 1.0 _[glycan:Y:tree:root=0;cuts=7;retained=Fuc1Hex3HexNAc3;site=N5]^2
    1215.9869 1.0 _[glycan:Y:tree:root=0;cuts=5;retained=Fuc1Hex3HexNAc3;site=N5]^2
    1218.4797 1.0 _[peptide:b8;glycan=HexNAc1;site=N5]^1

Y ions of a glycopeptide contain the peptide, and Y0 is the peptide without glycan. In
glycopeptide spectra, the names of glycan ions also end with the site.
