from pyopenms import AASequence, DecoyGenerator


def test_debruijn_decoys_preserve_shared_repeats_and_seed():
    proteins = [AASequence.fromString("APEPTIDEG"), AASequence.fromString("KPEPTIDEA")]

    def generate():
        generator = DecoyGenerator()
        generator.startDeBruijn(2, 4711)
        for protein in proteins:
            generator.addProteinToDeBruijn(protein)
        generator.finalizeDeBruijn()
        return [generator.deBruijn(protein).toString() for protein in proteins]

    first = generate()
    assert first == generate()
    assert len(first[0]) == len(proteins[0].toString())
    assert len(first[1]) == len(proteins[1].toString())
    assert first[0][3:8] == first[1][3:8]
