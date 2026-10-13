"""Modification-aware RNA digestion coordinates and candidate regressions."""
import pyopenms as oms


def test_cleavage_sensitive_parent_terminal_candidate():
    digestion = oms.RNaseDigestion()
    digestion.setEnzyme("RNase_T1")
    digestion.setMissedCleavages(0)
    sequence = oms.NASequence.fromString("AG")
    assert digestion.getFragmentPositions(sequence) == [(0, 2)]
    mod = oms.RibonucleotideDB.getInstance().getRibonucleotide("m2,2G")
    groups = digestion.inferCleavageSensitiveMods({mod})
    assert len(groups.combined()) == 1
    products = digestion.digestWithCleavageSensitiveMods(sequence, groups, 1)
    assert {(p.fragment.toString(), tuple(p.position)) for p in products} == {
        ("AG", (0, 2)), ("A[m2,2G]", (0, 2))
    }



def test_coordinate_products_and_copy_constructors():
    digestion = oms.RNaseDigestion()
    digestion.setEnzyme("RNase_T1")
    sequence = oms.NASequence.fromString("AGUC")
    products = digestion.digestProducts(sequence)
    assert [(p.fragment.toString(), tuple(p.position)) for p in products] == [
        ("AGp", (0, 2)), ("UC", (2, 2))
    ]
    copy = oms.RNaseDigestionProduct(products[0])
    products[0].position = (8, 2)
    assert tuple(copy.position) == (0, 2)
    assert oms.RNaseDigestionProduct().fragment.empty()
    groups = digestion.inferCleavageSensitiveMods({
        oms.RibonucleotideDB.getInstance().getRibonucleotide("m2,2G")
    })
    copied_groups = oms.RNaseCleavageSensitiveModGroups(groups)
    del groups
    assert len(copied_groups.combined()) == 1
    assert not oms.RNaseCleavageSensitiveModGroups().combined()
