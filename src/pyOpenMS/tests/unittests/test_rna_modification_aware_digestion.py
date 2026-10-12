"""Modification-aware RNA digestion coordinates and candidate regressions."""
import gc

import pytest
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


@pytest.mark.parametrize("input_kind", ["set", "generator"])
def test_combined_borrows_database_entries_after_inputs_are_collected(input_kind):
    digestion = oms.RNaseDigestion()
    digestion.setEnzyme("RNase_T1")
    db = oms.RibonucleotideDB.getInstance()
    if input_kind == "set":
        modifications = {db.getRibonucleotide("m2,2G")}
        groups = digestion.inferCleavageSensitiveMods(modifications)
        modifications.clear()
    else:
        modifications = (db.getRibonucleotide(code) for code in ["m2,2G"])
        groups = digestion.inferCleavageSensitiveMods(modifications)
    del modifications
    gc.collect()

    # Each temporary combined set is collected before the next call. It must never own the DB entry.
    for _ in range(3):
        assert {mod.getCode() for mod in groups.combined()} == {"m2,2G"}
        gc.collect()
    assert db.getRibonucleotide("m2,2G").getCode() == "m2,2G"


def test_inferred_copies_and_products_outlive_mutated_python_modifications():
    digestion = oms.RNaseDigestion()
    digestion.setEnzyme("RNase_T1")
    input_sequence = oms.NASequence.fromString("A[m2,2G]")
    modification = input_sequence[1]  # NASequence getters return Python-owned copies.
    modifications = {modification}
    groups = digestion.inferCleavageSensitiveMods(modifications)
    copied_groups = oms.RNaseCleavageSensitiveModGroups(groups)
    products = digestion.digestWithCleavageSensitiveMods(
        oms.NASequence.fromString("AG"), copied_groups, 1
    )
    modified_product = next(p for p in products if p.fragment.toString() == "A[m2,2G]")
    copied_product = oms.RNaseDigestionProduct(modified_product)
    copied_fragment = oms.NASequence(copied_product.fragment)

    modifications.clear()
    modification.setCode("python_input_was_mutated")
    modification.setMonoMass(modification.getMonoMass() + 50.0)
    # Editing an accepted input copy after inference must not edit the stored group or the products.
    assert {mod.getCode() for mod in copied_groups.combined()} == {"m2,2G"}
    assert copied_product.fragment.toString() == "A[m2,2G]"
    del modification, modifications, input_sequence, groups, copied_groups, products, modified_product
    gc.collect()

    assert copied_product.fragment.toString() == "A[m2,2G]"
    assert tuple(copied_product.position) == (0, 2)
    del copied_product
    gc.collect()
    assert copied_fragment.toString() == "A[m2,2G]"
    assert copied_fragment.getMonoWeight() == oms.NASequence.fromString("A[m2,2G]").getMonoWeight()


def test_inference_rejects_unknown_and_edited_modification_definitions():
    digestion = oms.RNaseDigestion()
    digestion.setEnzyme("RNase_T1")
    db = oms.RibonucleotideDB.getInstance()
    unknown = oms.Ribonucleotide(db.getRibonucleotide("m2,2G"))
    unknown.setCode("unregistered_python_rnase_modification")
    with pytest.raises(ValueError, match="not owned by RibonucleotideDB"):
        digestion.inferCleavageSensitiveMods({unknown})

    edited = oms.Ribonucleotide(db.getRibonucleotide("m2,2G"))
    edited.setMonoMass(edited.getMonoMass() + 50.0)
    with pytest.raises(ValueError, match="differs from the RibonucleotideDB entry"):
        digestion.inferCleavageSensitiveMods({edited})
