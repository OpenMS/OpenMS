"""RNA decoy residue transformations retain terminal chemistry."""
import pytest
import pyopenms as oms


def test_rna_decoy_reversal_and_identity():
    sequence = oms.NASequence.fromString("pAUGCp")
    assert sequence.getReversed().toString() == "pCGUAp"
    assert sequence.getReversed(True, True).toString() == "pAGUCp"
    assert sequence.getSequenceIdentity(sequence) == 1
    assert oms.NASequence().getReversed().empty()
    with pytest.raises(ValueError):
        sequence.getSequenceIdentity(oms.NASequence.fromString("A"))


def test_rna_decoy_shuffle_is_seeded_and_bounded():
    sequence = oms.NASequence.fromString("pAUGCAUGCAUGCAUGCp")
    first = sequence.getShuffled(42, 0.5, 30, True, True)
    second = sequence.getShuffled(42, 0.5, 30, True, True)
    assert first.toString() == second.toString()
    assert first.size() == sequence.size()
    assert first.getSequenceIdentity(sequence) <= 0.5
    assert first.getFivePrimeMod().getCode() == sequence.getFivePrimeMod().getCode()
    assert first.getThreePrimeMod().getCode() == sequence.getThreePrimeMod().getCode()
    assert first[0].getCode() == sequence[0].getCode()
    assert first[first.size() - 1].getCode() == sequence[sequence.size() - 1].getCode()
    assert oms.NASequence.fromString("AU").getShuffled(42).toString() == "AU"
