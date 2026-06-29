"""
Deprecation-shim tests for the pandas DataFrame compatibility wrappers
(OpenMS issue #9460, Section 2: pyOpenMS nanobind rewrite).

``pyopenms.peptide_identifications_to_df()`` and
``pyopenms.update_scores_from_df()`` are deprecated module-level wrappers
(``pyopenms/_dataframes_compat.py``) kept for backward compatibility after the
DataFrame helpers were renamed onto ``PeptideIdentificationList``. Each must

  * still resolve and work, and
  * emit a ``DeprecationWarning`` naming the ``to_df()`` /
    ``update_scores_from_df()`` replacement,

then forward to the replacement. The general dataframe tests call these wrappers
but never assert the warning; this pins the deprecate-and-forward contract.
"""

import pytest

import pyopenms

pd = pytest.importorskip("pandas")


def _make_peptide_id_list():
    """A small PeptideIdentificationList with two identified spectra."""
    pep_list = pyopenms.PeptideIdentificationList()
    for rt, mz, score, seq, td in [
        (100.0, 500.25, 50.0, "PEPTIDE", "target"),
        (200.0, 600.30, 75.0, "ANOTHER", "decoy"),
    ]:
        pid = pyopenms.PeptideIdentification()
        pid.setRT(rt)
        pid.setMZ(mz)
        pid.setScoreType("Mascot")
        pid.setIdentifier("test_id")

        hit = pyopenms.PeptideHit()
        hit.setScore(score)
        hit.setCharge(2)
        hit.setSequence(pyopenms.AASequence.fromString(seq))
        hit.setMetaValue("target_decoy", td)
        pid.setHits([hit])

        pep_list.append(pid)
    return pep_list


def test_peptide_identifications_to_df_warns_and_forwards():
    peps = _make_peptide_id_list()

    # the deprecated wrapper must emit a DeprecationWarning naming the replacement
    with pytest.warns(DeprecationWarning, match="to_df"):
        df = pyopenms.peptide_identifications_to_df(peps)

    # ... and forward to PeptideIdentificationList.to_df() (identical result)
    assert isinstance(df, pd.DataFrame)
    ref = peps.to_df()
    assert df.shape == ref.shape
    assert list(df.columns) == list(ref.columns)


def test_update_scores_from_df_warns_and_forwards():
    peps = _make_peptide_id_list()
    df = peps.to_df()
    assert "Mascot" in df.columns          # the score column
    df["Mascot"] = [10.0, 20.0]            # overwrite the scores

    # the deprecated wrapper must emit a DeprecationWarning naming the replacement
    with pytest.warns(DeprecationWarning, match="update_scores_from_df"):
        pyopenms.update_scores_from_df(peps, df, "Mascot")
    via_wrapper = [peps[i].getHits()[0].getScore() for i in range(peps.size())]

    # ... and forward the write-back, matching the direct (non-deprecated) call
    peps_direct = _make_peptide_id_list()
    df_direct = peps_direct.to_df()
    df_direct["Mascot"] = [10.0, 20.0]
    peps_direct.update_scores_from_df(df_direct, "Mascot")
    via_direct = [peps_direct[i].getHits()[0].getScore() for i in range(peps_direct.size())]

    assert via_wrapper == via_direct == [10.0, 20.0]


if __name__ == "__main__":
    test_peptide_identifications_to_df_warns_and_forwards()
    test_update_scores_from_df_warns_and_forwards()
    print("All deprecated-dataframe-shim tests passed!")
