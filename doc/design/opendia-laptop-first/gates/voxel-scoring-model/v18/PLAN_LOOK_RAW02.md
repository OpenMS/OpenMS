# V18 look: charge-sibling co-elution, raw02 replication (written before raw02 is read)

Raw01 look (descriptive, labels used only to summarise): `sibling_look.py`. For each D11FWT precursor winner,
take the other charge states of the same library peptide in the population and record the best per-fold
decoy-normalised full D11F score among their candidate windows with apex within 3 s (`sib_z_near`).

Raw01 numbers that motivate this: in own-score band decoy percentile [0.99, 0.9999) with sib_z_near > 2,
human 78, entrapment 0, decoy 1 (D 54). ENT-vs-DEC AUC of sib_z_near 0.508, D-vs-DEC 0.632.

Raw02 rule, same script with raw=raw02, nothing else changed (threshold 2 and 3 s fixed from raw01):
1. In that band with sib_z_near > 2: entrapment <= 3 and decoy <= 3, human >= 40.
2. ENT-vs-DEC AUC of sib_z_near (all winners with a sibling) <= 0.52.
Pass means: worth a frozen V18 plan. It is not a recovery claim and does not change any gate.
