# PAIRDIAG_PDEC: where the residual 0.54 pair asymmetry sits (diagnostic, PLAN_PDEC D design data)

Written 19 Nov 01:21. This uses only the PLAN_PDEC D outputs, after the D verdict was committed. It changes no verdict.
The script is `pairdiag.py` and its output is in `PAIRDIAG.txt`.

**Quantity:** P(entrapment target beats its own predicted decoy), among pairs whose better member has q ≤ 0.05. It is
pooled over both arms and 5 seeds, n = 11,714. Overall: **0.545**.

| split | P(win) (n) |
|---|---|
| **decoy and target have the same fragment count** | **0.535** (9,985) |
| decoy has fewer fragments than its target | 0.559 (893) |
| decoy has more fragments than its target | **0.648** (836) |
| target has fewer than 6 fragments | **0.645** (1,163) |
| target has 6 fragments | 0.534 (10,551) |
| charge 2 / 3 / 4 | 0.542 / 0.545 / 0.561 |
| length ≤ 9 / 10-12 / 13-16 / 17+ | 0.537 / 0.528 / 0.552 / 0.565 |
| internal K/R no / yes | 0.546 / 0.542 |

## Reading
1. **Mismatched fragment counts are the largest single contributor.**
   - 15 % of pairs have a decoy with a different number of fragments than its target. They win at 0.56-0.65.
   - The 1,163 targets with fewer than 6 fragments win 0.645. Such a target's decoy usually has 6 fragments, so the pair
     is scored on unequal footing.
   - Why that favours the target is not tested here. A likely reason is that the scorer's features depend on the
     fragment count, so a k < 6 target against a 6-fragment decoy is not an exchangeable comparison.
2. **Even with equal counts, P is 0.535 (±0.005).** That is still above 0.50, so fragment count is not the whole story.
   The rest is not located. Candidates: longer peptides (0.565) and the RT-shift quartile pattern, which is
   non-monotone and not interpretable yet.
3. **Implication:** a decoy should get exactly as many fragments as its target, neither more nor fewer. That is what
   `OpenSwathDecoyGenerator -predict_decoys` does by construction ("as many as its target has"), and `-predict_keep_short`
   relaxes it only downwards.
   - So PLAN_ODP, which uses OpenMS's decoys, avoids the "decoy more" class, whose P is 0.648.
   - DIALibGen predicts decoys independently, so in PLAN_PDEC 7.1 % of pairs had a decoy with more fragments.
   - This is a testable prediction for PLAN_ODP: its pair win should sit below PLAN_PDEC's 0.54.

This is not a criterion and nothing here is used to retune PLAN_PDEC.
