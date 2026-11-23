# BYTECHECK_FASTINT2: second interference resolver speed-up is output-identical (23 Nov 2026)

**Code.** OpenMS `3b1d90ed67` (local, not pushed), on top of `17430ba8eb`. Only `OpenDIAInterference.cpp` changed:
- `sharedFragments` uses one sorted merge of both fragment lists instead of a binary search per fragment;
- the nearest window per claimant precursor is now kept on the fly, so the code no longer sorts every partner window.

**Container.** `opendia-pd10-fastint2.sif`, sha256 `81aa2d76…`; the OpenDIA md5 `625e7e9a…` matches the local build. `ctest -R '^OpenDIA' -E TOPP`: 14/14 pass.

**Replay** (resolver input dumped from D r≈1.56 V7M6, 3.88 M windows, laptop, 8 pinned threads):
- the output hash `d36cec9099cc9a0a` and all link statistics are identical;
- 67-79 s before, 59-60 s after.

**Byte-check on codon** (`fastint2-20261123`, same references as `gbmfix-20261123`): **9/9 SAME**.

| Cell | Reference | Files | Wall | Peak RSS | Resolver, pass 1 / pass 2 |
|---|---|---|---:|---:|---:|
| D156 (D r≈1.56 def) | `gbmdef-20261122/D_def_s0` | 3/3 SAME | 5:36 | 9.3 GB | 7.9 / 6.3 s |
| D7V (D r7 V7M6) | `v7m6t-20261122/D_r7_V7M6_1` | 3/3 SAME | 25:13 | 44.4 GB | 230 / 187 s |
| H7 (H r7 def) | `entlib-gbt-20261122/runs/H_gbt_s0` | 3/3 SAME | 23:27 | 25.7 GB | 67 / 64 s |

**Reading.**
- The D7V resolver time fell from 402 / 335 s (RESULT_V7M6T2, `fastint`) to 230 / 187 s. That saves about 320 s on a run of about 2,600 s.
- The nodes differ (standard codon here, Intel 6336Y there), so this is not a matched timing.
- Any V7M6 timing claim needs a new frozen plan.
