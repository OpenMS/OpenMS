# Second dataset: Agilent 6550 Q-TOF (PXD069532, broccoli, 3 + 3 runs)

Inputs: the six ProteoWizard mzML and the FASTA of
archive.openms.de/openms/testfiles/20260717-ProteomicsLFQ-Multi-Instrument/PXD069532-Agilent-6550/ (sha256 below).
Entrapment database: the dataset's 56,805 targets (without its rev_ decoys) → frozen make_entrapment.py → 42657f57…
IDs as for HF-X (#10585 recipe), but 20 ppm precursor tolerance (Q-TOF, as the FragPipe reference;
scripts/id_one_agilent.sh). ProteomicsLFQ as for HF-X (scripts/lfq_arm_agilent.sh, scripts/run_agilent.sh).

| | off | donor_max_q 0.01 | donor_max_q 0.001 |
|---|---|---|---|
| PSMs at Percolator q < 0.01 per run | 505–719 | | |
| precursors in PIP-ECHO | 618 | 618 (369 transferable) | 618 (280 transferable) |
| PIP-ECHO transfer candidates / decoys | 6 / 0 | | |
| transfers kept | 0 | 0 | 0 |

PIP-ECHO needs at least 20 decoy transfers (min_decoys) to estimate the transfer FDR. Without them it drops all
transfers ("The transfer FDR could NOT be controlled ... all 6 transferred (acceptor) feature(s) are being DROPPED").
So no arm transfers anything, and the dataset says nothing about the donor gate. FragPipe reports 498 MBR values on
411 ions for these runs (FragPipe-LFQ-MBR.zip, combined_ion.tsv). The SCIEX TripleTOF set of the same folder was not
run: FragPipe reports only 312 MBR values there (4,540 ions), too few for an entrapment count.

sha256:
e23994b6558e6ca42ea8df6e448bb3e354fa6635ba98dbed80e01eadfb90926e  202505616.mzML
2b421311ee41dca29443a9eb1768f9e2716e607ea42fa04cebdb1835a91a709a  202505617.mzML
dee61d7ec8a3377ada82de78e4377f10b349c22b6fef1e2629e687a254e7d024  202505618.mzML
73da8340cea40896196c5eb0d78a209c9a612dbb7346e36ae9a4b480367bb298  202505622.mzML
ecf4afea9daf99fbf4642bee012bd5148d17892334b9a2716b91614339b53e25  202505623.mzML
2204851948346e585de6e8a7aa256a28a6516f710060ab563b7d1686391d35c6  202505624.mzML
267c80d917b25b7ca560a3bc2a286af0cb1f79b259c9cedd5f0d2684b84338ab  2026-07-17-decoys-contam-GCF_000695525.1_BOL_protein.fas
42657f579b5d8b47cc0f27d1796856d2cb94975f1c93e6bcda4833145fdf133a  BOL_entrap_td.fasta
