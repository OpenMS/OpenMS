# RNA sequence data

## human_background_rna.fasta

The built-in background database of NucleicAcidSearchEngine (`-background:enable` without `-background:database`).
NASE searches these sequences, with decoys it generates for them, alongside the target database. Spectra whose best
match is a background sequence (or one of its decoys) are left out of the results and of the false discovery rate
calculation. This keeps spectra of abundant human RNAs that the target database lacks from matching decoys or false
targets (the "subset-neighbor search" of Lin et al., J Proteome Res 2021, 20:4153).

Contents (1769 sequences, accessions starting with `BG_`):

| Source | Entries |
| --- | --- |
| GtRNAdb hg38 mature tRNAs (`hg38-mature-tRNAs.fa`, 432 tRNAs) | each with 5' phosphate and the 3' ends CCA, CC, C and none |
| NCBI NC_012920.1 (mitochondrial genome, rCRS) | the 22 annotated tRNAs (5' phosphate, 3' CCA); the genome on both strands |
| NCBI NR_046235.3, NR_023363.1 | 45S pre-rRNA, 5S rRNA |
| NCBI NR_002715.1, NR_001445.2, NR_004430.2, NR_002716.3, NR_003925.1, NR_002756.2, NR_004394.1, NR_004391.1, NR_004392.1, NR_004393.1, NR_001571.2, NR_026703.1, NR_002312.1, NR_003051.3, NR_002999.2 | 7SL, 7SK, U1, U2, U4, U5A, U6, Y1, Y3, Y4, Y5, vault 1-1, RNase P (H1), RNase MRP and scaRNA20 RNAs |

The background is for human samples. It does not contain sample-specific sequence variants or modified forms that
the search settings cannot express; such molecules can still match decoys or false targets.

To rebuild the file (e.g. after a GtRNAdb update), fetch the inputs listed in `build_human_rna_background.py` and run

    python3 build_human_rna_background.py <download dir> human_background_rna.fasta
