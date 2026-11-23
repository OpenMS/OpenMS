# CUDA infrastructure retry, 23 November 2026

## Attempt 1, retained failure

Frozen89ae34a351 with routing amendments through9cb53f57e0. Accepted jobs41517072 (GPU prediction),41517073 (CPU preparation),41517074 (15-search array),41517075 (evaluation). Prediction ran on codon-gpu-005, A10080GB, driver580.178.04. It digested the full unsampled FASTA into2,692,562 eligible peptides,8,649,692 precursors and164,893,268 eligible fragment transitions. These are input/eligibility counts, not successful predicted library output or identifications.

Prediction failed exit13 before inference: ONNX CUDA provider requires `cudaLibraryGetKernel` but bundled CUDA12.6 libcudart does not provide it. No CPU fallback was permitted. Preparation and searches were dependency-cancelled, evaluation failed because no library exists. Archive all first-attempt logs, accounting and directory artifacts under attempt1 before resubmission. They remain in end-to-end cost/failure reporting. No completed search or quantitative/FDR result exists.

## Frozen infrastructure-only retry

Keep the same DIALibGen container/model/binary, same A100/16CPU/384GiB/12h, FASTA, raw data, eligibility, fragments6-20, decoys, seeds and endpoints. Install only NVIDIA's pinned `nvidia-cuda-runtime-cu12==12.8.90` wheel into `cuda_runtime128` and prepend its library directory to the GPU container's LD_LIBRARY_PATH. Hash the actual libcudart binary in `cuda_runtime.sha256` and verify it inside the job before prediction. Do not drop require_gpu or change any statistical rule. Model predictions still run on CUDA; any numerical sensitivity is retained as infrastructure provenance, not assumed byte-identical.

NVIDIA's CUDA12.8 runtime documentation lists the required library-management API. This repair changes runtime linkage, not the scientific recipe. Commit/push this amendment and updated frozen checksums before retry. No additional retry is silently authorized if this fix fails. Preserve attempt1 and new job IDs separately. No search result is replaced or best seed selected.
