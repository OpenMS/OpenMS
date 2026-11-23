# Scheduler-only pre-execution amendment, 23 November 2026

Initial frozen commit89ae34a351 was pushed before submission. The first submission was rejected by Codon's job-submit policy: explicit GPU partition selection is forbidden and GPU routing must use `--gres`. No job was accepted, no library was generated, and no search/evaluation outcome exists.

Remove only `#SBATCH --partition=gpu` from predict_gpu.sbatch. Retain one A100,16 CPU,384GiB,12h and every scientific/generation/evaluation option unchanged. Let site policy route `--gres=gpu:a100:1`. Re-freeze checksums, commit and push before retry. Preserve original remote submission.log, write retry log separately. No resource increase, changed criterion or replaced seed.

A subsequent `sbatch --test-only` also rejected explicit `bigmem` selection before any job was submitted. Remove `#SBATCH --partition=bigmem` from prepare.sbatch, retaining512GiB and all CPU/decoy settings. Site policy routes by `--mem`. Standard search/evaluation requests passed scheduler-only validation. This second amendment is also committed/pushed before any computation.
