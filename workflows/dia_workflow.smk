report: "report.workflow.rst"

# Obtain run_ids from centroided DIA mzML files in dia_data folder
run_ids, = glob_wildcards("data_dia/{run}.mzML.gz")

rule all:
    input: 
        "results/dialignr/merged.osw"

rule openswath:
    input:
        mzml="data_dia/{run}.mzML.gz",
        pqp="data_library/SpyogenesAssayLibrary_decoy.pqp"
    params:
        irt="data_library/strep_iRT_small.TraML"
    output:
        osw=report("results/openswath/{run}.osw", caption="report/openswath.rst", category="Step 1", subcategory="{run}"),
        chrom="results/xics/{run}.sqMass"
    singularity:
        "docker://openms/executables:latest"
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8192
    shell:
        "localtmp=$SLURM_TMPDIR && "
        "cache=$localtmp/$(cat /proc/sys/kernel/random/uuid)/ && mkdir -p $cache && "
        "cp {input.pqp} $cache/$(basename {input.pqp}) && cp {params.irt} $cache/$(basename {params.irt}) && cp {input.mzml} $cache/$(basename {input.mzml}) && "
        "OpenSwathWorkflow -in $cache/$(basename {input.mzml}) -tr $cache/$(basename {input.pqp}) -tr_irt $cache/$(basename {params.irt}) -out_osw $cache/$(basename {output.osw}) -out_chrom $cache/$(basename {output.chrom}) -sort_swath_maps -threads {threads} -readOptions cacheWorkingInMemory -tempDirectory $cache -min_upper_edge_dist 1 -mz_extraction_window 75 -mz_extraction_window_unit ppm -mz_extraction_window_ms1 35 -mz_extraction_window_ms1_unit ppm -Scoring:DIAScoring:dia_extraction_window 75 -Scoring:DIAScoring:dia_extraction_unit ppm -enable_ms1 true -rt_extraction_window 600 -extra_rt_extraction_window 100 -mz_correction_function quadratic_regression_delta_ppm -Scoring:TransitionGroupPicker:background_subtraction exact -Scoring:TransitionGroupPicker:PeakIntegrator:baseline_type vertical_division_min -Scoring:Scores:use_ms1_mi true -Scoring:Scores:use_mi_score true -Scoring:TransitionGroupPicker:compute_total_mi -Scoring:Scores:use_total_mi_score -batchSize 5500 -min_rsq 0.6 && "
        "mv $cache/$(basename {output.osw}) {output.osw} && "
        "mv $cache/$(basename {output.chrom}) {output.chrom} && "
        "rm -rf $cache"

rule merge_osw:
    input:
        template=rules.openswath.input.pqp,
        osws=expand("results/openswath/{run}.osw", run=run_ids),
    output:
        report("results/pyprophet/merged.osw", caption="report/merge_osws.rst", category="Step 2")
    singularity:
        "docker://pyprophet/pyprophet:2.1.10"
    shell:
        "localtmp=$SLURM_TMPDIR && "
        "cache=$localtmp/$(cat /proc/sys/kernel/random/uuid)/ && mkdir -p $cache && "
        "cp {input} $cache && "
        "pyprophet merge --template $cache/$(basename {input.template}) --out $cache/$(basename {output}) $cache/*.osw && "
        "mv $cache/$(basename {output}) {output} && "
        "rm -rf $cache"

rule pyprophet_scoring:
    input:
        rules.merge_osw.output,
    output:
        report("results/pyprophet/scored/merged.osw", caption="report/pyprophet_scoring.rst", category="Step 3")
    singularity:
        "docker://pyprophet/pyprophet:2.1.10"
    threads:10
    shell:
        "localtmp=$SLURM_TMPDIR && "
        "cache=$localtmp/$(cat /proc/sys/kernel/random/uuid)/ && mkdir -p $cache && "
        "cp {input} $cache && "
        "pyprophet score --classifier=LDA --in $cache/$(basename {input}) --level=ms1ms2 --xeval_num_iter=3 --ss_initial_fdr=0.2 --ss_iteration_fdr=0.01 --threads={threads} && "
        "mv $cache/$(basename {output}) {output} && "
        "mv $cache/*.pdf results/pyprophet/scored/ && "
        "rm -rf $cache"

rule pyprophet_contexts:
    input:
        rules.pyprophet_scoring.output
    output:
        report("results/osw/merged.osw", caption="report/pyprophet_contexts.rst", category="Step 4")
    singularity:
        "docker://pyprophet/pyprophet:2.1.10"
    shell:
        "localtmp=$SLURM_TMPDIR && "
        "cache=$localtmp/$(cat /proc/sys/kernel/random/uuid)/ && mkdir -p $cache && "
        "cp {input} $cache && "
        "pyprophet peptide --context=experiment-wide --in $cache/$(basename {input}) && "
        "pyprophet peptide --context=global --in $cache/$(basename {input}) && "
        "pyprophet protein --context=global --in $cache/$(basename {input}) && "
        "mv $cache/$(basename {output}) {output} && "
        "mv $cache/*.pdf results/osw/ && "
        "rm -rf $cache"

rule dialignr:
    input:
        rules.pyprophet_contexts.output
    output:
        tsv="dialignr.tsv",
        osw=report("results/dialignr/merged.osw", caption="report/dialignr.rst", category="Step 5")
    params:
        "results/"
    singularity:
        "docker://singjust/dialignr:alignment_map"
    threads: 30
    shell:
        "localtmp=$SLURM_TMPDIR && "
        "cache=$localtmp/$(cat /proc/sys/kernel/random/uuid)/ && mkdir -p $cache && "
        #"for f in xics/*.sqMass; do mv -v --  '$f' '${{f%.sqMass}}.chrom.sqMass';done && "
        "bash rename_xics.sh && "
        "cp -R {params}/* $cache && "
        "Rscript --verbose /alignTargetedRunsParallel_cli.R --dataPath=$cache --oswMerged=TRUE --params=context:experiment-wide,maxPeptideFdr:0.05,maxFdrQuery:0.1,baseSubtraction:TRUE,globalAlignment:loess --applyFun=BiocParallel::bplapply --regBioCP='BiocParallel::register(BiocParallel::MulticoreParam(workers={threads},progressbar=TRUE))' &> dialignr.log && "
        # "mv dialignr.tsv {output.tsv} && "
        "mv $cache/$(basename {output.osw}) {output.osw} && "
        "rm -rf $cache"