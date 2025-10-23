def get_fastp_extra_params(wildcards):
    """Build fastp extra parameters from config."""
    cfg = config["check_quality"]["fastp"]
    params = [
        f"-q {cfg['qualified_quality_phred']}",
        f"-u {cfg['unqualified_percent_limit']}",
        f"-n {cfg['n_base_limit']}",
        f"-e {cfg['average_qual']}",
        f"-l {cfg['length_required']}",
        f"--length_limit {cfg['length_limit']}",
    ]

    # Conditional parameters
    if cfg["trim_poly_g"]:
        params.extend(["-g", f"--poly_g_min_len {cfg['poly_g_min_len']}"])
    if cfg["cut_tail"]:
        params.extend(
            [
                "-3",
                f"--cut_tail_window_size {cfg['cut_tail_window_size']}",
                f"--cut_tail_mean_quality {cfg['cut_tail_mean_quality']}",
            ]
        )
    if cfg["cut_front"]:
        params.extend(
            [
                "-5",
                f"--cut_front_window_size {cfg['cut_front_window_size']}",
                f"--cut_front_mean_quality {cfg['cut_front_mean_quality']}",
            ]
        )
    if cfg["low_complexity_filter"]:
        params.extend(["-y", f"-Y {cfg['complexity_threshold']}"])
    if cfg["overrepresentation_analysis"]:
        params.append("-p")

    return " ".join(params)


def get_fastp_input(wildcards):
    """Get fastp input files based on read type."""
    if IS_PAIRED:
        return [
            f"resources/reads/raw/{wildcards.sample_id}_1.fastq.gz",
            f"resources/reads/raw/{wildcards.sample_id}_2.fastq.gz",
        ]
    else:
        return [f"resources/reads/raw/{wildcards.sample_id}.fastq.gz"]


def get_fastp_output(wildcards):
    """Get fastp output files based on read type."""
    if IS_PAIRED:
        return [
            f"resources/reads/trimmed/{wildcards.sample_id}_1.fastq.gz",
            f"resources/reads/trimmed/{wildcards.sample_id}_2.fastq.gz",
        ]
    else:
        return [f"resources/reads/trimmed/{wildcards.sample_id}.fastq.gz"]


# FastQC rules - must remain conditional due to different wildcard patterns
if IS_PAIRED:

    rule run_fastqc:
        input:
            "resources/reads/raw/{sample_id}_{read}.fastq.gz",
        output:
            html="results/plots/check_quality/run_fastqc/{sample_id}_{read}.html",
            zip="results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
        log:
            "results/logs/check_quality/run_fastqc/{sample_id}_{read}.log",
        benchmark:
            "results/benchmarks/check_quality/run_fastqc/{sample_id}_{read}.tsv"
        params:
            "--quiet",
        threads: config["check_quality"]["threads"]
        resources:
            mem_mb=config["check_quality"]["memory_limit_mb"],
            runtime=30,  # 30 minutes
        wrapper:
            "v7.8.1/bio/fastqc"

else:

    rule run_fastqc:
        input:
            "resources/reads/raw/{sample_id}.fastq.gz",
        output:
            html="results/plots/check_quality/run_fastqc/{sample_id}.html",
            zip="results/plots/check_quality/run_fastqc/{sample_id}_fastqc.zip",
        log:
            "results/logs/check_quality/run_fastqc/{sample_id}.log",
        benchmark:
            "results/benchmarks/check_quality/run_fastqc/{sample_id}.tsv"
        params:
            "--quiet",
        threads: config["check_quality"]["threads"]
        resources:
            mem_mb=config["check_quality"]["memory_limit_mb"],
            runtime=30,  # 30 minutes
        wrapper:
            "v7.8.1/bio/fastqc"


# Fastp rule - unified with input functions (but output must still be conditional)
if IS_PAIRED:

    rule run_fastp:
        input:
            sample=get_fastp_input,
        output:
            trimmed=[
                "resources/reads/trimmed/{sample_id}_1.fastq.gz",
                "resources/reads/trimmed/{sample_id}_2.fastq.gz",
            ],
            failed="resources/reads/trimmed/{sample_id}.failed.fastq",
            html="results/plots/check_quality/run_fastp/{sample_id}.html",
            json="results/plots/check_quality/run_fastp/{sample_id}.json",
        log:
            "results/logs/check_quality/run_fastp/{sample_id}.log",
        benchmark:
            "results/benchmarks/check_quality/run_fastp/{sample_id}.tsv"
        params:
            extra=get_fastp_extra_params,
        threads: config["check_quality"]["threads"]
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000,
            runtime=lambda wildcards, attempt: attempt * 120,
        wrapper:
            "v7.8.1/bio/fastp"

else:

    rule run_fastp:
        input:
            sample=get_fastp_input,
        output:
            trimmed=["resources/reads/trimmed/{sample_id}.fastq.gz"],
            failed="resources/reads/trimmed/{sample_id}.failed.fastq",
            html="results/plots/check_quality/run_fastp/{sample_id}.html",
            json="results/plots/check_quality/run_fastp/{sample_id}.json",
        log:
            "results/logs/check_quality/run_fastp/{sample_id}.log",
        benchmark:
            "results/benchmarks/check_quality/run_fastp/{sample_id}.tsv"
        params:
            extra=get_fastp_extra_params,
        threads: config["check_quality"]["threads"]
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000,
            runtime=lambda wildcards, attempt: attempt * 120,
        wrapper:
            "v7.8.1/bio/fastp"


# MultiQC rules - unified with input functions
def get_multiqc_inputs(wildcards):
    """Generate input list for MultiQC based on configuration and read type."""
    inputs = []

    # Always include fastp outputs
    inputs.extend(
        expand(
            "results/plots/check_quality/run_fastp/{sample_id}.json",
            sample_id=samples.index,
        )
    )

    # Add fastqc outputs if enabled
    if config["check_quality"]["run_fastqc"]:
        if IS_PAIRED:
            inputs.extend(
                expand(
                    "results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
                    sample_id=samples.index,
                    read=reads,
                )
            )
        else:
            inputs.extend(
                expand(
                    "results/plots/check_quality/run_fastqc/{sample_id}_fastqc.zip",
                    sample_id=samples.index,
                )
            )

    return inputs


rule run_multiqc:
    input:
        get_multiqc_inputs,
    output:
        "results/plots/check_quality/run_multiqc/report.html",
    log:
        "results/logs/check_quality/run_multiqc.log",
    benchmark:
        "results/benchmarks/check_quality/run_multiqc.tsv"
    params:
        extra="",
        use_input_files_only=True,
    resources:
        mem_mb=2000,
        runtime=30,
    wrapper:
        "v7.8.1/bio/multiqc"


def get_multiqc_after_salmon_inputs(wildcards):
    """Generate input list for MultiQC after salmon based on configuration and read type."""
    inputs = []

    # Always include fastp outputs
    inputs.extend(
        expand(
            "results/plots/check_quality/run_fastp/{sample_id}.json",
            sample_id=samples.index,
        )
    )

    # Add fastqc outputs if enabled
    if config["check_quality"]["run_fastqc"]:
        if IS_PAIRED:
            inputs.extend(
                expand(
                    "results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
                    sample_id=samples.index,
                    read=reads,
                )
            )
        else:
            inputs.extend(
                expand(
                    "results/plots/check_quality/run_fastqc/{sample_id}_fastqc.zip",
                    sample_id=samples.index,
                )
            )

    # Add salmon outputs
    inputs.extend(
        expand(
            [
                "resources/reads/quantified_salmon/{sample_id}/lib_format_counts.json",
                "resources/reads/quantified_salmon/{sample_id}/aux_info/meta_info.json",
                "resources/reads/quantified_salmon/{sample_id}/libParams/flenDist.txt",
            ],
            sample_id=samples.index,
        )
    )

    return inputs


rule run_multiqc_after_salmon:
    input:
        get_multiqc_after_salmon_inputs,
    output:
        "results/plots/check_quality/run_multiqc_after_salmon/report.html",
    log:
        "results/logs/check_quality/run_multiqc_after_salmon.log",
    benchmark:
        "results/benchmarks/check_quality/run_multiqc_after_salmon.tsv"
    params:
        extra="",
        use_input_files_only=True,
    resources:
        mem_mb=2000,
        runtime=30,
    wrapper:
        "v7.8.1/bio/multiqc"
