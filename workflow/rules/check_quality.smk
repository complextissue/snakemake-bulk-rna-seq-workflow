# Quality control rules that automatically adapt to paired-end or single-end reads

if IS_PAIRED:

    # Paired-end read rules
    rule run_fastqc:
        input:
            "resources/reads/raw/{sample_id}_{read}.fastq.gz",
        output:
            html="results/plots/check_quality/run_fastqc/{sample_id}_{read}.html",
            zip="results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
        log:
            "results/logs/check_quality/run_fastqc/{sample_id}_{read}.log",
        params:
            "--quiet",
        threads: config["check_quality"]["threads"]
        resources:
            mem_mb=config["check_quality"]["memory_limit_mb"],
        wrapper:
            "v7.8.1/bio/fastqc"

    rule run_fastp:
        input:
            sample=[
                "resources/reads/raw/{sample_id}_1.fastq.gz",
                "resources/reads/raw/{sample_id}_2.fastq.gz",
            ],
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
        params:
            extra=lambda wildcards, config=config: " ".join(
                [
                    f"-q {config['check_quality']['fastp']['qualified_quality_phred']}",
                    f"-u {config['check_quality']['fastp']['unqualified_percent_limit']}",
                    f"-n {config['check_quality']['fastp']['n_base_limit']}",
                    f"-e {config['check_quality']['fastp']['average_qual']}",
                    f"-l {config['check_quality']['fastp']['length_required']}",
                    f"--length_limit {config['check_quality']['fastp']['length_limit']}",
                    "-g" if config["check_quality"]["fastp"]["trim_poly_g"] else "",
                    (
                        f"--poly_g_min_len {config['check_quality']['fastp']['poly_g_min_len']}"
                        if config["check_quality"]["fastp"]["trim_poly_g"]
                        else ""
                    ),
                    "-3" if config["check_quality"]["fastp"]["cut_tail"] else "",
                    (
                        f"--cut_tail_window_size {config['check_quality']['fastp']['cut_tail_window_size']}"
                        if config["check_quality"]["fastp"]["cut_tail"]
                        else ""
                    ),
                    (
                        f"--cut_tail_mean_quality {config['check_quality']['fastp']['cut_tail_mean_quality']}"
                        if config["check_quality"]["fastp"]["cut_tail"]
                        else ""
                    ),
                    "-5" if config["check_quality"]["fastp"]["cut_front"] else "",
                    (
                        f"--cut_front_window_size {config['check_quality']['fastp']['cut_front_window_size']}"
                        if config["check_quality"]["fastp"]["cut_front"]
                        else ""
                    ),
                    (
                        f"--cut_front_mean_quality {config['check_quality']['fastp']['cut_front_mean_quality']}"
                        if config["check_quality"]["fastp"]["cut_front"]
                        else ""
                    ),
                    (
                        "-y"
                        if config["check_quality"]["fastp"]["low_complexity_filter"]
                        else ""
                    ),
                    (
                        f"-Y {config['check_quality']['fastp']['complexity_threshold']}"
                        if config["check_quality"]["fastp"]["low_complexity_filter"]
                        else ""
                    ),
                    (
                        "-p"
                        if config["check_quality"]["fastp"][
                            "overrepresentation_analysis"
                        ]
                        else ""
                    ),
                ]
            ),
        threads: config["check_quality"]["threads"]
        wrapper:
            "v7.8.1/bio/fastp"

else:

    # Single-end read rules
    rule run_fastqc:
        input:
            "resources/reads/raw/{sample_id}.fastq.gz",
        output:
            html="results/plots/check_quality/run_fastqc/{sample_id}.html",
            zip="results/plots/check_quality/run_fastqc/{sample_id}_fastqc.zip",
        log:
            "results/logs/check_quality/run_fastqc/{sample_id}.log",
        params:
            "--quiet",
        threads: config["check_quality"]["threads"]
        resources:
            mem_mb=config["check_quality"]["memory_limit_mb"],
        wrapper:
            "v7.8.1/bio/fastqc"

    rule run_fastp:
        input:
            sample=["resources/reads/raw/{sample_id}.fastq.gz"],
        output:
            trimmed=["resources/reads/trimmed/{sample_id}.fastq.gz"],
            failed="resources/reads/trimmed/{sample_id}.failed.fastq",
            html="results/plots/check_quality/run_fastp/{sample_id}.html",
            json="results/plots/check_quality/run_fastp/{sample_id}.json",
        log:
            "results/logs/check_quality/run_fastp/{sample_id}.log",
        params:
            extra=lambda wildcards, config=config: " ".join(
                [
                    f"-q {config['check_quality']['fastp']['qualified_quality_phred']}",
                    f"-u {config['check_quality']['fastp']['unqualified_percent_limit']}",
                    f"-n {config['check_quality']['fastp']['n_base_limit']}",
                    f"-e {config['check_quality']['fastp']['average_qual']}",
                    f"-l {config['check_quality']['fastp']['length_required']}",
                    f"--length_limit {config['check_quality']['fastp']['length_limit']}",
                    "-g" if config["check_quality"]["fastp"]["trim_poly_g"] else "",
                    (
                        f"--poly_g_min_len {config['check_quality']['fastp']['poly_g_min_len']}"
                        if config["check_quality"]["fastp"]["trim_poly_g"]
                        else ""
                    ),
                    "-3" if config["check_quality"]["fastp"]["cut_tail"] else "",
                    (
                        f"--cut_tail_window_size {config['check_quality']['fastp']['cut_tail_window_size']}"
                        if config["check_quality"]["fastp"]["cut_tail"]
                        else ""
                    ),
                    (
                        f"--cut_tail_mean_quality {config['check_quality']['fastp']['cut_tail_mean_quality']}"
                        if config["check_quality"]["fastp"]["cut_tail"]
                        else ""
                    ),
                    "-5" if config["check_quality"]["fastp"]["cut_front"] else "",
                    (
                        f"--cut_front_window_size {config['check_quality']['fastp']['cut_front_window_size']}"
                        if config["check_quality"]["fastp"]["cut_front"]
                        else ""
                    ),
                    (
                        f"--cut_front_mean_quality {config['check_quality']['fastp']['cut_front_mean_quality']}"
                        if config["check_quality"]["fastp"]["cut_front"]
                        else ""
                    ),
                    (
                        "-y"
                        if config["check_quality"]["fastp"]["low_complexity_filter"]
                        else ""
                    ),
                    (
                        f"-Y {config['check_quality']['fastp']['complexity_threshold']}"
                        if config["check_quality"]["fastp"]["low_complexity_filter"]
                        else ""
                    ),
                    (
                        "-p"
                        if config["check_quality"]["fastp"][
                            "overrepresentation_analysis"
                        ]
                        else ""
                    ),
                ]
            ),
        threads: config["check_quality"]["threads"]
        wrapper:
            "v7.8.1/bio/fastp"


# MultiQC rules - adapt input based on run_fastqc config and read type
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
    params:
        extra="",
        use_input_files_only=True,
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
    params:
        extra="",
        use_input_files_only=True,
    wrapper:
        "v7.8.1/bio/multiqc"
