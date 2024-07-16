rule run_fastqc:
    input:
        "resources/reads/{sample_id}_{read}.fastq.gz",
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
        "v3.13.6/bio/fastqc"


rule run_fastp:
    input:
        sample=["resources/reads/{sample_id}_1.fastq.gz", "resources/reads/{sample_id}_2.fastq.gz"],
    output:
        trimmed=["resources/reads/trimmed/{sample_id}_1.fastq.gz", "resources/reads/trimmed/{sample_id}_2.fastq.gz"],
        failed="resources/reads/trimmed/{sample_id}.failed.fastq",
        html="results/plots/check_quality/run_fastp/{sample_id}.html",
        json="results/plots/check_quality/run_fastp/{sample_id}.json",
    log:
        "results/logs/check_quality/run_fastp/{sample_id}.log",
    params:
        extra="",
    threads: config["check_quality"]["threads"]
    wrapper:
        "v3.13.6/bio/fastp"


rule run_multiqc:
    input:
        [
            *expand(
                "results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
                sample_id=samples.index,
                read=reads,
            ),
            *expand(
                "results/plots/check_quality/run_fastp/{sample_id}.json",
                sample_id=samples.index,
            ),
        ],
    output:
        "results/plots/check_quality/run_multiqc/report.html",
    log:
        "results/logs/check_quality/run_multiqc.log",
    params:
        extra="",
        use_input_files_only=True,
    wrapper:
        "v3.13.6/bio/multiqc"


rule run_multiqc_after_salmon:
    input:
        [
            *expand(
                "results/plots/check_quality/run_fastqc/{sample_id}_{read}_fastqc.zip",
                sample_id=samples.index,
                read=reads,
            ),
            *expand(
                "results/plots/check_quality/run_fastp/{sample_id}.json",
                sample_id=samples.index,
            ),
            *expand(
                "resources/reads/quantified/{sample_id}/lib_format_counts.json",
                sample_id=samples.index,
            ),
            *expand(
                "resources/reads/quantified/{sample_id}/aux_info/meta_info.json",
                sample_id=samples.index,
            ),
            *expand(
                "resources/reads/quantified/{sample_id}/libParams/flenDist.txt",
                sample_id=samples.index,
            ),
        ],
    output:
        "results/plots/check_quality/run_multiqc_after_salmon/report.html",
    log:
        "results/logs/check_quality/run_multiqc_after_salmon.log",
    params:
        extra="",
        use_input_files_only=True,
    wrapper:
        "v3.13.6/bio/multiqc"
