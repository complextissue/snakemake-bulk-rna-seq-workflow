"""
Rule for downloading raw sequencing data from ENA
"""


def construct_ena_url(wildcards):
    """
    Construct ENA HTTPS URL from SRR accession
    Format: https://ftp.sra.ebi.ac.uk/vol1/fastq/SRRXXX/0YZ/SRRZZZZZZ/SRRZZZZZZ_R.fastq.gz
    where XXX = first 6 chars, YZ = last 2 digits, ZZZZZZ = full accession, R = read number
    """
    srr = wildcards.sample_id
    # Extract the numeric part (e.g., "27174816" from "SRR27174816")
    srr_num = srr.replace("SRR", "")
    # Get last 2 digits for the folder structure (e.g., "16")
    last_two = "0" + srr_num[-2:]
    # Get first 6 characters of SRR (e.g., "SRR271")
    first_six = srr[:6]

    # Determine if paired-end or single-end
    if IS_PAIRED:
        read = wildcards.read
        return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{first_six}/{last_two}/{srr}/{srr}_{read}.fastq.gz"
    else:
        return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{first_six}/{last_two}/{srr}/{srr}.fastq.gz"


if IS_PAIRED:

    rule download_reads_ena:
        """
        Download paired-end reads from ENA
        """
        output:
            "resources/reads/raw/{sample_id}_{read}.fastq.gz",
        params:
            url=construct_ena_url,
        log:
            "results/logs/download_reads/{sample_id}_{read}.log",
        resources:
            download_slots=1,
        retries: 3
        shell:
            """
            wget --tries=3 --waitretry=5 --retry-connrefused --timeout=60 \
                -O {output} {params.url} > {log} 2>&1 || \
            (echo "Failed to download {params.url} after retries" >> {log} && exit 1)
            """

else:

    rule download_reads_ena:
        """
        Download single-end reads from ENA
        """
        output:
            "resources/reads/raw/{sample_id}.fastq.gz",
        params:
            url=construct_ena_url,
        log:
            "results/logs/download_reads/{sample_id}.log",
        resources:
            download_slots=1,
        retries: 3
        shell:
            """
            wget --tries=3 --waitretry=5 --retry-connrefused --timeout=60 \
                -O {output} {params.url} > {log} 2>&1 || \
            (echo "Failed to download {params.url} after retries" >> {log} && exit 1)
            """


rule download_all_reads:
    """
    Download all reads for samples in metadata
    """
    input:
        expand(
            (
                "resources/reads/raw/{sample_id}_{read}.fastq.gz"
                if IS_PAIRED
                else "resources/reads/raw/{sample_id}.fastq.gz"
            ),
            sample_id=samples.index,
            read=reads if IS_PAIRED else [],
        ),
    default_target: True
