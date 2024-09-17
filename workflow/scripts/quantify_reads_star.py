__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

from snakemake.script import snakemake
from snakemake.shell import shell


extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = [snakemake.input.fq1] if isinstance(snakemake.input.fq1, str) else snakemake.input.fq1
fq2 = snakemake.input.get("fq2")
if fq2:
    fq2 = [snakemake.input.fq2] if isinstance(snakemake.input.fq2, str) else snakemake.input.fq2
    assert len(fq1) == len(fq2), "input-> equal number of files required for fq1 and fq2"
input_str_fq1 = ",".join(fq1)
input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
input_str = " ".join([input_str_fq1, input_str_fq2])


if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand gunzip -c"
elif fq1[0].endswith(".bz2"):
    readcmd = "--readFilesCommand bunzip2 -c"
else:
    readcmd = ""

out_unmapped = snakemake.output.get("unmapped", "")
if out_unmapped:
    out_unmapped = "--outReadsUnmapped Fastx"

index = snakemake.input.get("idx")
if not index:
    index = snakemake.params.get("idx", "")


if "--outSAMtype BAM SortedByCoordinate" in extra:
    stdout = "BAM_SortedByCoordinate"
elif "BAM Unsorted" in extra:
    stdout = "BAM_Unsorted"
else:
    stdout = "SAM"

shell(
    "STAR "
    " --runThreadN {snakemake.threads}"
    " --genomeDir {index}"
    " --readFilesIn {input_str}"
    " {readcmd}"
    " {extra}"
    " {out_unmapped}"
    " --outTmpDir {snakemake.output.output_path}/STARtmp"
    " --outFileNamePrefix {snakemake.output.output_path}/"
    " --outStd {stdout}"
    " > {snakemake.output.aln}"
)
