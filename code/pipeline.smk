from snakemake.utils import validate

validate(config)

SAMPLES=config["SAMPLES"]
REF=config["REF_GENOME"]
BED=config["BED"]
AGENT_DIR=config["AGENT_DIR"]

rule all:
    input:
        expand("{sample}.LocatIt.bam",sample=SAMPLES)

rule TrimFastq:
    input:
        fq_one="{sample}_R1_001.fastq.gz",
        fq_two="{sample}_R3_001.fastq.gz"
    output:
    shell:
        "java -jar ${{AGENT_DIR}}/SurecallTrimmer.jar "
        "-fq1 {input.fq_one} -fq2 {input.fq_two} "
        "-xt -out_loc ${{OUT_PATH}}"

rule AlignFastq:
    input:
        "{sample}_R1_001.fastq.gz",
        "{sample}_R3_001.fastq.gz"
    output:
        "{sample}.bam"
    shell:
        "bwa mem -t 2 -M -I 200,100,600,1 "
        "-R @RG\tID:1\tSM:Tissue\tPL:ILLUMINA\tLB:XT "
        "{REF_GENOME} "
        "{input} "

rule ExtractUMI:
    input:
        bam="{sample}.bam",
        umi_fastq="{sample}_R2_.fastq.gz",
        bed=REF_BED
    output:
        "{sample}.LocatIt.bam"
    shell:
        "java –Xmx12G -jar ${AGENT_DIR}/LocatIt.jar "
        "[-X temp_directory] [-t temp_directory for keeps] "
        "-L [-PM:xm,Q:xq..] [-C] [-i] "
        " -IB -OB [-r] [-d NN] [-2] "
        "[-q 25][-m 2][-c 2500][-H sam_header_file] "
        "-b {input.bed} "
        "-o {output} "
        "{input.bam} "
        "{input.umi_fastq} "
