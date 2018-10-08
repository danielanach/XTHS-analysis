SAMPLES=config["SAMPLES"]
REF_GENOME=config["REF_GENOME"]
BED=config["BED"]
AGENT_DIR=config["AGENT_DIR"]
FGBIO_DIR="/home/dnachman/XTHS-analysis/code/fgbio"
OUT_DIR=config["OUT_DIR"]

rule all:
    input:
        expand(OUT_DIR + "/bam/{sample}.umi.bam",sample=SAMPLES)

rule TrimFastq:
    input:
        fq_one=OUT_DIR + "/fastq/{sample}_R1_001.fastq.gz",
        fq_two=OUT_DIR + "/fastq/{sample}_R3_001.fastq.gz"
    output:
        OUT_DIR + "/trimmed/{sample}_R1_001.fastq.gz",
        OUT_DIR + "/trimmed/{sample}_R3_001.fastq.gz",
    shell:
        "java -jar {AGENT_DIR}/SurecallTrimmer_v4.0.1.jar "
        "-fq1 {input.fq_one} -fq2 {input.fq_two} "
        "-xt "
        "-qualityTrimming 10 -minFractionRead 50 "
        "-out_loc {OUT_DIR}/trimmed/"

rule AlignFastq:
    input:
        OUT_DIR + "/trimmed/{sample}_R1_001.fastq.gz",
        OUT_DIR + "/trimmed/{sample}_R3_001.fastq.gz"
    output:
        OUT_DIR + "/bam/{sample}.bam"
    shell:
        "bwa mem -t 2 -M -R \"@RG\\tID:tissue\\tSM:tissue\" "
        "{REF_GENOME} {input} | "
        "samtools view -bh - > {output}"

rule Sort:
    input:
        OUT_DIR + "/bam/{sample}.bam"
    output:
        bam=OUT_DIR + "/bam/{sample}.sort.bam"
    shell:
        "samtools sort {input} > {output.bam}"

rule Index:
    input:
        OUT_DIR + "/bam/{sample}.sort.bam"
    output:
        OUT_DIR + "/bam/{sample}.sort.bam.bai"
    shell:
        "samtools index {input}"


rule AnnotateBam:
    input:
        bam=OUT_DIR + "/bam/{sample}.sort.bam",
        umi_fastq=OUT_DIR + "/fastq/{sample}_R2_001.fastq.gz"
    output:
        OUT_DIR + "/bam/{sample}.umi.bam"
    shell:
        "java -Xmx60g -XX:+UseParallelGC -jar {FGBIO_DIR}/target/scala-2.12/fgbio-0.7.0-9b76546-SNAPSHOT.jar "
        "AnnotateBamWithUmis -i {input.bam} "
        "-f {input.umi_fastq} "
        "-o {output} > fgbio_umi.txt"
        "\n"
        "samtools index {output}"


# rule ExtractUMI:
#     input:
#         bam=OUT_DIR + "/bam/{sample}.sort.bam",
#         umi_fastq=OUT_DIR + "/fastq/{sample}_R2_001.fastq.gz"
#     output:
#         OUT_DIR + "/bam/{sample}.LocatIt.bam"
#     shell:
#         "java –Xmx10G -jar {AGENT_DIR}/LocatIt_v4.0.1.jar "
#         "-X /scratch/dana/temp/ -t /scratch/dana/bam/intermediate/ "
#         "-PM:MI,Q:cQ,I:cD -C -IB -OB -r -d 1 "
#         "-b {BED} "
#         "-o {output} "
#         "-l {BED} "
#         "{input.bam} "
#         "{input.umi_fastq} "
#         "\n"
#         "samtools index {output}"
