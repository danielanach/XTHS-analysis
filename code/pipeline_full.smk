SAMPLES=config["SAMPLES"]
LANE=config["LANE"] #Need to support multiple lane
RG=config["RG_NAME"]

FQ_ONE=config["FQ_READ_ONE"]
FQ_TWO=config["FQ_READ_TWO"]
FQ_UMI=config["FQ_READ_UMI"]

A_FWD=config["ADAPTER_FWD"]
A_REV=config["ADAPTER_REV"]

REF=config["REF_GENOME"]
BED=config["BED"]
INTERVAL_LIST=config["INTERVAL_LIST"]

FGBIO_JAR=config["FGBIO_JAR"]
PICARD_JAR=config["PICARD_JAR"]

THREADS=config["THREADS"]
RAM=config["RAM"]

OUT_DIR=config["PROJECT_DIR"]

rule all:
    input:
        expand(OUT_DIR + "/metrics/{sample}.familysize.txt", sample=SAMPLES),
        expand(OUT_DIR + "/metrics/{sample}.dedup.metrics.txt", sample=SAMPLES),
        expand(OUT_DIR + "/metrics/{sample}.HS.metrics.txt", sample=SAMPLES),
	expand(OUT_DIR + "/metrics/{sample}.HS.metrics.per_target.txt", sample=SAMPLES),
        expand(OUT_DIR + "/bam/final/{sample}.sort.consensus.dedup.bam",sample=SAMPLES)

rule TrimFastq:
    input:
        fq_one=OUT_DIR + "/fastq/{sample}_" + LANE + "_" + FQ_ONE + "_001.fastq.gz",
        fq_two=OUT_DIR + "/fastq/{sample}_" + LANE + "_" + FQ_TWO + "_001.fastq.gz"
    output:
        fq_one=OUT_DIR + "/trimmed/{sample}.trimmed_" + FQ_ONE + ".fastq.gz",
        fq_two=OUT_DIR + "/trimmed/{sample}.trimmed_" + FQ_TWO + ".fastq.gz"
    shell:
        "cutadapt "
        "-a {A_FWD} -A {A_REV} "
        "-m 40 "
	"-o {output.fq_one} -p {output.fq_two} "
        "{input.fq_one} {input.fq_two}"

rule AlignFastq:
    input:
        fq_one=OUT_DIR + "/trimmed/{sample}.trimmed_" + FQ_ONE + ".fastq.gz",
        fq_two=OUT_DIR + "/trimmed/{sample}.trimmed_" + FQ_TWO + ".fastq.gz"
    output:
        OUT_DIR + "/bam/{sample}.bam"
    shell:
        "bwa mem "
        "-t {THREADS} -M "
        "-R \"@RG\\tID:{RG}\\tSM:{RG}\" "
        "{REF} {input.fq_one} {input.fq_two} "
        "| samtools view -bh - > {output}"

rule Sort:
    input:
        OUT_DIR + "/bam/{sample}.bam"
    output:
        OUT_DIR + "/bam/{sample}.sort.bam"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} SortSam "
        "I={input} O={output} "
        "SO=queryname"

rule SetMateInformation:
    input:
        OUT_DIR + "/bam/{sample}.sort.bam"
    output:
        OUT_DIR + "/bam/{sample}.mate.bam"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {FGBIO_JAR} "
        "SetMateInformation "
        "-i {input} -o {output} "

rule AnnotateBam:
    input:
        bam=OUT_DIR + "/bam/{sample}.mate.bam",
        fq_umi=OUT_DIR + "/fastq/{sample}_" + LANE + "_" + FQ_UMI + "_001.fastq.gz"
    output:
        OUT_DIR + "/bam/{sample}.umi.bam"
    shell:
        "java -Xmx{RAM}g -jar {FGBIO_JAR} "
        "AnnotateBamWithUmis "
        "-i {input.bam} -f {input.fq_umi} "
        "-o {output} "

rule GroupReadsByUmi:
    input:
        OUT_DIR + "/bam/{sample}.umi.bam"
    output:
        bam=OUT_DIR + "/bam/{sample}.grouped.bam",
        metrics=OUT_DIR + "/metrics/{sample}.familysize.txt"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {FGBIO_JAR} "
        "GroupReadsByUmi "
        "-i {input} -o {output.bam} "
        "-f {output.metrics} "
	    "-m 10 -s Edit -e 5 "

rule CallMolecularConsensus:
    input:
        OUT_DIR + "/bam/{sample}.grouped.bam"
    output:
        OUT_DIR + "/bam/{sample}.consensus.bam"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {FGBIO_JAR} "
        "CallMolecularConsensusReads "
        "-i {input} -o {output} "
        "--min-reads=1 "

rule BamToFastq:
    input:
        OUT_DIR + "/bam/{sample}.consensus.bam"
    output:
        fq_one=OUT_DIR + "/fastq/{sample}.consensus_" + FQ_ONE + ".fastq",
        fq_two=OUT_DIR + "/fastq/{sample}.consensus_" + FQ_TWO + ".fastq"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} SamToFastq "
        "I={input} "
        "FASTQ={output.fq_one} "
        "SECOND_END_FASTQ={output.fq_two} "

rule AlignConsensusFastq:
    input:
        fq_one=OUT_DIR + "/fastq/{sample}.consensus_" + FQ_ONE + ".fastq",
        fq_two=OUT_DIR + "/fastq/{sample}.consensus_" + FQ_TWO + ".fastq"
    output:
        OUT_DIR + "/bam/{sample}.consensus.aligned.bam"
    shell:
        "bwa mem "
        "-t {THREADS} -M "
        "-R \"@RG\\tID:{RG}\\tSM:{RG}\" "
        "{REF} {input.fq_one} {input.fq_two} "
        "| samtools view -bh - > {output}"

rule SortConsensus:
    input:
        bam=OUT_DIR + "/bam/{sample}.consensus.aligned.bam"
    output:
        bam=OUT_DIR + "/bam/final/{sample}.consensus.aligned.sort.bam"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} SortSam "
        "I={input.bam} O={output.bam} "
        "SO=queryname"

rule MarkDuplicates:
    input:
        OUT_DIR + "/bam/final/{sample}.consensus.aligned.sort.bam"
    output:
        bam=OUT_DIR + "/bam/final/{sample}.consensus.dedup.bam",
        metrics=OUT_DIR + "/metrics/{sample}.dedup.metrics.txt"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} "
        "MarkDuplicates "
        "I={input} O={output.bam} "
    	"M={output.metrics} "
    	"REMOVE_SEQUENCING_DUPLICATES=true "

rule CollectHsMetrics:
    input:
        OUT_DIR + "/bam/final/{sample}.consensus.dedup.bam"
    output:
        metrics=OUT_DIR + "/metrics/{sample}.HS.metrics.txt",
	    per_t_metrics=OUT_DIR + "/metrics/{sample}.HS.metrics.per_target.txt"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} "
        "CollectHsMetrics "
	    "I={input} O={output.metrics} "
    	"R={REF} "
    	"TARGET_INTERVALS={INTERVAL_LIST} "
    	"BAIT_INTERVALS={INTERVAL_LIST} "
    	"PER_TARGET_COVERAGE={output.per_t_metrics} "
    	"COVERAGE_CAP=10000"

rule SortDuplicated:
    input:
        OUT_DIR + "/bam/final/{sample}.consensus.dedup.bam"
    output:
        OUT_DIR + "/bam/final/{sample}.sort.consensus.dedup.bam"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} SortSam "
        "I={input} O={output} "
        "SO=coordinate"
