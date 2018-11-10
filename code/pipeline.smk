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

RAM=config["RAM"]

OUT_DIR=config["PROJECT_DIR"]

rule all:
    input:
        expand(OUT_DIR + "/metrics/{sample}.familysize.txt", sample=SAMPLES),
        expand(OUT_DIR + "/metrics/{sample}.dedup.metrics.txt", sample=SAMPLES),
        expand(OUT_DIR + "/metrics/{sample}.HS.metrics.txt", sample=SAMPLES),
        expand(OUT_DIR + "/metrics/{sample}.HS.metrics.per_target.txt", sample=SAMPLES),
        expand(OUT_DIR + "/bam/{sample}.consensus.dedup.bam",sample=SAMPLES)

rule TrimFastq:
    input:
        fq_one=OUT_DIR + "/fastq/{sample}_" + LANE + "_" + FQ_ONE + "_001.fastq.gz",
        fq_two=OUT_DIR + "/fastq/{sample}_" + LANE + "_" + FQ_TWO + "_001.fastq.gz"
    output:
        fq_one=OUT_DIR + "/trimmed/{sample}.trimmed_" + FQ_ONE + ".fastq.gz",
        fq_two=OUT_DIR + "/trimmed/{sample}.trimmed_" + FQ_TWO + ".fastq.gz"
    message: "Trimming adapters from fastqs on the following files {input}."
    shell:
        "cutadapt "
        "-a {A_FWD} -A {A_REV} -m 40 "
        "-o {output.fq_one} -p {output.fq_two} "
        "{input.fq_one} {input.fq_two}"

rule AlignFastq:
    input:
        fq_one=OUT_DIR + "/trimmed/{sample}.trimmed_" + FQ_ONE + ".fastq.gz",
        fq_two=OUT_DIR + "/trimmed/{sample}.trimmed_" + FQ_TWO + ".fastq.gz"
    output:
        OUT_DIR + "/bam/{sample}.bam"
    threads: 2
    message: "Aligning fastqs with {threads} threads on the following files {input}."
    shell:
        "bwa mem "
        "-t {threads} -M "
        "-R \"@RG\\tID:{RG}\\tSM:{RG}\" "
        "{REF} {input.fq_one} {input.fq_two} | "
        "samtools view -bh - | "
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} "
        "SortSam "
        "I=/dev/stdin O={output} SO=queryname"

rule SetMateInformation:
    input:
        OUT_DIR + "/bam/{sample}.bam"
    output:
        OUT_DIR + "/bam/{sample}.mate.bam"
    message: "Setting mate information on the following files {input}."
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
    message: "Annotating bam with UMIs with {threads} threads on the following files {input}."
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {FGBIO_JAR} "
        "AnnotateBamWithUmis "
        "-i {input.bam} -f {input.fq_umi} "
        "-o {output} "

rule GroupReadsByUmi:
    input:
        OUT_DIR + "/bam/{sample}.umi.bam"
    output:
        bam=OUT_DIR + "/bam/{sample}.grouped.bam",
        metrics=OUT_DIR + "/metrics/{sample}.familysize.txt"
    message: "Grouping reads by UMI on the following files {input}.\n" +
             "Minimum mapping quality = 10.\n" +
             "The allowable number of edits between UMIs = 2."
    threads: 2
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {FGBIO_JAR} "
        "--compression=0 "
        "GroupReadsByUmi "
        "-i {input} -o /dev/stdout "
        "-f {output.metrics} "
        "-m 10 -s Edit -e 2 | "
        "samtools view -b "
        "-@ {threads} "
        "-o {output.bam} "

rule CallConsensus:
    input:
        OUT_DIR + "/bam/{sample}.grouped.bam"
    output:
        OUT_DIR + "/bam/{sample}.consensus.bam"
    message: "Calling consensus on the following files {input}.\n" +
             "Minimum number of reads per UMI = 1."
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {FGBIO_JAR} "
        "CallMolecularConsensusReads "
        "-i {input} -o {output} "
        "-R {RG} --min-reads=1 "

rule AlignConsensusBam:
    input:
        OUT_DIR + "/bam/{sample}.consensus.bam"
    output:
        OUT_DIR + "/bam/{sample}.consensus.align.bam"
    threads: 2
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} SamToFastq "
        "I={input} "
        "FASTQ=/dev/stdout "
        "INTERLEAVE=true NON_PF=true | "
        "bwa mem "
        "-t {threads} -p "
        "-R \"@RG\\tID:{RG}\\tSM:{RG}\" "
        "{REF} /dev/stdin | "
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} "
        "MergeBamAlignment "
        "ALIGNED=/dev/stdin "
        "UNMAPPED={input} "
        "O={output} "
        "R={REF} "
        "ALIGNED_READS_ONLY=true "
        "INCLUDE_SECONDARY_ALIGNMENTS=false "
        "PRIMARY_ALIGNMENT_STRATEGY=BestMapq "
        "SORT_ORDER=queryname "

rule FilterConsensusReads:
    input:
        OUT_DIR + "/bam/{sample}.consensus.align.bam"
    output:
        OUT_DIR + "/bam/{sample}.consensus.filter.bam"
    threads: 2
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {FGBIO_JAR} "
        "--compression=0 "
        "FilterConsensusReads "
        "-i {input} "
        "-o /dev/stdout "
        "-r {REF} "
        " -M 1 -E 0.05 -e 0.5 -N 20 -n 0.2 | "
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} "
        "SortSam "
        "I=/dev/stdin O={output} "
        "SO=coordinate"

rule ClipOverlappingReads:
    input:
        OUT_DIR + "/bam/{sample}.consensus.filter.bam"
    output:
        OUT_DIR + "/bam/{sample}.consensus.clip.bam"
    threads: 2
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {FGBIO_JAR} "
        "--compression=0 "
        "ClipBam "
        "-i {input} -o /dev/stdout "
        "-r {REF} "
        "--clip-overlapping-reads=true | "
        "samtools view -b "
        "-@ {threads} "
        "-o {output} -"

rule MarkDuplicates:
    input:
        OUT_DIR + "/bam/{sample}.consensus.clip.bam"
    output:
        bam=OUT_DIR + "/bam/{sample}.consensus.dedup.bam",
        metrics=OUT_DIR + "/metrics/{sample}.dedup.metrics.txt"
    shell:
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} "
        "MarkDuplicates "
        "I={input} O=/dev/stdout "
    	"M={output.metrics} "
    	"REMOVE_SEQUENCING_DUPLICATES=true | "
        "java -Xmx{RAM}g -XX:-UseParallelGC -jar {PICARD_JAR} "
        "SortSam "
        "I=/dev/stdin O={output.bam} "
        "SO=coordinate"

rule CollectHsMetrics:
    input:
        OUT_DIR + "/bam/{sample}.consensus.dedup.bam"
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
    	"COVERAGE_CAP=10000 "
