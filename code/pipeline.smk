AGENT_DIR="/Users/DanielaNachmanson/work/AGeNT"
REF_DIR=""

rule all:
    input:
        "sample.bam"

rule FastqToBam:
    input:
    output:
    shell:
    " picard CreateUnmappedBam "



rule ExtractUMI:
    input:
        bam: "sample.bam"
        fastq: "sample.fastq"
        bed: ""
    output:

    shell:
        "java –Xmx12G -jar LocatIt.jar "
        "[-X temp_directory] [-t temp_directory for keeps] "
        "-L [-PM:xm,Q:xq..] [-C] [-i] "
        " -IB -OB [-r] [-d NN] [-2] "
        "[-q 25][-m 2][-c 2500][-H sam_header_file] "
        "-b {input.bed} "
        "-o output_file_name "
        "input_bam_file_name "
        "index2_fastq_file_1[.gz] "
        "[...index2_fastq_file_N[.gz]] "

rule plot:
    input:
        "raw/{dataset}.csv"
    output:
        "plots/{dataset}.pdf"
    shell:
        "somecommand {input} {output}"
