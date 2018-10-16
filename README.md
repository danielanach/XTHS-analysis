# XTHS analysis pipeline

Analysis pipeline for Agilent SureSelectXT HS sequencing data.

Shareable Project powered by <https://spro.io>

PROJECT_DIR/    
|── fastq    
│   |── sample.R1.fastq    
│   |── sample.R2.fastq    
│   |── sample.R3.fastq     
│── trimmed    
│   |── sample.trimmed.R1.fastq    
│   │── sample.trimmed.R2.fastq    
│   └── sample.trimmed.R3.fastq    
|── bam    
│   |── sample.bam    
│   |── sample.sort.bam    
│   |── sample.sort.mate.bam     
│   |── sample.umi.bam    
│   |── sample.umi.grouped.bam    
│   |── sample.consensus.bam     
│   |── sample.consensus.aligned.bam    
│   └── sample.consensus.aligned.clipped.bam    
|── metrics    
│   |── sample.family_size.txt    
│   └── sample.hs.txt    
└── variant_calls    
    |── sample.snps.vcf    
    |── sample.snvs.vcf    
    |── sample.indels.vcf    
    |── sample.sv.vcf    
    └── sample.cn.vcf    
