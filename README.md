# XTHS analysis pipeline

Analysis pipeline for Agilent SureSelectXT HS sequencing data.

**NOTE: this pipeline is NOT ready for use**

# Additional software:
fgbio

```
git clone https://github.com/fulcrumgenomics/fgbio.git
sbt assembly
cd fgbio/target/scala-2.12
mv fgbio-<version>.jar fgbio.jar
```

VarDict
```
git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git
```

Shareable Project powered by <https://spro.io>

PROJECT_DIR/    
│── fastq    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.R1.fastq    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.R2.fastq    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── sample.R3.fastq     
│── trimmed    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.trimmed.R1.fastq    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.trimmed.R2.fastq    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── sample.trimmed.R3.fastq    
│── bam    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.bam    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.sort.bam    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.sort.mate.bam     
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.umi.bam    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.umi.grouped.bam    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.consensus.bam     
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.consensus.aligned.bam    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── sample.consensus.aligned.clipped.bam    
│── metrics    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.family_size.txt    
│&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── sample.hs.txt    
└── variant_calls    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.snps.vcf    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.snvs.vcf    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.indels.vcf    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;│── sample.sv.vcf    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── sample.cn.vcf    
