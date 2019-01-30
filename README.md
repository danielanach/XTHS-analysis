# **NOTE: this pipeline is NOT ready for use**

# XTHS analysis pipeline

Analysis pipeline for Agilent SureSelectXT HS sequencing data.


# Additional software:
fgbio

```
git clone https://github.com/fulcrumgenomics/fgbio.git
sbt assembly
cd fgbio/target/scala-2.12
mv fgbio-<version>.jar fgbio.jar
```

gradle

```
wget https://services.gradle.org/distributions/gradle-4.10.2-all.zip
sudo mkdir /opt/gradle
sudo unzip -d /opt/gradle gradle-4.10.2-bin.zip
ls /opt/gradle/gradle-4.10.2
bin  docs  getting-started.html  init.d  lib  LICENSE  media  NOTICE  samples  src
```

VarDictJava

```
git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git
cd VarDictJava
./gradlew clean installDist
```

You can add VarDictJava on PATH by adding this line to .bashrc:

```
export PATH=/path/to/VarDict/bin:$PATH
```

snpSift
```
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

download dbSNP database
```
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
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
