---
bibliography: references.bib
---

# FastqToDiagnostics

Repository for code and documentation for the Malian Data Science and Bioinformatics Network (MD-BioNet).

## Introduction

Here we describe steps to perform a simple small nucleotide variation calling pipeline for Whole Exome sequencing reads.

The pipline performs read quality control using *fastp*, alignemt using *bwa mem* algorithm, duplicate marking using *Picard* implemented in the *GATK* suite. Allignment statistics are collected using *samtools stats* and *picard*. Coverage information over the target capture region is determined *mosdepth*. Variant calling is performed using *deepvariant* WES model.

![Pipeline](Diagram.png)

## Setup the environment

``` bash
# Install micromamba
echo "installing micromamba" ;
# Linux Intel (x86_64):
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba ;
# Linux/bash:
./bin/micromamba shell init -s bash -p ~/micromamba  # this writes to your .bashrc file
# sourcing the bashrc file incorporates the changes into the running session.
# better yet, restart your terminal!
echo  "export PATH=\${PATH}:$PWD/bin/" >> ~/.bashrc ;
source ~/.bashrc ; 
# micromamba needs to be install
micromamba create -n FastqToDiagnostics \
-c bioconda -c conda-forge awscli bwa samtools mosdepth fastp gatk4 bcftools multiqc
```

## Download refernce genome and indexes

``` bash
micromamba activate FastqToDiagnostics
## sequences
aws s3 --no-sign-request --region eu-west-1 sync \
    s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/ \
    ./references/Homo_sapiens/GATK/GRCh37/
## indexes 
aws s3 --no-sign-request --region eu-west-1 \
    sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/BWAIndex/ \
    ./references/Homo_sapiens/GATK/GRCh37/
## indexes
```

## Download test sequencing rreads : GIAB sample

``` bash
## Intervals
mkdir intervals
wget -c -P ./intervals \
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz
# unzip and strip "chr" string
gunzip -d  ./intervals/nexterarapidcapture_expandedexome_targetedregions.bed.gz
sed -i 's/chr//' ./intervals/nexterarapidcapture_expandedexome_targetedregions.bed
## reads
mkdir reads
wget -c -P ./reads https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget -c -P ./reads https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
## 
```

## Reads Quality control using Fastp

Read quality control and trimming is performed using fastp [@chen2018]

``` bash
mkdir results/
SamplePref=./reads/NIST7035_TAAGGCGA_L001
SampleName=$(basename ${SamplePref})
time fastp -i ${SamplePref}_R1_001.fastq.gz \
      -I ${SamplePref}_R2_001.fastq.gz \
      --detect_adapter_for_pe \
      -L -Q \
      -j results/${SampleName}_fastp.json \
      -h results/${SampleName}_fastp.html \
      -o results/${SampleName}_fastp_R1_001.fastq.gz \
      -O results/${SampleName}_fastp_R2_001.fastq.gz
```

``` bash
# multiqc  for nicer graphs
cd results
multiqc -f .
cd ..
```

## Map the reads to the reference Genome

``` bash
Genome=./references/Homo_sapiens/GATK/GRCh37/human_g1k_v37_decoy.fasta
SamplePref=./reads/NIST7035_TAAGGCGA_L001
SampleName=$(basename ${SamplePref})
R1=results/${SampleName}_fastp_R1_001.fastq.gz
R2=results/${SampleName}_fastp_R1_001.fastq.gz
time bwa mem -t 6 -R  "@RG\tID:${SampleName}\tSM:${SampleName}\tPL:ILLUMINA" \
$Genome  ${R1} ${R2} | samtools view -bS --threads 6 \
| samtools sort --threads 6 > results/${SampleName}.bam  2> results/${SampleName}.bam.log
```

## Mark duplicates

``` bash
# index the alignment
time samtools index --thread 5 results/${SampleName}.bam
# Mark the duplicates with picard
time gatk MarkDuplicatesSpark  \
-I results/${SampleName}.bam \
-O results/${SampleName}.dm.bam \
--remove-all-duplicates false \
-M results/${SampleName}.picard.duplicate.metrics.txt \
-L ./intervals/nexterarapidcapture_expandedexome_targetedregions.bed \
--spark-master local[5]
# index the new bam
time samtools index --thread 5 results/${SampleName}.dm.bam
# multiqc
cd results/ ; multiqc -f . ; cd ..
```

## Collect alignment statistics

``` bash
Genome=./references/Homo_sapiens/GATK/GRCh37/human_g1k_v37_decoy.fasta
# mapping quality 
time samtools stats  \
    results/${SampleName}.dm.bam > results/${SampleName}.samtools.flagstats.txt
# allignment metrics
time gatk CollectAlignmentSummaryMetrics \
          R=${Genome} \
          I=results/${SampleName}.dm.bam  \
          O=results/${SampleName}.picard.allignment.metrics.txt
cd results/ ; multiqc -f . ; cd ..
```

## Collect coverage statistics

``` bash
Interval=./intervals/nexterarapidcapture_expandedexome_targetedregions.bed
# mostdepth to compute target coverage
time mosdepth --by  ${Interval} \
    --threads 5 \
    --thresholds 1,10,15,20,30,50,75,100 \
    results/${SampleName}.dm \
    results/${SampleName}.dm.bam
cd results/ ; multiqc -f . ; cd ..
```

## Variant calling using Deepvariant

This requires docker to be installed on the system.

For illustration purposes we do the calling on a smaller set of the exome

``` bash
# pseudo autosomal regions if do calling on the sex
# chromosome
wget -c -P $PWD/intervals https://storage.googleapis.com/deepvariant/case-study-testdata/GRCh37_PAR.bed
Interval=./intervals/nexterarapidcapture_expandedexome_targetedregions.bed
# run deep variant
BIN_VERSION="1.6.0"
time docker run --rm \
  -v $PWD/results:/input \
  -v $PWD/intervals:/intervals \
  -v $PWD/references:/ref \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref=/ref/Homo_sapiens/GATK/GRCh37/human_g1k_v37_decoy.fasta \
  --reads=/input/${SampleName}.dm.bam \
  --output_vcf=/input/${SampleName}.deepvariant.vcf.gz \
  -num_shards=5 \
  --logging_dir=/input/logs \
  --regions=/intervals/demo_calling.bed
```

## Variant filtration

``` bash
# Extract passed variants
bcftools view -f PASS -Oz -o \
    results/${SampleName}.deepvariant.pass.vcf.gz \
    results/${SampleName}.deepvariant.vcf.gz 
bcftools index -tf results/${SampleName}.deepvariant.pass.vcf.gz
# tabulate some usefull statistics
bcftools stats \
    results/${SampleName}.deepvariant.pass.vcf.gz \
    > results/\${SampleName}.deepvariant.pass.bcftools.stats.txt 
## Multiqc
multiqc cd results/ ; multiqc -f . ; cd ..
```

## Variant annotation and prioritization

Several tools exist for this

The classics

-   VEP : https://grch37.ensembl.org/Multi/Tools/VEP/Ticket?tl=W9TV20C3YLQmIiFE

-   SNPEff

-   Annovar : https://wannovar.wglab.org/done/429351/9u0NoDYBRcjQxGu4/index.html

More integrated tools

OpenCravat

Varsome (commercial)

Enliter (commercial )

Geneius (commercial )

## Validate this pipeline

## References