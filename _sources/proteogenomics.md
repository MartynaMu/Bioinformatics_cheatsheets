Proteogenomics for neoantigen discovery
================
Martyna Siewiera
2025-09-26

## Overview

## Computational set-up

Windows 10 with wsl2, docker, nextflow, miniconda

## Data input:

• .fastq files from RNA sequencing of paired tumor and nascent adjacent
tissue or PBMC

• Reference genome FASTA from ensemble, 113 release, soft-masked,
primary assembly

## Tools and environments:

• GATK4

• Samtools …

• Bcftools …

• RNA-seq pipeline …

## Notes:

First process the truncated .fastq files to ensure the pipeline is
working, use seqtk, for example:

```bash
gzip -d <1.FASTQ.GZ>
gzip -d <2.FASTQ.GZ>

seqtk -s100 <1.FASTQ> 100000 > sub1.FASTQ
seqtk -s100 <2.FASTQ> 100000 > sub2.FASTQ

gzip *.FASTQ
```

## Step by step:

-----

### 1. FASTQC and MULTIQC

It's better to evaluate the sequencing quality before running the whole
pipeline. Then you can skip those later.

```bash
conda activate MOOC

mkdir -p QC_Reports

fastqc *_1.fastq *_2.fastq --ourdir QC_Reports

multiqc .
```

-----

### 2. Sample sheet preparation

A list of fastq files and their association is required for the RNASEQ
pipeline. Use 'create-samplesheet.R' script, that goes like this:

```r
library(tidyverse)
path_to_project <- "D:/NSCLC_cohort" #directory of the project
setwd(path_to_project)

files <- list.files("data", full.names = TRUE) #directory of where the fastq files are
samples <- files %>%
  str_sub(end = -12, start = 6) %>%
  unique()

samplesheet <- data.frame(sample = samples,
                          fastq_1 = keep(files, str_detect(files, "_1")),
                          fastq_2 = keep(files, str_detect(files, "_2")),
                          strandedness = "auto") # or reverse if confirmed

write_csv(samplesheet, "samplesheet.csv")
```

-----

### 3. RNA-seq Nextflow pipeline with optional reference mapping with STAR

Open docker app, in terminal go to the project directory and run:

```bash
sudo nextflow run \
  nf-core/rnaseq \
  --input samplesheet.csv \
  --fasta /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  --gtf /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.gtf \
  --star_index /mnt/d/Reference/STAR_reference/index/star \
  --outdir results \
  --email mamuszcz@gmail.com \
  --aligner star_salmon \
  --stringtie_ignore_gtf true \ #will use new transcripts resulted from stringtie
  --save_align_intermeds true \
  --skip_fastqc true \
  --extra_star_align_args '--outTmpDir ~/star_temp' \
  -profile docker \
  -c star.config

```

star.config file looks as follows:

```json
process {
  resourceLimits = [
    cpus: 7,  # leave some for other processes
    memory: '27.GB', # with wsl2 set to 30 out of 32 GB
    time: 24.h
  ]
  withName: STAR_GENOMEGENERATE { # in case you wanna index with STAR
    ext.args = [
      '--genomeChrBinNbits 15',
      '--genomeSAsparseD 2'
    ].join(' ').trim()
  }
}
```

#### Evaluate new transcripts:

```bash
gffcompare -G -r /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.113.gtf results/star_salmon/stringtie/SAMPLE.transcripts.gtf
```

You may want to add new transcripts to the fresh run and calculate their
abundance

-----

### 4. Create fasta dictionary

```bash
gatk CreateSequenceDictionary -R /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
```

-----

### 5. Variant calling

With RNA-seq data, we must splt N cigars which splits reads by splice
junctions that Mutect2 isn't aware of

```bash
gatk SplitNCigarReads \
  -R /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  -I sampleA.markdup.sorted.bam \
  -O sampleA.markdup.sorted.split.bam \
  --create-output-bam-index true
```

#### Somatic variants call - Mutect2

```bash
gatk Mutect2 \
  -R /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  -I sampleB.markdup.sorted.split.bam \
  -I sampleA.markdup.sorted.split.bam \
  -normal sampleA \
  --dont-use-soft-clipped-bases true \
  --output sample_somatic_variants.vcf \
  --min-base-quality-score 20 # optional
  # --genotype-germline-sites true
```

Now filter Mutect2 calls

```bash
gatk FilterMutectCalls \
  -R /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  -V sample_somatic_variants.vcf \
  -O sample_somatic_variants_FILT.vcf
```

#### Germline variant calling

Call tumor and benign tissue one by one:

```bash
gatk HaplotypeCaller \
  -Xmx26g \
  -R /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  -I sampleA.markdup.sorted.split.bam \
  -O sampleA.g.vcf.gz \
  -ERC GVCF \ # generates a GVCF for later joint genotyping for tumor and healthy tissue 
  --dont-use-soft-clipped-bases # i.a. dont use introns (?)
  
```

Combine GVCF files for common SNV \<- germline

```bash
gatk CombineGVCFs \
  -R /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  --variant sampleA.g.vcf.gz \
  --variant sampleB.g.vcf.gz \
  -O sample_combined.g.vcf.gz
```

Genotype the variants:

```bash
gatk GenotypeGVCFs \
  -R /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  -V sample_combined.g.vcf.gz \
  -O sample_germline_variants.vcf
```

#### Variant annotation and filtering

How many variants pass the filter?

```bash
grep -o 'SITE' sample_somatic_variants_FILT.vcf | wc -l
```

You can save only the records that PASS the filter:

```bash
bcftools view -f PASS sample_somatic_variants_FILT.vcf -o sample_somatic_variants_FILT2.vcf 
```

Annotate with SnpEff or ... to estimate the effects of the variant.

```bash
snpEff -Xmx26g GRCh38.p14 -cvsStats sample_vc_summary.csv -v sample_somatic_variants_FILT.vcf > sample_somatic_variants_ANN.vcf
# or without csv summary
snpEff -Xmx26g GRCh38.p14 -v sample_somatic_variants_FILT.vcf > sample_somatic_variants_ANN.vcf
```

```bash
SnpSift -Xmx26g annotate \
  -v /mnt/d/Reference/ClinVar/clinvar.vcf.gz \
  sample_somatic_variants_ANN.vcf > sample_somatic_variants_ANN_CLINVAR.vcf
```

-----

### 6. Filter for non-synonymous variants (i.a. non-silent, resulting in a change in AA sequence)

```bash
bcftools view \
  -i "INFO/ANN~'missense-variants\|stop_gained\|stop_lost\|start_lost'" sample_somatic_variants_ANN.vcf \
  -o sample_somatic_variants_ANN_NS.vcf
```

-----

### 7. Find mutations which were found as germline AND somatic

Prepare for intersection:

```bash
bgzip -c sampleA.vcf > sampleA.vcf.gz
tabix -p vcf sampleA.vcf.gz
```

Find and save shared variants between germline and somatic:

```bash
bcftools isec \
  -n=2 \
  -O v \
  -o sample_shared_variants.vcf \
  sample_somatic_variants.vcf.gz \
  sample_germline_variants.vcf.gz 
```

Identify and save somatic only:

```bash
bcftools isec \
  -n=1 \
  -c none \
  -O v \
  -o sample_somatic_only_variants.vcf \
  sample_somatic_variants.vcf.gz \
  sample_germline_variants.vcf.gz
```

Identify and save germline only:

```bash
bcftools isec \
  -n=1 \
  -c none \
  -O v \
  -o sample_germline_only_variants.vcf \
  sample_germline_variants.vcf.gz \
  sample_somatic_variants.vcf.gz
```

#### Reasonable read depth for variants:

-   somatic: 100 dp
-   germline polymorphism: 30 dp

-----

### 8. Incorporate variants into nucleotide FASTA reference

```bash
gatk FastaAlternateReferenceMaker \
   -R /mnt/d/Reference/113-Release/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
   -O output.fasta \
   -V input.vcf \
```

-----

### 9. Translate nucleotide FASTA with variants to FASTA AA sequence with 6 frames (three forward and three reverse frames)

Acc. to perplexity, use either Online tools:

• <https://www.ebi.ac.uk/jdispatcher/st/emboss_transeq>

• <https://www.bioinformatics.org/sms2/translate.html>

Offline:

• biostring R package `translate()`

• Biopython Bio.Seq (Seq) and Bio.Alphabet (IUPAC) modules -
`translate()`

