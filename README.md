# Reference Sequence and Annotation Files

Supplementary genome reference files in addition to the basic FASTA and GTF.

Disclaimer: Reference genomes are complicated (see: [which human reference genome to use](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)).
Use at your own risk.

---

## Initializing Reference Genome Directory

Genomics bioinformatic tools generally require a reference file/index, but the exact requirements can vary.
This repository includes notes for setting up a reference genome directory that can be used for common analysis types (RNA-seq, ChIP-seq, whole genome/exome sequencing) as well as some pre-generated supplementary files.
The protocol is inspired by Illumina iGenomes, a collection of sequence and annotation files for commonly analyzed genomes.

The two primary files that are required:

* `genome.fa` - genome sequence in FASTA format
* `genes.gtf` - gene annotations in GTF format

These can usually be downloaded from genome browsers like Ensembl, UCSC, or NCBI.
The best source will vary depending on the species of interest.
There are also species-specific resources like FlyBase or WormBase.

To avoid potential downstream problems, be aware of:

* contig names (such as "chr1") should be identical in the FASTA (the full line after `>`) and the GTF (column 1)
* contig names should not have spaces
* the GTF should only contain contigs that are available in the FASTA
* the GTF should include `gene_name` or `gene_id` attributes (in column 9) for all records
* the GTF should include `gene`, `transcript`, and `exon` features (column 3)

## Generating Common Supplementary Files

Generate FASTA index `genome.fa.fai` using samtools:

```bash
samtools faidx genome.fa
```

Generate sequence dictionary `genome.dict` using Picard (used by Picard and GATK):

```bash
java -Xms16G -Xmx16G -jar /path/to/picard.jar CreateSequenceDictionary \
  REFERENCE=genome.fa OUTPUT=genome.dict
```

Generate `chrom.sizes` (used by bedtools and UCSC utilities):

```bash
cut -f 1,2 genome.fa.fai > chrom.sizes
```

Generate genome BED file using bedtools:

```bash
bedtools makewindows -g chrom.sizes -w 999999999 \
  | LC_ALL=C sort -k1,1 \
  > genome.bed
```

Generate Bowtie 2 index (also compatible with Bowtie 1.2.3+):

```bash
mkdir Bowtie2
ln -s ../genome.fa Bowtie2/genome.fa
bowtie2-build --threads 4 Bowtie2/genome.fa Bowtie2/genome
```

Generate BWA index:

```bash
mkdir BWA
ln -s ../genome.fa BWA/genome.fa
bwa index -a bwtsw BWA/genome.fa
```

Generate gene predictions in refFlat format (must be gzipped for Picard CollectRnaSeqMetrics):

```bash
# gtfToGenePred is part of UCSC Genome Browser Utilities
gtfToGenePred -genePredExt -geneNameAsName2 genes.gtf refFlat.tmp
paste <(cut -f 12 refFlat.tmp) <(cut -f 1-10 refFlat.tmp) > refFlat.txt
rm -f refFlat.tmp
gzip refFlat.txt
```

Generate STAR index:

```bash
mkdir STAR
STAR --runMode genomeGenerate --sjdbOverhang 100 --runThreadN 4 \
  --genomeDir STAR --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf
```

Generate kallisto index:

```bash
kallisto index --make-unique \
  --index kallisto.idx \
  transcripts.fa
```

Generate Salmon index:

```bash
salmon index --type quasi --threads 4 \
  --transcripts transcripts.fa \
  --index Salmon
```

Generate  Salmon index for decoy-aware selective alignment (requires Salmon 0.14+):

```bash
bash /path/to/SalmonTools/scripts/generateDecoyTranscriptome.sh \
  -j 4 -g genome.fa -a genes.gtf -t transcripts.fa \
  -o Salmon
salmon index --type quasi --threads 4 \
  --transcripts Salmon/gentrome.fa \
  --decoys Salmon/decoys.txt \
  --index Salmon
```
