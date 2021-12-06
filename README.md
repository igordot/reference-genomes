# Reference Sequence and Annotation Files

Supplementary genome reference files in addition to the basic FASTA and GTF as well as notes and commentary.

Disclaimer: Reference genomes are complicated and there are many caveats, even for commonly used ones (see: [which human reference genome to use](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)).
Use at your own risk.

---

## Initializing Reference Genome Directory

Genomics bioinformatic tools generally require a reference file/index, but the exact requirements can vary.
This repository includes notes for setting up a reference genome directory that can be used for common analysis types (RNA-seq, ChIP-seq, whole genome/exome sequencing) as well as some pre-generated supplementary files.
The organization scheme is inspired by Illumina iGenomes, a collection of sequence and annotation files for commonly analyzed genomes.
The files are placed in separate directories based on the genome reference version, such as `hg38` or `mm10`.
Within each genome directory, the files are named based on the type.

The two primary files that are required:

* `genome.fa` - genome sequence in FASTA format
* `genes.gtf` - gene annotations in GTF format

These can usually be downloaded from genome browsers like Ensembl, UCSC, or NCBI.
The best source will vary depending on the species of interest.
There are also species-specific resources like FlyBase or WormBase.
The differences can be minor (different contig names), but they could be more substantial (thousands of additional genes). 

## Sanity Testing

There are a few issues to be aware of to avoid potential downstream problems.
Both `genome.fa` and `genes.gtf` have to be valid FASTA and GTF files, but they should also be compatible with each other.

If a FASTA index and a sequence dictionary can be successfully generated (see below), the FASTA is likely compatible with all tools.
One exception is that contig names (such as "chr1") with spaces will cause issues for some applications.
Check contig names with:

```bash
grep "^>" genome.fa
```

Contig names should be identical in the FASTA and the GTF.
Create a list of contigs in both files and check that the GTF does not have more contigs than the FASTA:

```bash
cat genome.fa.fai | cut -f 1 | sort > contigs.fa.txt
cat genes.gtf | cut -f 1 | sort | uniq > contigs.gtf.txt
wc -l contigs.*.txt
```

This counts the total number of contigs from each file.

Check that all the contig names in the GTF are present in the FASTA:

```bash
diff --side-by-side --suppress-common-lines contigs.fa.txt contigs.gtf.txt
```

Only the non-matching contigs will be output.
If the FASTA has additional contigs, that should be fine.

The GTF should include `gene_name` and `gene_id` attributes for all records.
Check the number of genes with and without those attributes:

```bash
cat genes.gtf | grep -F -v "transcript_id" | grep -F -w "gene" | wc -l
cat genes.gtf | grep -F -v "transcript_id" | grep -F -w "gene" | grep -F "gene_name" | wc -l
cat genes.gtf | grep -F -v "transcript_id" | grep -F -w "gene" | grep -F "gene_id" | wc -l
```

The three numbers should be identical.

The GTF should include `gene`, `transcript`, and `exon` features (column 3) for all genes.

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
kallisto index --make-unique --index kallisto.idx transcripts.fa  
```

Generate Salmon index:

```bash
salmon index --type quasi --threads 4 \
  --transcripts transcripts.fa \
  --index Salmon
```

Generate Salmon (version 1.2+) index for decoy-aware selective alignment:

```bash
mkdir Salmon
cat transcripts.fa genome.fa > Salmon/gentrome.fa
gzip -v Salmon/gentrome.fa
# a default k of 31 is optimized for reads >75bp, consider a smaller k with shorter reads
salmon index --threads 4 --kmerLen 23 \
  --transcripts Salmon/gentrome.fa.gz \
  --decoys Salmon/decoys.txt \
  --gencode \
  --index Salmon
```
