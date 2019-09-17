
`genome.fa`

```bash
wget -O genome.fa.gz \
  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip genome.fa.gz
```

`genome.fa.fai`

```bash
samtools faidx genome.fa
```

`genome.dict`

```bash
java -Xms8G -Xmx8G -jar ${PICARD_ROOT}/picard.jar CreateSequenceDictionary REFERENCE=genome.fa OUTPUT=genome.dict
```

`chrom.sizes`

```bash
cut -f 1,2 genome.fa.fai > chrom.sizes
```

`genome.bed`

```bash
bedtools makewindows -g chrom.sizes -w 999999999 | LC_ALL=C sort -k1,1 > genome.bed
```

`genes.gtf`

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz
```
