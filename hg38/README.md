
`genome.fa`

```bash
wget -O genome.fa.gz \
  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip genome.fa.gz
```

`genes.gtf`

```bash
wget -O genes.gtf.gz \
  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz
gunzip genes.gtf.gz
```
