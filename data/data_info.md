# Data Information

Data for bands and chromosomes was retrieved from RefSeq. Their GRCh37 & 38 gff.gz files were obtained from: [NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml). These files were filtered in the command line for particular features (specifics below).

## Filtering

### Chromosomes

```shell
gunzip -c 'Human Genome RefSeq.gff.gz' | awk '$3 == "region" && $1 ~ /^NC/ && $2 == "RefSeq" {print}' > chromosomes.txt
```

**once parsed there is an additional entry that doesn't match the 23 that I've removed*

### Bands

```shell
gunzip -c 'Human Genome RefSeq.gff.gz' | awk '$3 == "region" && $1 ~ /^NW/ && $2 == "RefSeq" {print}' > bands.txt
```
