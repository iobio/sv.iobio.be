# Data Information

Data for bands and chromosomes was retrieved from RefSeq. Their GRCh37 & 38 gff.gz files were obtained from: [NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml). These files were filtered in the command line for particular features (specifics below).

## Filtering

### Chromosomes

```shell
gunzip -c 'Human Genome RefSeq.gff.gz' | awk '$3 == "region" && $1 ~ /^NC/ && $2 == "RefSeq" {print}' > chromosomes.txt
```

**once parsed there is an additional entry that doesn't match the 23 that I've removed*

### Bands

Bands retrieved from UCSC data repositories as text files. [Docs](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=cytoBand&hgta_table=cytoBand&hgta_doSchema=describe+table+schema)

### Centromeres

Centromeres were also retrieved from UCSC data repositories. [Docs](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=centromeres&hgta_table=centromeres&hgta_doSchema=describe+table+schema)

hg19 did not provide a centromeres only reference but was listed in the "gaps" reference [Docs](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema)
