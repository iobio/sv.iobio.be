# Getting the data

You'll need a few data sources before building the Docker image. Everything is
stored in `./data`, combining with some files that are already kept in the git
repo.


## Dose Sensitivity

If these links die, try going to this page: https://search.clinicalgenome.org/kb/downloads#section_dosage

```bash
curl https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv -o data/ClinGen_gene_hg19.tsv
curl https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv -o data/ClinGen_gene_hg38.tsv
curl https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv -o data/ClinGen_region_hg19.tsv
curl https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv -o data/ClinGen_region_hg38.tsv
```


## geneinfo sqlite DB

See https://github.com/iobio/geneinfo.db for building instructions. Copy the
file to `data/geneinfo.db/gene.iobio.db`.


## HPO sqlite DB

See https://github.com/iobio/pheno_matcher_be_rust/blob/main/src/hpoAssociations/populateDB.py
and the containing folder for build instructions. Copy the file to
`data/hpo.db`.


## SVAFotate files

Currently you need to get these from Tom in the Quinlan lab:

```
data/SVAFotate_popAFs_GRCh37.sorted.v4.1.bed.gz
data/SVAFotate_popAFs_GRCh37.sorted.v4.1.bed.gz.tbi
data/SVAFotate_popAFs_GRCh38.sorted.v4.1.bed.gz
data/SVAFotate_popAFs_GRCh38.sorted.v4.1.bed.gz.tbi
```

Your final data directory should look like this:

```
data
├── centromeres_hg38.txt
├── chromosomes_hg19.txt
├── chromosomes_hg38.txt
├── ClinGen_gene_hg19.tsv
├── ClinGen_gene_hg38.tsv
├── ClinGen_region_hg19.tsv
├── ClinGen_region_hg38.tsv
├── cytoBand_hg19.txt.gz
├── cytoBand_hg38.txt.gz
├── data_info.md
├── gaps_ref_hg19.txt.gz
├── geneinfo.db
│   └── gene.iobio.db
├── hpo.db
├── SVAFotate_popAFs_GRCh37.sorted.v4.1.bed.gz
├── SVAFotate_popAFs_GRCh37.sorted.v4.1.bed.gz.tbi
├── SVAFotate_popAFs_GRCh38.sorted.v4.1.bed.gz
└── SVAFotate_popAFs_GRCh38.sorted.v4.1.bed.gz.tbi
```


# Building the image

Additional details in [deploy_notes.md].

You'll need Docker and Apptainer (forked from Singularity. CHPC still uses
Singularity, but commands are the same).

```bash
docker build -t backend.sv.iobio .
```

```bash
apptainer build backend.sv.iobio.sif docker-daemon://backend.sv.iobio
```

# Running

To run in CHPC:

```bash
module load singularity
singularity exec backend.sv.iobio.sif node /app.js
```

You can optionally provide a port:

```bash
singularity exec backend.sv.iobio.sif node /app.js 7478
```
