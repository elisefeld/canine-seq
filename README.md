# canine-seq

This project automates the analysis of paired-end bulk RNA-seq data in unaligned CRAM format.  

## Features

- `samtools`: FASTQ extraction from CRAM files
- `fastqc`: quality control metrics
- `fastp`: Adapter and quality trimming
- `HISAT2`: Genome alignment
- `featureCounts`: Gene-level quantification
- `kallisto` Pseudoalignment and gene-level quantification. Disabled by default.
- `MultiQC`: Aggregated quality control reports

---
## Usage

### Input Structure

Raw data and the sample sheet must be stored in a project-specific directory under `data/`. Reference files must be stored under a reference directory named using the species, build, and release values specified in `config.yaml`.

The sample sheet must contain a column named `sample_id` containing sample names that exactly match the CRAM filenames (excluding extensions).

```text
canine_seq/
├── config/
│   └── config.yaml
├── data/
│   ├── <PROJECT>/
│   │   ├── sample.cram
│   │   ├── sample.cram.crai
│   │   └── <SAMPLE_SHEET>
│   │
│   └── reference/
│       └── <SPECIES>_<BUILD>_<RELEASE>/
│           ├── <FASTA>
│           ├── <CDNA>
│           └── <GTF>
│
└── runs/
    └── <RUN_NAME>/
        ├── logs/
        ├── report/
        └── results/
```

### Example Configuration

**config.yaml**

```yaml
project:
  name: my_project

run:
  run_name: my_run
  sample_sheet: sample_sheet.csv
  ref_species: Canis_lupus_familiarisgsd
  ref_build: UU_Cfam_GSD_1.0
  ref_release: 116

ref:
  fasta: genome.fa
  cdna: cdna.fa
  gtf: genes.gtf.gz
```

**sample_sheet.csv**

```csv
sample_id
SAMPLE001
```

### Example Directory Layout

```text
canine_seq/
├── data/
│   ├── my_project/
│   │   ├── SAMPLE001.cram
│   │   ├── SAMPLE001.cram.crai
│   │   └── sample_sheet.csv
│   │
│   └── reference/
│       └── Canis_lupus_familiarisgsd_UU_Cfam_GSD_1.0_116/
│           ├── genome.fa
│           ├── cdna.fa
│           └── genes.gtf.gz
│
└── runs/
    └── my_run/
        ├── logs/
        ├── report/
        └── results/
```

---

## Dependencies
- Snakemake
- Python ≥ 3.9
- pandas
- samtools
- FastQC
- fastp
- HISAT2
- kallisto
- featureCounts (Subread)
- MultiQC

It is recommended to run this program inside a conda enviornment using a conda package manager and `config/env.yaml` .


## Installation

Clone the repository:

```bash
git clone https://github.com/elisefeld/canine-seq.git
cd project
```

Detailed installation instructions are available in:

```text
INSTALL.md
```
---
## Running the Workflow

Perform a dry run:

```bash
snakemake -n -p
```

Run the pipeline:

```bash
snakemake --cores 8
```

Generate a workflow DAG:

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

---

## Output

Upon completion, the workflow will generate:

```text
runs/<RUN_NAME>/
├── logs/
├── report/
│   └── multiqc_report.html
└── results/
    ├── fastqc/
    ├── fastp/
    ├── flagstat/
    ├── kallisto/
    ├── featurecounts/
    ├── bam/
    └── hisat2/
```

## License

This project is distributed under the MIT license. See license.md for more details. 

## Citations

If you use this workflow in a publication, please cite the following tools:

- Snakemake
- samtools
- FastQC
- fastp
- HISAT2
- kallisto
- featureCounts (Subread)
- MultiQC

Please also cite the reference genome and gene annotation resources used in your analysis.