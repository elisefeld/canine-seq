# canine-seq

canine-seq is a snakemake workflow for processing bulk paired-end RNA-seq data stored in unaligned CRAM format. The workflow performs quality control, adapter trimming, alignment, gene-level quantification and generation of multiQC reports. Reads that fail to align to the host genome are aligned against a viral genome for detection of viral transcripts. 

## Features

- `samtools`: FASTQ extraction from CRAM files
- `fastqc`: quality control metrics
- `fastp`: Adapter and quality trimming
- `HISAT2`: Genome alignment
- `featureCounts`: Gene-level quantification
- `MultiQC`: Aggregated quality control reports


![Rulegraph of canine-seq workflow](docs/workflow.png "Workflow")
---

## Installation
Clone the repository:

```bash
git clone https://github.com/elisefeld/canine-seq.git
cd project
```

Detailed installation instructions are available in [install.md](docs/INSTALL.md).

### Dependencies
- Snakemake
- Python ≥ 3.14
- pandas
- samtools
- FastQC
- fastp
- HISAT2
- featureCounts (Subread)
- MultiQC

It is recommended to run this program inside a conda environment using a conda package manager and `config/env.yaml`.

---

## Usage

### Input

Raw data and the sample sheet must be stored in a project-specific directory under `data/`. Reference files for both host and virus must be stored under a reference directory  specified by the name variable `config.yaml`.

The sample sheet must contain a column named `sample_id` containing sample names that exactly match the CRAM filenames (excluding extensions). An example sample sheet is provided [here](config/example_sheet.csv).

```text
canine_seq/
├── config/
|   ├── env.yaml
│   └── config.yaml
├── data/
│   ├── <PROJECT>/
│   │   ├── sample.cram
│   │   ├── sample.cram.crai
│   │   └── <SAMPLE_SHEET>
│   │
│   └── reference/
│       └── <REF_MAIN_NAME>/
│           ├── <FASTA>
│           └── <GTF>
|        └── <REF_VIRUS_NAME>/
|           ├── <FASTA>
│           └── <GTF>
└── runs/
    └── <RUN_NAME>/
        └── logs/
          ├── main/
          └── virus/
        └── report/
          ├── main/
          └── virus/
        └── igv/
          ├── main/
          └── virus/
        └── results/
          ├── main/
          └── virus/
```

### Example Configuration

**config.yaml**

```yaml
project:
  name: my_project

run:
  run_name: test_run
  sample_sheet: sample_sheet.csv

ref:
  main:
    name: host_reference
    fasta: genome.fa
    gtf: genes.gtf

  virus:
    name: virus_reference
    fasta: virus.fa
    gtf: virus.gtf

threads:
  samtools: 8
  fastqc: 2
  fastp: 8
  hisat2: 16
  featurecounts: 4
```

**sample_sheet.csv**

```csv
sample_id
SAMPLE001
```

### Example Directory Layout

```text
canine_seq/
├── config/
|   ├── env.yaml
│   └── config.yaml
├── data/
│   ├── my_project/
│   │   ├── SAMPLE001.cram
│   │   ├── SAMPLE001.cram.crai
│   │   └── sample_sheet.csv
│   │
│   └── reference/
│       └── host_reference/
│           ├── genome.fa
│           └── genes.gtf
|        └── virus_reference/
|           ├── virus.fa
│           └── virus.gtf
└── runs/
    └── test_run/
        └── logs/
          ├── main/
          └── virus/
        └── report/
          ├── main/
          └── virus/
        └── igv/
          ├── main/
          └── virus/
        └── results/
          ├── main/
          └── virus/
```

### Threads
The number of threads can be specified for each tool using `config.yaml`.
```yaml
threads:
  samtools: 6
  fastqc: 1
  fastp: 4
  kallisto: 4
  hisat2: 8
  featurecounts: 2
```
---

## Running the Workflow

Perform a dry run:

```bash
snakemake -n -p
```

Generate a workflow DAG:

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

Generate a workflow rulegraph:

```bash
snakemake --rulegraph | dot -Tpdf > rulegraph.pdf
```

Run the pipeline:

```bash
snakemake --cores 8
```



## Output

Upon completion, the workflow will generate:

```text
runs/
└── <RUN_NAME>/
    ├── logs/
    │   ├── main/
    │   └── virus/
    │
    ├── report/
    │   ├── main/
    │   │   └── multiqc_report.html
    │   │
    │   └── virus/
    │       └── multiqc_virus_report.html
    │
    ├── results/
    │   ├── main/
    │   │   ├── fastp/
    │   │   ├── fastqc/
    │   │   ├── featurecounts/
    │   │   ├── flagstat/
    │   │   ├── hisat2/
    │   │   └── sort/
    │   │
    │   └── virus/
    │       ├── fastqc/
    │       ├── featurecounts/
    │       ├── fastq/
    │       ├── hisat2/
    │       ├── sort/
    │       └── unmapped/
    │
    └── igv/
        └── main/
            ├── *.bam
            └── *.bam.bai
```

`featurecounts` will contain the unmerged gene counts for each sample. 

---
## License

This project is distributed under the MIT license. See the [license](LICENSE) for more details. 

## Citations

If you use this workflow in a publication, please cite the following tools:

- Snakemake
- samtools
- FastQC
- fastp
- HISAT2
- featureCounts (Subread)
- MultiQC

Please also cite the reference genome and gene annotation resources used in your analysis.