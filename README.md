# canine-seq

This project automates the analysis of bulk RNA-seq data (.cram format) from canine subjects. 

## Overview

This project provides a reproducible workflow for:

- Converting .cram to .bam
- 
- 

The workflow is implemented using Snakemake and is designed to support reproducible analysis.

## Features

- Reproducible workflow execution
- Automated dependency management
- Modular pipeline design
- Configurable parameters
- Built-in quality control and reporting

## Workflow

```text
Input Data
    ↓
Preprocessing
    ↓
Analysis
    ↓
Quality Control
    ↓
Results
```

## Quick Start

### Dependencies
- snakemake v
- a conda package manager (miniforge3 was used for this project)
 
Clone the repository:

```bash
git clone https://github.com/elisefeld/canine-seq.git
cd project
```

Install dependencies:

```bash
# See INSTALL.md for detailed setup
```

Run the workflow:

```bash
snakemake --cores 8
```

## Repository Structure

```text
.
├── config/         # configuration files
├── docs/           # project documentation
├── logs/           # workflow logs
├── resources/      # raw and reference data
├── scripts/        # analysis scripts
├── workflow/       # snakefile and rules
└── results/        # generated outputs
```

## Configuration

Pipeline settings are controlled through:

```text
config/config.yaml
```

Modify the configuration file before running analyses.

## Input Data

This project currently only supports already aligned data in cram format. Add a folder in /resources/data/raw named with the project name. It should include the raw data and a sample sheet.

```text
/resources/data/raw/PROJECT_NAME
├── sample_sheet.csv
├── file_1.cram
└── file_1.cram.crai
```

## Outputs

The workflow generates:

| Output | Description |
|----------|-------------|
| results/summary.csv | Final summary table |
| results/figures/ | Generated plots |
| results/reports/ | QC reports |


## Installation

Detailed installation instructions are available in
[INSTALL.md](INSTALL.md).

## License

This project is distributed under the MIT license. See license.md for more details. 



