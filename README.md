# ilv_automatisierung_part_1_group_4


# Pig Transcriptome Analysis Pipeline (GSE159583)

This Nextflow pipeline is designed to reproduce the transcriptomic analysis of the study:  
**"Exploring With Transcriptomic Approaches the Underlying Mechanisms of an Essential Oil-Based Phytogenic in the Small Intestine and Liver of Pigs"** (Le Coz et al., 2021).

## Project Overview

-- TOADD PROJECT OVERVIEW --

## Pipeline Features

- Flexible input options: Download data from SRA or use local FASTQ files
- Quality control of raw RNA-Seq data using FastQC
- Dual Processing Approach:
  - Alignment-based quantification using STAR and featureCounts (Use for high Computing Power environments)
  - Alignment-free quantification using Kallisto (Recommended If RAM Acess < 16GB)
- Post-processing analysis including PCA, correlation heatmaps, and differential expression analysis using DESeq2

## Prerequisites

### Usage of Docker 

- **Nextflow**
- **Docker** 
- **Internet Connection** 

### Without Docker

If Docker is not used the following software must be installed:

- **SRA Toolkit**
- **FastQC**
- **STAR**
- **Kallisto**
- **Subread (featureCounts)**

Additionally the following python packages are required for the post-processing analysis:

- **pandas**
- **numpy**
- **matplotlib**
- **seaborn**
- **scikit-learn**
- **scipy**

## Installation & Setup

1. **Clone the Repository**

```bash
git clone https://git.fh-campuswien.ac.at/c2410542006/ilv_automatisierung_part_1_group_4.git
cd ilv_automatisierung_part_1_group_4
```   
2. **Add the Data (if local data is used)**

```bash
mkdir data
# Place your FASTQ files in the 'data' directory
# Example: data/Sample1_1.fastq.gz, data/Sample1_2.fastq.gz
```

## Basic Usage

### Running the Pipeline

```bash
nextflow run main.nf 
```

## Advanced Usage

### Downloading Data from SRA (if needed)

```bash
nextflow run main.nf --download_sra true
```

### Choosing the Quantification Tool

1. STAR Path:

```bash
nextflow run main.nf --tool STAR
```
#### Using a Pre-built STAR Index

```bash
nextflow run main.nf --tool STAR --star_index path/to/your/star_index_directory
```

2. Kallisto Path:

```bash
nextflow run main.nf --tool KALLISTO
```

## Parameters

### General Parameters

- `--tool`: Choose between `STAR` or `KALLISTO` for the quantification method.
- `--reads`: Path to the directory containing FASTQ files (default: `data/`).
- `--output_dir`: Path to the output directory (default: `results/`).
- `--docker`: Use Docker containers for the pipeline (default: true).

### Reference Data

- `--genome_url`: URL to download the reference genome FASTA file.
- `--gtf_url`: URL to download the GTF annotation file.
- `--transcriptome_url`: URL to download the transcriptome FASTA file.

### STAR Specific Parameters

- `--star_index`: Path to a pre-built STAR index directory (if not provided, the index will be built).
- `--max_forks`: Maximum number of process instances to be executed in parallel for STAR alignment (default: 5).


## Expected Output

The output files will be organized in the `results` directory, structured as follows:

```
results/
├── alignment/
│   ├── Sample1_Aligned.out.bam
│   └── ...
├── analysis/
│   ├── deseq2/
│   │   ├── de_results/
│   │   └── figures/
│   │   └── ...
│   └── qc/
│       ├── figures/
│       └── ...
│       
├── counts/
│   ├── kallisto/
│   │   ├── Sample1/
│   │   │   └── abundance.tsv
│   │   └── ...
│   └── STAR/
│       ├── Sample1_featurecounts.txt
│       └── ...
├── fastqc/
│   ├── Sample1_fastqc.html
│   ├── Sample1_fastqc.zip
│   └── ...
└── index/
    ├── kallisto_index/
    └── star_index/
``` 

| Verzeichnis | Inhalt |
| :--- | :--- |
| `fastqc/` | Qualitätsberichte der Rohdaten (HTML) |
| `alignment/` | BAM-Dateien (only if STAR is used) |
| `counts/` | Gen-Counts (STAR) or abundance tables (Kallisto) |
| `analysis/qc/` | PCA-plots, correlation heatmaps and count matrices |
| `analysis/deseq2/` | Volcano-plots and DE-tables |

## Software & Versions

| Tool | Docker Image / Tag | Version | 
| :--- | :--- | :--- |
| **SRA Tools** | `sra-tools:3.0.3--h87f3376_0` | 3.0.3 |
| **FastQC** | `fastqc:0.11.9--0` | 0.11.9 | 
| **STAR** | `star:2.6.1d--0` | 2.6.1d | 
| **Kallisto** | `kallisto:0.46.2--h4f7b962_1` | 0.46.2 | 
| **Subread** | `subread:2.1.1--h577a1d6_0` | 2.1.1 | 
| **Scanpy / Python** | `scanpy:1.4.4--py_0` | 1.4.4 |