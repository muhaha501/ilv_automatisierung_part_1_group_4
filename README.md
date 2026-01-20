# ilv_automatisierung_part_1_group_4


# Pig Transcriptome Analysis Pipeline (GSE159583)

This Nextflow pipeline is designed to reproduce the transcriptomic analysis of the study:  
**"Dietary phytonutrients and oxidized vegetable oil improve the immune response and liver metabolism in pigs"** (Le Coz et al., 2021).

## Project Overview

-- TOADD PROJECT OVERVIEW --

## Pipeline Features

- Quality control of raw RNA-Seq data using FastQC
- Dual Processing Approach:
  - Alignment-based quantification using STAR and featureCounts (Use for high Computing Power environments)
  - Alignment-free quantification using Kallisto (If RAM Acess < 16GB)


## Prerequisites

- **Nextflow**
- **Docker** (Derzeit noch nicht umgesetzt da FH Server kein Docker unterstützt)


## Installation & Setup

1. **Clone the Repository**

```bash
git clone [https://git.fh-campuswien.ac.at/c2410542006/ilv_automatisierung_part_1_group_4.git](https://git.fh-campuswien.ac.at/c2410542006/ilv_automatisierung_part_1_group_4.git)
cd ilv_automatisierung_part_1_group_4
```   
2. **Add the Data**

```bash
mkdir data
# Place your FASTQ files in the 'data' directory
# Example: data/Sample1_1.fastq.gz, data/Sample1_2.fastq.gz
```

## Usage

1. STAR Path:

```bash
nextflow run main.nf --tool STAR
```
### Using a Pre-built STAR Index

```bash
nextflow run main.nf --tool STAR --star_index path/to/your/star_index_directory
```

2. Kallisto Path:

```bash
nextflow run main.nf --tool Kallisto
```

## Expected Output

The output files will be organized in the `results` directory, structured as follows:



```
results/
├── fastqc/
│   ├── Sample1_fastqc.html
│   └── Sample1_fastqc.zip
├── counts/
│   ├── Sample1_featurecounts.txt
│   └── ...
└── kallisto_counts/
    ├── Sample1_abundance.tsv
    └── ...
``` 


## ToDo 


- Automate the download of the data from GEO (GSE159583) within the pipeline. (Or additional File)
- Further analysis steps such as differential expression analysis and visualization have to be added in future versions of this pipeline.