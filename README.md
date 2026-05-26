# NGS Amplicon Analysis Pipeline for CRISPR Genome Editing Assessment

Bioinformatics pipeline for analyzing next-generation sequencing (NGS) amplicon data from CRISPR-Cas genome editing experiments in plants, using the **Hi-TOM** (High-Throughput Tracking of Mutations) barcoding system.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Innocent250/ngs-analysis-pipeline/blob/main/NGS_Analysis_Pipeline.ipynb)

## Overview

This notebook implements an end-to-end analysis workflow for assessing CRISPR-mediated genome editing outcomes using Illumina amplicon sequencing data.

## What This Pipeline Does

- Processes Illumina MiSeq paired-end amplicon sequencing reads
- Merges paired-end reads with using overlap sequences.
- Demultiplexes merged reads by Hi-TOM barcode pairs (up to 96 samples)
- Quantifies editing outcomes (indel frequencies, allele distributions) using CRISPResso2
- Supports multiplexed target analysis from barcoded amplicon libraries

## Pipeline Steps

```
                      WET LAB                                    COMPUTATIONAL
 ┌─────────────────────────────────────┐      ┌──────────────────────────────────────────┐
 │  Genomic DNA                        │      │  FASTQ R1 + R2                           │
 │       │                             │      │       │                                  │
 │       ▼                             │      │       ▼                                  │
 │  PCR Round 1 (target-specific       │      │  Read Merging (fastp)                    │
 │   + bridging sequences)             │      │       │                                  │
 │       │                             │      │       ▼                                  │
 │       ▼                             │      │  Demultiplexing (barcode splitting)      │
 │  PCR Round 2 (Hi-TOM barcoding      │      │       │                                  │
 │   + Illumina adapters)              │      │       ▼                                  │
 │       │                             │      │  CRISPResso2 Batch Analysis              │
 │       ▼                             │      │       │                                  │
 │  Pool + cleanup + sequencing        │      │       ▼                                  │
 │  (Illumina MiSeq 2x250bp)           │      │  Indel frequencies, allele tables,       │
 │                                     │      │  mutation profiles per sample            │
 └─────────────────────────────────────┘      └──────────────────────────────────────────┘
```

## Hi-TOM Library Preparation

This pipeline is designed for data generated using the **Hi-TOM** barcoding system ([Liu et al., 2019](https://doi.org/10.1007/s11427-019-9570-x)). Below is a summary of the primer design and barcoding strategy.

### Primer Architecture

The Hi-TOM system uses two rounds of PCR:

**Round 1 — Target-specific amplification with bridging sequences:**

```
Forward: 5'- ggagtgagtacggtgtgc -[target-specific F primer]- 3'
Reverse: 5'- gagttggatgctggatgg -[target-specific R primer]- 3'
                 └── bridging sequence (same for all targets)
```

**Round 2 — Barcoding + Illumina adapter addition:**

```
Forward: 5'- [Illumina adapter] -[F barcode]- ggagtgagtacggtgtgc - 3'
Reverse: 5'- [Illumina adapter] -[R barcode]- gagttggatgctggatgg - 3'
```

### Illumina Adapter Sequences example

| Read | Adapter Sequence (5' → 3') |
|------|----------------------------|
| Forward | `ACACTCTTTCCCTACACGACGCTCTTCCGATCT` |
| Reverse | `GACTGGAGTTCAGACGTGTGCTCTTCCGATCT` |

### Barcode Table (96 unique combinations)

Forward barcodes (columns 1-12 on a 96-well plate):

| ID | Barcode |
|----|---------|
| Hi-TOM-F-1 | `gcttGCGTt` |
| Hi-TOM-F-2 | `gcttGTAGt` |
| Hi-TOM-F-3 | `gcttACGCt` |
| Hi-TOM-F-4 | `gcttCTCGt` |
| Hi-TOM-F-5 | `gcttGCTCt` |
| Hi-TOM-F-6 | `gcttAGTCt` |
| Hi-TOM-F-7 | `gcttCGACt` |
| Hi-TOM-F-8 | `gcttGATGt` |
| Hi-TOM-F-9 | `gcttATACt` |
| Hi-TOM-F-10| `gcttCACAt` |
| Hi-TOM-F-11| `gcttGTGCt` |
| Hi-TOM-F-12| `gcttACTAt` |

Reverse barcodes (rows A-H on a 96-well plate):

| ID | Barcode |
|----|---------|
| Hi-TOM-R-A | `ctgtGCGTt` |
| Hi-TOM-R-B | `ctgtGTAGt` |
| Hi-TOM-R-C | `ctgtACGCt` |
| Hi-TOM-R-D | `ctgtCTCGt` |
| Hi-TOM-R-E | `ctgtGCTCt` |
| Hi-TOM-R-F | `ctgtAGTCt` |
| Hi-TOM-R-G | `ctgtCGACt` |
| Hi-TOM-R-H | `ctgtGATGt` |

Forward barcodes are preceded by `gctt` and reverse by `ctgt` to enhance specificity. The 12 forward x 8 reverse combinations give 96 unique sample barcodes.

### Amplicon Size Considerations

- Azenta Amplicon-EZ sequences up to **500 bp** (including bridging and barcode sequences, excluding Illumina adapters) with 2x250 bp reads
- Target amplicons of **250-400 bp** for optimal read merging (maximizes overlapping sequence)
- The overlapping region in the middle is used for quality-aware merging of forward and reverse reads

### Pooling Guidelines

| Sample type | Recommended pool size | Rationale |
|------------|----------------------|-----------|
| Plant tissue (T0 lines) | Up to 48 samples | Higher coverage per sample |
| Protoplast assays | Up to 24 samples | Even higher coverage needed |
| Maximum capacity | 96 samples | Limited by unique barcode combinations |

Reducing the number of samples per pool increases sequencing depth per sample.

### Ready-to-Use Barcode Reference

A complete 96-well plate barcode CSV template is included at [`barcodes/hi_tom_barcodes.csv`](barcodes/hi_tom_barcodes.csv). To use it:
1. Copy the CSV
2. Replace the `Sample` column with your actual sample names
3. Remove rows for unused wells
4. Use this file as input for the demultiplexing step in the notebook

## Methods & Tools Used

- **[fastp](https://github.com/opengene/fastp)** for paired-end read merging with quality-aware overlap detection
- **[CRISPResso2](https://github.com/pinellolab/CRISPResso2)** for CRISPR editing quantification
- **Python** (pandas, NumPy, BioPython) for data processing
- **R / ggplot2** integration for publication-quality visualization
- **Illumina MiSeq 2x250bp** amplicon sequencing data


## How to Use

### Option 1: Google Colab (Recommended)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Innocent250/ngs-analysis-pipeline/blob/main/NGS_Analysis_Pipeline.ipynb)

### Option 2: Run Locally
```bash
git clone https://github.com/Innocent250/ngs-analysis-pipeline.git
cd ngs-analysis-pipeline
pip install -r requirements.txt
jupyter notebook NGS_Analysis_Pipeline.ipynb
```

### Input Requirements

1. **Paired-end FASTQ files** (`.fastq` or `.fastq.gz`) from Illumina MiSeq
2. **Barcode CSV file** with columns: `Sample`, `Barcode_L`, `Barcode_R`
3. **Amplicon reference sequence** for your target locus
4. **Guide RNA sequence** (20 nt, without PAM)

## Requirements

See `requirements.txt` for Python dependencies.

## References

- Liu, Q. et al. (2019). Hi-TOM: a platform for high-throughput tracking of mutations induced by CRISPR/Cas systems. *Sci China Life Sci*, 62, 1-7.
- Clement, K. et al. (2019). CRISPResso2 provides accurate and rapid genome editing sequence analysis. *Nature Biotechnology*, 37, 224-226.

## Author

**Innocent Byiringiro**
- PhD Candidate, Plant Sciences, University of Maryland
- Email: ibyiring@umd.edu
- Lab: [Qi Lab - Plant Genome Engineering](https://www.yipingqi.com/)
