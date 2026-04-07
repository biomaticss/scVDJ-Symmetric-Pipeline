# scVDJ Symmetric Pipeline: Single-Cell Antibody Repertoire Profiling

##  Background & Significance
High-throughput single-cell V(D)J sequencing is paramount for understanding adaptive immune responses and discovering therapeutic monoclonal antibodies (mAbs). Traditional bulk sequencing loses the critical native pairing of Heavy (VH) and Light (VK/VL) chains. 

This repository provides a custom, highly optimized bioinformatics pipeline for processing structurally barcoded single-cell antibody libraries. It is specifically designed to evaluate immune repertoires post-immunization, accurately reconstruct native VH-VK pairs, profile Somatic Hypermutations (SHMs), and isolate dominant, functional clonotypes for downstream affinity validation and recombinant expression.

##  Library Architecture & Assay Design
Our custom single-cell droplet microfluidic assay utilizes a symmetric physical linkage strategy to maintain VH-VK pairing. The amplicon architecture is strictly defined as:
`[BC1] - [VH] - [BC3] - [VK] - [BC2]`

* **`[BC3]` (Central Hub):** The shared barcode connecting heavy and light chains within the same droplet.
* **`[BC1]` & `[BC2]` (Flanking Identifiers):** Used for molecular indexing, over-amplification screening, and serving as physical validation tags for downstream Sanger sequencing.

To accurately resolve this long construct, the sequencing is split into two complementary groups:
1. **Topology Group (PE150):** Reads the linkages `[BC1]-[BC3]` and `[BC3]-[BC2]` to establish a bipartite molecular index.
2. **Repertoire Group (PE300):** Reads the sequence payloads `[VH]-[BC3]` and `[VK]-[BC3]` to capture the full V(D)J recombinations.

##  Pipeline Algorithm & Workflow

### Stage 1: Topological Assembly (Dual-Track Indexing)
* Extracts multi-segment barcodes using ultra-fast Regex pattern matching.
* Dynamically resolves PCR amplification biases and droplet-overload artifacts by evaluating UMI collision frequencies. 
* Outputs a clean, bi-directional `BC1-BC3-BC2` relational database.

### Stage 2: Physical Assembly
* Scans PE300 reads for valid `BC3` hubs.
* Symmetrically pairs VH and VK full-length sequences based on shared droplet hubs.

### Stage 3: High-Resolution Annotation & Pathology Diagnostics
* **MiXCR Integration:** Performs robust alignment against the IMGT reference database to identify V, D, and J gene segments.
* **Region Extraction:** Precisely delineates Framework Regions (FR1, FR2, FR3, FR4) and Complementarity-determining Regions (CDR1, CDR2, CDR3) at the amino acid level.
* **Defect Profiling:** Implements an advanced sliding-window diagnostic to detect internal frameshifts (`_`) and premature stop codons (`*`), strictly distinguishing genuine polymerase slippage from sequencing-edge truncation artifacts.

### Stage 4: Generalized L4 Clonal Clustering (DP Pruning)
To account for PCR point mutations and natural Somatic Hypermutation (SHM) during affinity maturation, the pipeline performs nested clustering:
* **L1-L3 Clustering:** Groups clones by strict identical nucleotide/amino acid sequences.
* **L4 Generalized Clustering:** Utilizes a **Levenshtein Distance Matrix** (Dynamic Programming) to collapse clonotypes sharing identical V/J usage and highly homologous CDR3 loops (default <= 10% edit distance tolerance). 

## Installation & Dependencies

The pipeline is written in Python 3 and heavily relies on multiprocessing for scalability.

### Prerequisites
1. **Python 3.8+**
2. **MiXCR** (Must be installed and accessible in the system `$PATH`) 
   * *Refer to [MiXCR official documentation](https://mixcr.com/) for licensing and installation.*

### Python Packages
```bash
pip install pandas numpy matplotlib seaborn regex
```

##  Usage

Execute the main pipeline using the standard command-line interface. Make sure all raw `.fastq` or `.fastq.gz` files are placed in a designated data directory.

```bash
python src/run.py \
    --sample "Project_mRNA_Vax_01" \
    --species "hsa" \
    --bc1_r1 data/raw/bc1_R1.fastq.gz \
    --bc1_r2 data/raw/bc1_R2.fastq.gz \
    --bc2_r1 data/raw/bc2_R1.fastq.gz \
    --bc2_r2 data/raw/bc2_R2.fastq.gz \
    --vh_fq data/raw/vh_reads.fastq.gz \
    --vk_fq data/raw/vk_reads.fastq.gz \
    --out_dir results/
```

### Key Parameters:
* `-s, --sample`: Base name for all output files.
* `-sp, --species`: Target species for VDJ alignment (`hsa` for human, `mmu` for mouse).
* `-o, --out_dir`: Destination directory for reports and metrics.

##  Outputs & Deliverables

Upon successful execution, the `results/` directory will contain:

1. **`*Client_Master_L4_Clonotypes.csv`**: The definitive output table. Contains fully paired, functionally validated, and L4-clustered clonotypes ranked by UMI frequency. Includes complete CDR/FR sequence breakdowns and structural annotations.
2. **`*Pipeline_Report.md`**: A comprehensive Markdown report detailing data retention rates, barcode collision statistics, and pathology hotspot mapping.
3. **`*Unfiltered_Library_Audit.csv`**: The raw, unpruned database containing all matched reads, including those flagged with stop codons or frameshifts, strictly for library quality control.
4. **`Visualizations/` (Directory)**:
   * `01_Top20_L4_Clones.png`: Barplot of the dominant clones.
   * `04_Accumulation_Curve.png`: Log-scaled clonal expansion trajectory.
   * `05_Mutation_Hotspots_Profiler.png`: A comprehensive diagnostic heatmap tracking truncation, frameshifts, and nonsense mutations across the 5'->3' antibody topology (FR1 to FR4).

##  License & Citation
*( MIT)*

If you utilize this pipeline in your research for antibody discovery or epitope mapping, please consider citing this repository.