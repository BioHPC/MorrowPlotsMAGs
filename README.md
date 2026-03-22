# Morrow Plots Metagenomics Pipeline

Scripts accompanying the manuscript:

> **First Metagenome-Assembled Genomes from the Historic Morrow Plots Reveal Management-Associated Dominance of Archaeal Ammonia Oxidizers**
>
> Viet D. Nguyen, Chushu Gao, Cory Gardner, Zhine Wang, Andrew J. Margenot, Laibin Huang\*, and Tae-Hyuk Ahn\*

## Repository Structure

```
├── Trimming_QC_scripts/        # Read quality control and trimming
│   ├── 01_fastqc.sh
│   └── 02_trim.sh
├── CoAssembly_Binning_Scripts/  # Co-assembly, binning, and annotation
│   ├── 01_make_groups.sh
│   ├── 02_coassemble.sh
│   ├── 03_map_and_depth.sh
│   ├── 04_binning.sh
│   ├── 05_binmap_DAS_Tool.sh
│   ├── 06_checkm.sh
│   ├── 07_rrna_trna_check.sh
│   └── 08_kraken2_group.sh
├── Tree_Script/                 # Phylogenetic tree construction
│   ├── 01_gtdbtk.sh
│   ├── 02_build_tree.sh
│   ├── 03_gtdbtk_functional_tree.sh
│   └── 04_iq_tree.sh
└── Figure_Scripts/              # Figure and table generation
    ├── Figure_2_Panel_B_and_C.R
    ├── Supplementary_Figure_S1.R
    ├── Supplementary_Table_S1.sh
    └── Table_1.sh
```

## Prerequisites

The pipeline requires the following tools. Version numbers reflect those used in the manuscript; other compatible versions may also work.

| Tool | Version | Purpose |
|------|---------|---------|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.5 | Read quality assessment |
| [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.36 | Adapter and quality trimming |
| [MultiQC](https://multiqc.info/) | — | QC report aggregation |
| [MEGAHIT](https://github.com/voutcn/megahit) | 1.2.9 | Metagenomic co-assembly |
| [BBMap](https://sourceforge.net/projects/bbmap/) | 39.33 | Read mapping |
| [samtools](http://www.htslib.org/) | — | BAM sorting and indexing |
| [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) | 2.12.1 | Genome binning |
| [MaxBin2](https://sourceforge.net/projects/maxbin2/) | 2.2.7 | Genome binning |
| [CONCOCT](https://github.com/BinPro/CONCOCT) | 1.1.0 | Genome binning |
| [DAS Tool](https://github.com/cmks/DAS_Tool) | 1.1.7 | Bin refinement |
| [CheckM](https://github.com/Ecogenomics/CheckM) | 1.2.4 | Genome quality assessment |
| [Barrnap](https://github.com/tseemann/barrnap) | 0.9 | rRNA gene prediction |
| [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) | 2.0.12 | tRNA gene prediction |
| [Kraken2](https://github.com/DerrickWood/kraken2) | — | Taxonomic classification |
| [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | 2.5.2 | MAG taxonomy (GTDB r226) |
| [FastTree](http://www.microbesonline.org/fasttree/) | — | Phylogenetic tree construction |
| [IQ-TREE](http://www.iqtree.org/) | 3.0.1 | Maximum-likelihood phylogenetics |
| [QUAST](https://github.com/ablab/quast) | — | Assembly quality statistics |
| [seqtk](https://github.com/lh3/seqtk) | — | FASTA/FASTQ manipulation |
| [seqkit](https://bioinf.shenwei.me/seqkit/) | — | Sequence toolkit |
| [GNU Parallel](https://www.gnu.org/software/parallel/) | — | Parallel job execution |
| [R](https://www.r-project.org/) | 4.5.2 | Statistical analysis and visualization |

R packages: `data.table`, `dplyr`, `stringr`, `ggplot2`, `patchwork`, `scales`, `tibble`

## Pipeline Overview

Scripts within each folder are numbered and should be executed sequentially. Edit the `BASE`, `RAW`, and `WORK` path variables at the top of each script to match your local directory layout before running.

---

### 1. Trimming and Quality Control (`Trimming_QC_scripts/`)

**`01_fastqc.sh`** — Performs initial quality assessment of raw paired-end metagenomic reads (Illumina NovaSeq) using FastQC v0.11.5 to evaluate per-base quality scores, sequence length distributions, GC content, and potential adapter contamination. FastQC analyses are performed on all forward (R1) and reverse (R2) reads for each sample.

**`02_trim.sh`** — Quality filtering and adapter trimming using Trimmomatic v0.36 in paired-end mode. Illumina adapter sequences are removed using the TruSeq3-PE-2 adapter reference file. A sliding window approach trims low-quality regions with a window size of 4 bases and a minimum average Phred score threshold of 30. Reads shorter than 70 bp after trimming are discarded. Only high-quality paired reads are retained for downstream analyses. Post-trimming quality is verified with FastQC, and results are aggregated with MultiQC.

> **Note:** The adapter file path is auto-detected from your Trimmomatic installation. If detection fails, set the `TRIMMOMATIC_ADAPT` environment variable to the full path of `TruSeq3-PE-2.fa`.

---

### 2. Co-assembly and Genome Binning (`CoAssembly_Binning_Scripts/`)

Co-assembly and genome reconstruction are performed using a stepwise workflow in which scripts are executed sequentially.

**`01_make_groups.sh`** — Groups trimmed reads by treatment, generating per-group read lists and a master table linking samples to co-assembly groups.

**`02_coassemble.sh`** — Co-assembles reads within each group using MEGAHIT v1.2.9 with a minimum contig length of 2,000 bp. Assemblies are generated independently for each treatment group using the `meta-large` preset optimized for complex metagenomic data.

**`03_map_and_depth.sh`** — Maps reads back to assembled contigs using BBMap to generate coverage profiles. Resulting BAM files are sorted and indexed using samtools, and contig-level depth information is computed using `jgi_summarize_bam_contig_depths` (from MetaBAT2) for downstream binning.

**`04_binning.sh`** — Performs genome binning using three complementary algorithms:
- **MetaBAT2** — run with a minimum contig length of 2,000 bp
- **MaxBin2** — run with a probability threshold of 0.8 and minimum contig length of 2,000 bp
- **CONCOCT** — applied after fragmenting contigs into 10 kbp segments and computing coverage tables

**`05_binmap_DAS_Tool.sh`** — Integrates and refines bin sets from the three binning methods using DAS Tool with the `blastp` search engine and a score threshold of 0.6, producing a non-redundant set of metagenome-assembled genomes (MAGs).

**`06_checkm.sh`** — Assesses genome quality using CheckM v1.2.4, estimating completeness and contamination based on lineage-specific marker genes and generating summary statistics for all MAGs.

**`07_rrna_trna_check.sh`** — Identifies ribosomal RNA genes using Barrnap v0.9 and transfer RNA genes using tRNAscan-SE v2.0.12, to evaluate gene content completeness and support genome quality assessment following MIMAG standards. 

**`08_kraken2_group.sh`** — Generates treatment-level taxonomic profiles by classifying pooled reads using Kraken2 against the `k2_standard_20230605` database (confidence threshold: 0.05; minimum hit groups: 2). Reads from all samples within each group are analyzed together to obtain group-level taxonomic composition for comparison with MAG-derived profiles.

---

### 3. Phylogenetic Tree Construction (`Tree_Script/`)

Phylogenetic analyses are performed using a stepwise workflow, with scripts executed sequentially.

**`01_gtdbtk.sh`** — Performs taxonomic classification and marker gene identification for all MAGs using GTDB-Tk v2.5.2 with the GTDB release 226 reference database. This step generates standardized taxonomic assignments and multiple sequence alignments (MSAs) for both bacterial (bac120) and archaeal (ar53) marker gene sets.

**`02_build_tree.sh`** — Collects and combines group-level MSAs across all co-assembly groups. For each domain, individual MSAs produced by GTDB-Tk are concatenated into a unified alignment, and phylogenetic trees are constructed using FastTree with the LG+Gamma model.

**`03_gtdbtk_functional_tree.sh`** — For targeted phylogenetic analyses (functional tree), a subset of genomes is selected and processed using GTDB-Tk de novo workflows. Resulting trees are cleaned and pruned to retain selected taxa for downstream visualization.

**`04_iq_tree.sh`** — Generates refined functional trees using IQ-TREE based on subset MSAs, with ModelFinder Plus (MFP) for automatic model selection and 1,000 ultrafast bootstrap + 1,000 SH-aLRT replicates for branch support.

---

### 4. Figure and Table Generation (`Figure_Scripts/`)

**`Table_1.sh`** — Compiles read filtering, co-assembly, and MAG recovery statistics (Table 1) from MultiQC, QUAST, and DAS Tool outputs.

**`Supplementary_Table_S1.sh`** — Generates assembly statistics, read mapping summaries, and bin counts for Supplementary Table S1.

**`Figure_2_Panel_B_and_C.R`** — Produces Figure 2 panels B (stacked bar plot of MAG phylum distribution across management groups) and C (boxplots of GC content and genome length across management groups) using CheckM and GTDB-Tk outputs.

**`Supplementary_Figure_S1.R`** — Generates Supplementary Figure S1, comparing phylum-level relative abundance between Kraken2 read-based classification and MAG-derived taxonomic profiles across treatment groups.

## Citation

If you use these scripts or data, please cite:

> Nguyen, V. D., Gao, C., Gardner, C., Wang, Z., Margenot, A. J., Huang, L.\*, & Ahn, T.-H.\* First Metagenome-Assembled Genomes from the Historic Morrow Plots Reveal Management-Associated Dominance of Archaeal Ammonia Oxidizers.
