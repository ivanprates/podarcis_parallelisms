# Podarcis parallelisms: Analysis scripts

Scripts associated with the manuscript:  
**"Repeated phenotypic evolution and its genomic correlates in Mediterranean wall lizards"**  
Ivan Prates et al.

## Directory Structure

### `scripts_R/`
R and shell scripts for phenotypic, colorimetric, and genomic analyses:
- `01a_phenotypic_analyses_2026-08-11.R`: Morphological and colorimetric trait variation, integration, and multivariate trajectories.
- `01b_one_brown_many_greens.R`: Comparative dorsal coloration analyses across seven *Podarcis* species (Figure S2).
- `02_phenotypic_parallelism_Ttests_2026-08.R`: Phenotypic trajectory parallelism tests and pairwise comparisons (Table S2).
- `03_TajimasD.R`: Tajima's D calculations across sliding windows.
- `04_intersect_baypass_and_TajimasD_outliers_03-26.sh`: Intersecting BayPass $XtX$ and Tajima's D outlier regions.
- `05_candidate_gene_sharing.R`: Identifying candidate loci and cross-species overlap.
- `06_candidate_gene_sharing_simulations.R`: Gene-length-weighted null simulations testing candidate gene sharing across species.
- `07_Fig4_2026-08-12.R`: Plotting candidate gene sharing patterns (Figure 4).
- `08_pigmentation_genes.R`: Curation and filtering of pigmentation- and chromatophore-related genes.
- `09_pigmentation_gene_overrepresentation.R`: Fisher's exact tests for pigmentation gene enrichment (Table S5).
- `10_chrom_12_candidate_region.R`: Differentiation and coverage analysis across the chromosome 12 candidate region (Figure S3).
- `11_candidate_gene_dissimilarity_vs_phenotypic_trajectories.R`: Generalized Dissimilarity Modeling (GDM) testing candidate gene vs. phenotypic dissimilarity (Figure S4).
- `12_Fig1.R`: Figure 1 map and layout.
- `13_FigS1_2026-08.R`: Genome-wide differentiation scan plots (Figure S1).
- `14_combine_sample_information_for_Table_S1.R`: Compiling sample and locality metadata (Table S1).

### `scripts_BayPass/`
SLURM job scripts for BayPass population genomic scans and $XtX$ outlier detection.

### `scripts_bioinformatics/`
SLURM job scripts and batch-submission wrappers for raw read quality control, trimming, mapping, variant calling, and VCF filtering.
