# singlecellRNAseq

This repository provides a step-by-step R pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data using Seurat and related packages. The workflow covers loading 10x Genomics matrices, quality control, clustering, sub-clustering, cell type annotation, and pathway analysis.

## Repository structure

- `codes/`: Ordered analysis scripts (step1, step2, ...) that implement the scRNA-seq workflow. Each script builds on outputs from earlier steps.
  - `step1-ReadDatamatrix.R`: Loads 10x Genomics matrices, initializes Seurat objects, and sources custom plotting and utility functions.
  - `step2-CheckCellMarkers.R`: Performs marker-based quality checks and sanity validations.
  - `step3-Assign_cell_type.R`: Assigns major cell types based on marker genes and clustering.
  - `step4-SubclusterAnalysis.R`: Runs sub-clustering within major cell types to refine populations.
  - `step5.x` scripts: Inspect specific lineages (T cells, dendritic cells, monocytes/macrophages, B cells, mast cells) and categorize states such as cycling or naïve/memory phenotypes.
  - `step6_map_subcelltype_back.R`: Maps refined subcell types back to the integrated object for downstream analysis.
  - `step7-subclusterAnalysis_for_modified_main_celltype.R` and `step8*` scripts: Iterate on sub-clustering and cell-type assignment with updated annotations.
  - `step9_map_subcelltype_back.R`: Updates the integrated object with the latest subcell annotations.
  - `step10_GO_and_KEGG_analysis_for_subclusters.R`: Performs pathway enrichment (GO and KEGG) on identified subclusters.
  - `step11.color_graphs.R`: Generates publication-ready color palettes and visualizations for clusters and annotations.
  - `step12.Monocle2_for_Cd4_Tcells.R` / `step12.Monocle2_for_Cd8_Tcells.R`: Runs trajectory analysis with Monocle2 for CD4 or CD8 T cell subsets.
- `functions/`: Reusable R utilities sourced by the pipeline (custom Seurat plotting helpers, gene symbol conversions, KEGG annotations, cell cycle genes, etc.).
- `README.md`: This document.

## Prerequisites

- R with access to common scRNA-seq packages. The pipeline explicitly loads:
  - `Seurat`, `tidyverse`, `ggplot2`, `sctransform`, `glmGamPoi`, `DoubletFinder`, `patchwork`, `clusterProfiler`, `org.Hs.eg.db`, `harmony`, `celldex`, `RColorBrewer`, `future`, `parallel`.
- Custom functions accessible in `functions/` (e.g., `custom_seurat_functions.R`, `PropPlot.R`, `SubClusterPropPlot.R`, `convertHumanGeneList.R`). Update the `source()` paths in the scripts if your directory layout differs.
- Input data organized as 10x Genomics expression matrices under `datamatrix/` (each sample in its own folder). The scripts expect to be run from a project working directory that also stores intermediate `.Rds` outputs.

## Usage

1. **Set project metadata**: In `step1-ReadDatamatrix.R`, fill in the project name, collaborator, and data location. Adjust the working directory (`validate_workdir(...)`) to point to your project folder.
2. **Install dependencies**: Ensure all required R packages and custom functions are available in your environment. The script will stop if any are missing.
3. **Load data and initialize objects**: Run `step1-ReadDatamatrix.R` to read 10x matrices, optionally cache Seurat objects (`sceList.raw.Rds`), and set up utilities for downstream plots.
4. **Perform QC and annotation**: Execute `step2-CheckCellMarkers.R` and `step3-Assign_cell_type.R` to review marker expression and assign major cell types.
5. **Refine clusters**: Use `step4-SubclusterAnalysis.R` and the `step5.x` lineage-focused scripts to explore subpopulations and states. Map refined annotations back to the integrated object with `step6_map_subcelltype_back.R` or `step9_map_subcelltype_back.R` as appropriate.
6. **Iterative improvements**: Apply `step7` and `step8*` scripts to rerun sub-clustering with updated labels or thresholds.
7. **Pathway and trajectory analysis**: Run `step10` for GO/KEGG enrichment on subclusters, `step11` to finalize color schemes/figures, and `step12` scripts for Monocle2 trajectories of CD4/CD8 T cells.

## Tips

- Many scripts assume relative paths from a consistent working directory; verify paths before execution.
- To speed up Seurat steps, you can enable parallelization with the commented `future` plan options inside `step1-ReadDatamatrix.R`.
- Intermediate `.Rds` files (e.g., `sceList.raw.Rds`) allow you to resume the pipeline without re-reading raw matrices.

## Support

If you encounter issues or need to adapt the pipeline to a new dataset, open an issue with details about your data structure, package versions, and the script/step where the problem occurs.
