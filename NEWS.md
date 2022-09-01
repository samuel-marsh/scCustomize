# scCustomize 0.8.0 (2022-06-XX)  
## Added
- Added `mito_name` parameter to `QC_Plots_Mito` to allow for custom specification of meta data column name that contains mitochondrial information.
- Added `QC_Plots_Combined_Vln()` function to return patchwork layout of 3 QC plots.
- Added Rhesus Macaque (macaca mulatta) to the accepted species list for `Add_Mito_Ribo_Seurat()` and `Add_Mito_Ribo_LIGER()` ([#28](https://github.com/samuel-marsh/scCustomize/issues/28)).
- Added `alpha_exp` and `alpha_na_exp` parameters to `FeaturePlot_scCustom` to allow for control of color scale transparency ([#21](https://github.com/samuel-marsh/scCustomize/issues/21)).
- `*_Highlight_Plot` functions can now plot multiple variables simultaneously using either one color for all variables or one color per variable ([#34](https://github.com/samuel-marsh/scCustomize/issues/34)).
- Added parameter `figure_plot` to `DimPlot_scCustom()`.  This removes axes and axes labels and adds axis legend on left bottom corner of plot ([#40](https://github.com/samuel-marsh/scCustomize/issues/40)).
- Added parameter `plot_legend` to `Stacked_VlnPlot`.  This solves issue with returning only one shared legend across all features being plotted ([#48](https://github.com/samuel-marsh/scCustomize/issues/48)).
- Added `Add_Cell_Complexity_Seurat` and `Add_Cell_Complexity_LIGER` functions to add cell QC complexity/novelty metric (log10(Genes) / log10(UMIs)).  
- Added 3 new CellBender functions `Read_CellBender_h5_Mat`, `Read_CellBender_h5_Multi_Directory`, `Read_CellBender_h5_Multi_File` to enable easy reading of new CellBender output files.
- Vignettes/Website updated with new function examples.  

## Changed
- **BREAKING CHANGE** Function name for iterative `VlnPlot` has been changed to `Iterate_VlnPlot_scCustom` to reflect that it now uses `VlnPlot_scCustom` to generate plots. 
- `QC_Plot_*` functions now use `VlnPlot_scCustom` internally to unify color scheme and rasterization parameters.
- `*_Highlight_Plot` functions no longer display "Unselected" in plot legend and use `DimPlot_scCustom` to generate plots ([#34](https://github.com/samuel-marsh/scCustomize/issues/34)).
- Updated Marsh et al., 2022 citation in vignettes.
- Have begun to move information, warning, and error messages to rlang/cli framework for clarity and style.

## Fixes
- Fixed DESCRIPTION file to specify colorway version upon installation ([#25](https://github.com/samuel-marsh/scCustomize/pull/25)).
- Fixed bug preventing `low_cutoff` from plotting via `QC_Plots_Mito`.
- Fixed bug in `Clustered_DotPlot` that prevented setting identity colors ([#29](https://github.com/samuel-marsh/scCustomize/issues/29)).
- Fixed bug in `FeaturePlot_scCustom` that returned NULL when setting `combine = FALSE` ([#31](https://github.com/samuel-marsh/scCustomize/issues/31)).
- Fixed bug in `Seq_QC_Plot_*` functions which resulted in groups being plotted out of order when specifying `plot_by` parameter.
- Fixed bug in `Seq_QC_Plot_*` functions that created color palette error when color palettes were not being used.
- Fixed bug in `DimPlot_scCustom` that caused mismatch of colors between plots when using `split.by` if one of the plots was missing 1 or more of the `group.by` levels ([#37](https://github.com/samuel-marsh/scCustomize/issues/37)).
- Fixed bug in `VlnPlot_scCustom` that caused raster warning messages to be displayed twice ([#42](https://github.com/samuel-marsh/scCustomize/issues/42)).
- Fixed bug in `Iterate_PC_Loading_Plots` that caused error when specifying current directory with `file_path = NULL` or `file_path = ""`
- Fixed bug in `DotPlot_scCustom` that prevented plotting of features in meta.data slot ([#44](https://github.com/samuel-marsh/scCustomize/issues/44)).
- Fixed error messaging/reporting in `Stacked_VlnPlot` when no supplied features were present.
- Fixed bug in `plotFactors_scCustom` that was ignoring provided file name.
- Fixed bug in `plotFactors_scCustom` that caused progress to only display progress up to 50% even when it was fully complete.
- Fixed bug in `Clustered_DotPlot` that resulted in error related to color palettes if number of clusters was greater than 36 ([#49](https://github.com/samuel-marsh/scCustomize/issues/49)).
- Fixed bug in `Add_Mito_Ribo_LIGER` that resulted custom column names (e.g. `mito_name = "pct.mt"`) being disregarded and also therefore issue with `overwrite` parameter. ([#51](https://github.com/samuel-marsh/scCustomize/issues/51)).
- Fixed bug in `Store_Misc_Info_Seurat` that prevented function from working.
- Typo/styling fixes.
 

# scCustomize 0.7.0 (2022-01-10)  
## Added
- Added `VlnPlot_scCustom` function.
- Added raster support to `Stacked_VlnPlot`
- Added `make_unique` parameter to `Extract_Top_Markers` function.
- Added `Clustered_DotPlot` function.
- Added Drosophila Melanogaster as default species option in `Add_Mito_Ribo_Seurat` and `Add_Mito_Ribo_LIGER`.

## Changed
- Now requires Seurat v4.0.6 (instead of v4.0.5) to support ability to rasterize points in `VlnPlot`.
- viridis color palette shortcuts now contain palettes with 30 colors (increased from 10).

## Fixes
- Fixed `Read_Metrics_10X` errors that occurred due to differing outputs depending on Cell Ranger version or type of assay.
- Added direct `importFrom` for `DefaultDimReduc` from SeuratObject package to avoid potential errors.
- Fixed typos/styling in function documentation.
 

# scCustomize 0.6.3 (2021-12-16)  
## Fixes
- Fixed `Read_Metrics_10X` errors that occured due to differing outputs depending on Cell Ranger version or type of assay.
- Added direct `importFrom` for `DefaultDimReduc` from SeuratObject to avoid potential errors.
 

# scCustomize 0.6.2 (2021-12-01)  
## Fixes
- Fixed barcode name duplication checks in `Merge_Sparse_Data_All`. ([#8](https://github.com/samuel-marsh/scCustomize/issues/8))
- Fixed package imports in DESCRIPTION to avoid installation errors.
- Fixed NULL check in `Read_Metrics_10X`, `Read10X_Multi_Directory`, and `Read10X_h5_Multi_Directory`.
 

# scCustomize 0.6.1 (2021-11-19)
## Added
- Added plot spacing control to `StackedVlnPlot` with parameters `plot_spacing` and `spacing_unit`. ([#6](https://github.com/samuel-marsh/scCustomize/issues/6))
- Added `scCustomize_Palette` function select palette to use (simplify internal code).

## Changes
- Changed citation info to reflect global DOI and not version DOI.

## Fixes
- Restore package color palette defaults to `Iterate_VlnPlot`.  
- Fix `Iterate_...` function checks for file path parameter if `file_path = NULL`.
  
# scCustomize 0.6.0 (2021-11-16)
## Added
- scCustomize is public!!  Version 0.6.0 is released!

## Changes
- Many function names have changed since private release see reference page/manual for updated function names.
