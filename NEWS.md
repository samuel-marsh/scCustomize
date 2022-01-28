# scCustomize 0.7.X (2022-XX-XX)  
## Added
- Added `mito_name` parameter to `QC_Plots_Mito` to allow for custom specification of meta data column with mitochondrial information.

## Changed
- `QC_Plot_*` functions now use `VlnPlot_scCustom` internally to unify color scheme and rasterization parameters.

## Fixes
- Fixed DESCRIPTION file to specify colorway version upon installation (#25).
- Fixed bug preventing `low_cutoff` from plotting via `QC_Plots_Mito`.
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
