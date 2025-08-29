# scCustomize 3.1.4 (2025-08-XX)  
## Added  
- Added cNMF vignette.  
- Added `reduction_name` and `reduction_key` parameters to `Read_Add_cNMF` to allow for control over reduction name (and ability to add results from more than one run of cNMF to the same Seurat object).  

## Changed  
- Internal function `yesno` will now provide Yes/No options in same order everytime to avoid mistaken input.  



## Fixes  
- Fixed errors in plotting with `Factor_Cor_Plot` that occured due to partial matrix plotting and row clustering.  
- Fixed column naming in `Top_Genes_Factor` for Seurat objects to align with format returned for LIGER objects.  



# scCustomize 3.1.3 (2025-08-25)  
## Added  
- Add support for LIGER objects using `Extract_Sample_Meta`.  
- Added new function `Dataset_Size_LIGER` to return data.frame containing cells per dataset in liger object in addition to accompanying sample meta data if desired.  
- Added new function `Get_Reference_LIGER` to return name of dataset matching selection criteria to use as reference for `rliger::quantileNorm()`.  
- Added new function `Add_MALAT1_Threshold` which implements QC procedure from Clarke & Bader (2024). bioRxiv \doi{doi.org/10.1101/2024.07.14.603469}.   
- Added new internal function `Check_Normalized` to confirm data within "data" layer of Seurat object is actually normalized (relevant only for V3/4 style objects/assays).  
- Added new parameters to `Add_Cell_QC_Metrics` (Seurat ONLY) to add module score of IEGs in addition to percent expressing.  
- Added new function `exAM_Scoring` to add module scores for exAM gene sets from Marsh et al., 2022 (\doi{10.1038/s41593-022-01022-8}).  
- Added new function `Proportion_Plot_per_Sample` to plot proportion of cells per sample across a specified condition.  
- Added new parameter `order_by_freq` to `Cluster_Stats_All_Samples`.  It is `TRUE` by default and returned data.frame is ordered by cluster frequency, setting FALSE orders data.frame by cluster order.  
- Updated `Convert_Assay` to allow V5 > V3 conversions even when normalized and scale data is absent ([#236](https://github.com/samuel-marsh/scCustomize/issues/236)).  
- Added `Read_Add_cNMF` to read and add results from cNMF as custom dimensionality reduction.  
- `Top_Genes_Factor` is now S3 generic and works with either Seurat or LIGER objects.  
- `Top_Genes_Factor` now supports `factor = "all"` which will return data.frame containing top X genes for all factors, 1 column per factor.  
- Added parameter `label_selected_features` to `Clustered_DotPlot` to allow for labeling only subset of plotted features.  
- Added new parameter to `Add_Cell_QC_Metrics` to add percentage of lncRNA counts per cell (see `add_lncRNA` parameter).  

  
  
## Changed  
**This release contains a number of BREAKING changes to parameter names:**  
  
- **BREAKING CHANGE** The parameter `num_genes` has been soft-deprecated in `Extract_Top_Markers`.  Please use `num_features` instead.  Using `num_genes` will warn user but still work until scCustomize v3.3.0.  
- **BREAKING CHANGE** The parameter `min_cells` and `min_features` have been soft-deprecated in `Create_CellBender_Merged_Seurat`.  Please use `min.cells` and `min.features` instead.  Using `min_cells` and `min_features` will warn user but still work until scCustomize v3.3.0.  
  
**This release contains a number of BREAKING changes to parameter names to harmonize across scCustomize and Seurat:**  
*Due to large number of functions affected the timeline for full deprecation of these parameters has been extended.  Old parameter names will issue warning but continue to work until v3.3.0.*  
  
- **BREAKING CHANGE** The `group_by` parameter has been soft-deprecated in `Plot_Median_Genes`, `Plot_Median_UMIs`, `Plot_Median_Mito`, `Plot_Median_Other`, `Plot_Cells_per_Sample`, `Percent_Expressing`, `DimPlot_LIGER`, and `Extract_Top_Markers`.  Please use `group.by` instead.  Using `group_by` will warn user but still work until scCustomize v3.3.0.  
- **BREAKING CHANGE** The `group_by_var` parameter has been soft-deprecated in `Proportion_Plot`, `Cluster_Stats_All_Samples`, `Median_Stats`, and `MAD_Stats`.  Please use `group.by` instead.  Using `group_by_var` will warn user but still work until scCustomize v3.3.0. 
- **BREAKING CHANGE** The `split_by` parameter has been soft-deprecated in `Percent_Expressing`, `DimPlot_LIGER`, and internal functions.  Please use `split.by` instead.  Using `split_by` will warn user but still work until scCustomize v3.3.0. 
  
**Non-breaking changes in this release:**  
  
- The following parameters in `plotFactors_scCustom` have been fully deprecated for LIGER objects >= V2: `reorder_datasets` and `reduction_label`.  
- Following prior deprecation warnings the following functions are now fully deprecated and replaced with updated functions: `Add_Cell_Complexity_LIGER`, `Add_Cell_Complexity_Seurat`, `Add_Cell_Complexity_Seurat`, `Add_Mito_Ribo_LIGER`, `Add_Mito_Ribo_Seurat`, `Gene_Present`, `Meta_Present_LIGER`, and `Split_FeatureScatter`.  
- Changed internal function `PercentAbove_Seurat` to match updates to Seurat to appropriately deal with NA values.  
- Changed functionality of several `QC_Plot*` functions to dynamically set nFeature or nCount variable name based on assay specified.  
- Changed default parameter value for `x_lab_rotate` in `Proportion_Plot` from FALSE to TRUE.  
- Changed param `selection.method` to `method` in `VariableFeaturePlot_scCustom` to account for deprecation in Seurat/SeuratObject.  
    

## Fixes  
- Fixed use of chicken as default species in some QC functions.  
- Fixed bug in `Read_Metrics_10X` that caused function failure.  
- Fixed bug in `Subset_LIGER` to ignore cluster column when subsetting based on other meta data variable.  
- Fixed bug in `QC_Histogram` that didn't allow it to properly plot feature data.  
- Fixed bug in `FeaturePlot_scCustom` that prevented `max.cutoff`/`min.cutoff` from being correctly passed when splitting plots ([#228](https://github.com/samuel-marsh/scCustomize/issues/228)).  
- Fixed check for file extension in `Iterate_PC_Loading_Plots`.  New internal function `check_extension` to ease these checks package-wide.  
- Fixed bug in behavior of `Extract_Top_Markers` when sorting the markers by "p_val_adj" that was selecting genes with highest p values instaed of lowest ([#229](https://github.com/samuel-marsh/scCustomize/issues/229)).  
- Fixed rotation of x-axis text in `Proportion_Plot`.  
- Added check for correct input format in `Extract_Top_Markers`.  
- Added check to `Plot_Median_Genes`, `Plot_Median_UMIs`, `Plot_Median_Mito`, and `Plot_Median_Other` to ensure that `group.by` and `sample_col` are different and provide informative error message if they are the same ([#233](https://github.com/samuel-marsh/scCustomize/issues/233)).  
- Fixed bug in `VariableFeaturePlot_scCustom` that prevented function from running.  
- Fixed `Factor_Cor_Plot` failure when using Seurat object due to lack of `reduction` parameter.  
- Fixed issue in `Clustered_DotPlot` when feature has zero expression in any cells with new `nan_error` parameter ([#178](https://github.com/samuel-marsh/scCustomize/issues/178)).  
- Fixed error in `Barcode_Plot` that prevented plotting when using newer versions of DropletUtils.  
- Code styling and typo fixes.  





# scCustomize 3.0.1 (2024-12-18)  
## Added  
- Added new parameters `output_width` and `output_height` to the `Iterate_*` family of plotting functions ([#217](https://github.com/samuel-marsh/scCustomize/issues/217)).  
  
  
## Changed  
    

## Fixes  
- Fixed bug in `Random_Cells_Downsample` that prevented setting identity using the `group.by` parameter.  
- Fixed bug in `Cell_Highlight_Plot` that didn't pass the reduction parameter properly ([#216](https://github.com/samuel-marsh/scCustomize/issues/216)).  
- Fixed bug when retrieving ensembl IDs for IEGs.  
- Fixed bug that prevented using `return_plots` in iterative plotting functions([#217](https://github.com/samuel-marsh/scCustomize/issues/217)).  




# scCustomize 3.0.0 (2024-12-05)  
## Added  
**Major Updates to Functionality with rliger Package:**  
*Added new utility functions to interact with liger v2.0.0+ object format change:*  
    - `Subset_LIGER` to quickly subset by cluster or other meta data variable.  
    - `Cells_by_Identities_LIGER` to extract list of barcodes sorted by values within given meta data column.  
    
*Extended the following Seurat/SeuratObject generic functions to work seamlessly with liger objects:*  
    - `Cells` to extract vector of all cells or list vectors of cells by dataset.  
    - `Features` to extract vector of all features or list vectors of features by dataset.  
    - `WhichCells` to extract vector or list of cells matching identity criteria.  
    - `Embeddings` to extract matrix containing dimensionality reduction embeddings or iNMF h.norm matrix.  
    - `Idents` and `Idents<-` to extract and set default identities/clusters.  
    
*Updated functions to interact with both old and new style liger objects:*  
    - `plotFactors_scCustom()`  
    - `Fetch_Meta`  
    - `Top_Genes_Factor`  
    - `Add_Mito_Ribo`  
    - `Add_Cell_Complexity`  
    - `DimPlot_LIGER`  
    - `Variable_Features_ALL_LIGER`  
    - `Feature_Present`  
    
*New functions compatible with old and new style liger objects:*  
    - Added new function `Add_Hemo` to add hemoglobin gene percentage for QC.  Also added as parameter to `Add_Cell_QC_Metrics`.  `Add_Hemo` supports all default species: (human, mouse, marmoset, zebrafish, rat, drosophila, rhesus macaque, and chicken) and works with both Seurat and liger objects.  
    
*New scCustomize generics to function across both Seurat and Liger objects:*  
    - `Add_Hemo` (see above).  
    - `Rename_Clusters` now S3 generic for setting new active.ident (Seurat) or defaultCluster (Liger).  
    
*New functions for Seurat and rliger v2.0.0+ only:*  
    - Added new function `Find_Factor_Cor` to return correlation matrix between factor gene loadings from liger or Seurat object.  
    - Added new function `Factor_Cor_Plot` to plot positive correlations from liger or Seurat object.  
    
*Updated functions to recommend new rliger equivalents for users with rliger v2.0.0+:*  
    - `as.LIGER`  
    - `as.Seurat`  

  
**General scCustomize Updates:**  
*New functions:*  
- Added new function `Add_Hemo` to add hemoglobin gene percentage for QC.  Also added as parameter to `Add_Cell_QC_Metrics`.  `Add_Hemo` supports all default species: (human, mouse, marmoset, zebrafish, rat, drosophila, and rhesus macaque) and works with both Seurat and liger objects.  
- Added new function `seq_zeros()` to create sequences with preceding zeros.  
- Added new function `Read_Metrics_CellBender` to read in the summary metrics csv file produced by CellBender.  Can either read all metrics files from parent directory of output folders or a single metrics file.  
- Added `Updated_MGI_Symbols` to check for update gene names/symbols in mouse data ([#202](https://github.com/samuel-marsh/scCustomize/issues/202)).  
- Added plotting function `Proportion_Plot` to plot pie chart or bar chart of proportion (or total counts) of cells in each identity class.  
- Added new function `Random_Cells_Downsample` to return either a vector or list with randomly downsampled cells for each identity class.  
- Added new function `Cells_per_Sample` to quickly return data.frame with just number of cells per sample.  
  
*Updated functions:*  
- Added new parameters `data_name` and `overwrite` to `Add_Alt_Feature_ID` to support new storage location.  
- Added `cells` parameter explicitly to `FeatureScatter_scCustom`.  
- Added Chicken (Gallus gallus) to default species for QC functions.  Thanks @dpearton; ([#176](https://github.com/samuel-marsh/scCustomize/issues/176)).  
- Added new plotting function `SpatialDimPlot_scCustom`.  Thanks for encouragement @puapinyoying @nina-hahn ([#160](https://github.com/samuel-marsh/scCustomize/issues/160)).  
- Added ability of `Read_Metrics_10X` to read a single metrics csv file and return data formatted the same way as when reading multiple files.  
- Added parameter `cutoff_line_width` to the `QC_Plot_*` family of plots to control line thickness of cutoff lines.  
- `Cluster_Stats_All_Samples` now returns data.frame with row order reflecting the frequency of cells.  
- `Add_Mito_Ribo` now supports datasets aligned to multi-species reference genomes ([#184](https://github.com/samuel-marsh/scCustomize/issues/184)).  
- Added parameter `add_prop_plot` to `DimPlot_scCustom` to return plot showing number or percent of cells per identity along with the DimPlot.  
- Added optional parameter `colors_use_assay2` to `FeaturePlot_DualAssay` which allows for specification of different palettes for the two plots ([#182](https://github.com/samuel-marsh/scCustomize/issues/182)).  
- Added new folder and scripts (see "data-raw/" on GitHub) detailing the creation of gene lists used in `Add_Cell_QC_Metrics`.  
- Added ensembl ID support for percent hemoglobin, msigdb, and IEG gene sets ([#186](https://github.com/samuel-marsh/scCustomize/issues/186)).  
- Add verbosity parameter to `Store_Misc_Info_Seurat` and `Store_Palette_Seurat`.  
- Explicitly reveal the `reduction` parameter in `Cluster_Highlight_Plot` and `Meta_Highlight_Plot` ([#198](https://github.com/samuel-marsh/scCustomize/issues/198)).  
- Added `show_row_names` `show_column_names`, `column_names_side`, `row_names_side`, `legend_position`, `legend_orientation`, `show_ident_legend`, and `show_ident_colors` parameters to `Clustered_DotPlot`.  Thanks for idea and code @johnminglu ([#199](https://github.com/samuel-marsh/scCustomize/issues/199)).  
- Updated `Split_Vector` to allow user to specify number of chunks or size of chunks for splitting vector.  
- Update `RenameClusters` with additional parameters to enable storage of both old idents and new idents in meta.data within the function.  
- Update `Add_Cell_QC_Metrics.Seurat` to explicitly reveal `list_species_names` parameter.  
- Added new vignette for spatial plotting.  
- Added new and expanded vignette on use of object QC functions for better clarity on these functions and their uses (previously was part of QC Plotting & Helpers/Utilities Vignettes).  Plotting elements of QC Plotting vignette are unchanged.  

  
## Changed  
- **BREAKING CHANGES** `Add_Top_Gene_Pct_Seurat` is now S3 generic that works with both Seurat and liger objects and has been renamed `Add_Top_Gene_Pct`.  
- `Add_Cell_QC_Metrics` is now S3 generic and works with both Seurat and liger objects.  
- Changed storage location for `Add_Alt_Feature_ID` to `@misc` slot of object for safer storage across object filtering.  
- Added error check in `as.anndata` to explicitly check for installation of anndata before starting conversion ([#162](https://github.com/samuel-marsh/scCustomize/issues/162)).  
- Updated `Plot_Median_Genes`, `Plot_Median_UMIs`, `Plot_Median_Mito`, `Plot_Median_Other`, `Plot_Cells_per_Sample` to understand "ident" as grouping variable.  
- Updated `Store_Misc_Info_Seurat` to use Seurat accessor/setter function `Seurat::Misc()`.  
- Updated documentation for `sample_names` in `Read_CellBender_h5_Multi_File` to clarify parameter behavior (related to ([#208](https://github.com/samuel-marsh/scCustomize/issues/208))).  
- Updated `Read_Metrics_10X` to support adjusts to metrics summary format and metric names in output from Cell Ranger v9+.  
- Some reorganization of R/ directory/scripts.  
   

## Fixes  
- Nebulosa plotting functions `Plot_Density_Custom` and `Plot_Density_Joint_Only` have been re-enabled for users with ggplot2 v3.5.0 following Nebulosa v1.12.1 update patch.     
- Fixed bug causing error in `Add_Cell_QC_Metrics` when `overwrite = TRUE` ([#165](https://github.com/samuel-marsh/scCustomize/issues/165)).  
- Fixed wrong description of parameter in manual entry for `DotPlot_scCustom` ([#158](https://github.com/samuel-marsh/scCustomize/issues/158)).  
- Fixed several potential errors in `as.anndata` from Seurat conversion that previously caused failures ([#168](https://github.com/samuel-marsh/scCustomize/issues/168)).  
- Fixed errors in `Create_Cluster_Annotation_File` if for file path and csv name errors.  
- Fixed error when using `plot_median` and more than one feature in `VlnPlot_scCustom` ([#169](https://github.com/samuel-marsh/scCustomize/issues/169)).  
- Fixed bug while collecting legends for `DimPlot_scCustom` due to changes in guides updated with ggplot2 v3.5.0 ([#171](https://github.com/samuel-marsh/scCustomize/issues/171)).  
- Fixed error in `Add_Sample_Meta` that still errored when setting `na_ok = TRUE`.  
- Fixed errors in `Plot_Median_*` family that caused issues when `group_by` parameter was NULL.  
- Fixed errors in `FeaturePlot_scCustom` when setting `combine = FALSE`.  
- Fixed bug in `DimPlot_scCustom` that could cause blank plot when rasterizing points.  
- Fixed bug in `MAD_Stats` that didn't respect `mad_num` parameter ([#183](https://github.com/samuel-marsh/scCustomize/issues/183)).  
- Fixed bugs in `MAD_Stats` that could cause issues if `mad_num` was less than or equal to 0 and returned error if setting `group_by_var` to "ident".  
- Replaced lingering instances of deprecated tidyr code .data[["var"]] with update `all_of`/`any_of` syntax.  
- Fixed issue that could occur with some meta data modifying functions due to column name collisions in internals of function ([#193](https://github.com/samuel-marsh/scCustomize/issues/193)).  
- Fixed issue that caused error when using `Cluster_Highlight_Plot` with `split.by` parameter ([#201](https://github.com/samuel-marsh/scCustomize/issues/201)).  
- Added check and informative error message to `Convert_Assay` ([#205](https://github.com/samuel-marsh/scCustomize/issues/205)).  
- Fixed issue with anndata conversion and Seurat V5 objects ([#195](https://github.com/samuel-marsh/scCustomize/issues/195)).  
- Fixed issue with `Updated_HGNC_Symbols` due to change in URL path for gene names ([#209](https://github.com/samuel-marsh/scCustomize/issues/209)).  
- Fixed bug in `DimPlot_scCustom` when `split.by` and `label.box = TRUE`.  
- Fixed bug in `DiscretePalette_scCustom` that didn't error when supplying invalid palette names.  
- Fixed bug in `DimPlot_LIGER` that provided uniformative error message when changing the default cluster ident.  
- Spelling and style fixes.  Thanks @kew24.  



# scCustomize 2.1.2 (2024-02-27)  
## Added  
- None.  


  
## Changed  
- None.   
   

## Fixes  
- Fixed breaking bug (and CRAN check errors) in `Plot_Density_Custom` and `Plot_Density_Joint_Only` due to error in Nebulosa following ggplot2 v3.5.0 update.  Functionality will be restored when Nebulosa package is updated.   


# scCustomize 2.1.1 (2024-02-23)  
## Added  
- None.  


  
## Changed  
- None.   
   

## Fixes  
- Fixed breaking bug in `as.liger` that prevented function from running properly.  




# scCustomize 2.1.0 (2024-02-21)  
## Added  
- Added `title_prefix` parameter to `Iterate_DimPlot_bySample` to unify with `Meta_Highlight_Plot`.  
- Added function `Split_Vector` to split vector in chunks of predetermined sizes.  
- Added new function `Updated_HGNC_Symbols` to update human gene symbols.  After first use does not require internet connection.  
- Added command logging to QC metric-related commands using `Seurat::LogSeuratCommand()`.  
- Added parameter `plot_legend` to `plotFactors_scCustom` to allow for removal to legend in factor loading plots.  
- Added new functionality to `Iterate_FeaturePlot_scCustom` to allow for plotting multiple plots per page when saving to single PDF document (see new parameters: `features_per_page` and `landscape`.  
- Added `LIGER_Features` utility function for LIGER objects (analogous to `Seurat::Features()`).  
- Added new generic `as.LIGER()` as enhanced method for conversion of Seurat objects or lists of Seurat objects to single LIGER object.  
- Added new generic `as.anndata()` to support conversion of Seurat and LIGER objects to anndata format.  
- Added function `Convert_Assay()` to perform easy conversions of Seurat assays between V3/4 (Assay) and V5 (Assay5) formats.  
- Added parameter `assay_type` to allow manual control of assay type when creating Seurat object from LIGER using `Liger_to_Seurat`.  Now `as.Seurat`.  
- Added param `grid_color` to `Clustered_DotPlot` to control the color of grid lines (default is no grid lines).  
- Added ability to split identities in `Clustered_DotPlot` by additional variable and maintain expression information.  
- Added `Split_Layers()` function for V5 objects.  
- Added `Add_Alt_Feature_ID` to add alternative feature ids to an Assay5 meta.data slot.  


  
## Changed  
- **BREAKING CHANGES** Several methods have been moved to S3 generics to function for both Seurat and LIGER objects using single function name and therefore some function names have changed.  Old functions will give deprecation warning and direct users to new functions.  
    - `Add_Mito_Ribo()` replaces `Add_Mito_Ribo_Seurat` and `Add_Mito_Ribo_LIGER`.  
    - `Add_Cell_Complexity()` replaces `Add_Cell_Complexity_Seurat` and `Add_Cell_Complexity_LIGER`.  
- **BREAKING CHANGES** `Meta_Present_LIGER` has been deprecated and wrapped inside of `Meta_Present`.  
- **SOFT-DEPRECATION** The function `Liger_to_Seurat()` has been soft-deprecated.  It is replaced by new extension of Seurat generic `as.Seurat` with added support for Liger objects, using all the same parameters as `Liger_to_Seurat`.  Full deprecation will occur in v2.2.0.  
- **SOFT-DEPRECATION** The function `Gene_Present` has been soft-deprecated.  It is replaced by `Feature_Present` which functions identically but better reflects that features present may also be proteins.  Full deprecation will occur in v2.2.0.  
- Parameter `legend` in `Iterate_DimPlot_bySample` has been inverted to `no_legend` to match `Meta_Highlight_Plot` parameters.  
- Updated `Liger_to_Seurat()` for compatibility with Seurat V5 structure ([#153](https://github.com/samuel-marsh/scCustomize/issues/153)).  Now part of `as.Seurat`.
- Default color palette change from v2.0.0 when number of groups is between 3-8 has been reverted.  Polychrome palette is default when number of groups is between 3-36.  
- In preparation of upcoming overhaul of rliger package, added package version checks to current rliger functions in order to prevent breaking errors.  Next update v2.2.0 will add cross-functionality between rliger package versions ([#161](https://github.com/samuel-marsh/scCustomize/issues/161)).  
   

## Fixes  
- General typo and style fixes.  
- Fixed point size check in some QC functions to avoid unnecessary error message.  
- Fixed redundant warning messages in `Stacked_VlnPlot` due to rasterization defaults.  
- Fixed issue setting `alpha_na_exp` appropriately in `FeaturePlot_scCustom`.  
- Fixed issue setting `alpha_exp` between Seurat versions 4 and 5 ([#144](https://github.com/samuel-marsh/scCustomize/issues/144)).  
- Fix duplicate legends in `DimPlot_scCustom` when levels are missing from a split plot.  
- Fixed bug in `FeaturePlot_scCustom` that could cause plots to be mislabeled when using `split.by` and depending on the order of features provided ([#150](https://github.com/samuel-marsh/scCustomize/issues/150)).  
- Fixes issue with automatic point size calculation for Seurat Objects.  
- Added check for presence of dimensionality reduction in `DimPlot_LIGER` ([#153](https://github.com/samuel-marsh/scCustomize/issues/153)).  
- Fixed bug in `Add_Mito_Ribo_LIGER` that caused it to return value of 0 for all cells (Now part of renamed `Add_Mito_Ribo` S3 generic).  
- Fixed legend display is `Clustered_DotPlot` to display percentage instead of proportion to match legend text.  
- Fixed `Percent_Expressing` error when `group_by = "ident"`.  
- Fixed error that caused features in non-default assays to be returned as not found when attempting to plot.  
- Fixed error in `DotPlot_scCustom` that didn't correctly pass `group.by` when plotting ([#158](https://github.com/samuel-marsh/scCustomize/issues/158)).  




# scCustomize 2.0.1 (2023-11-17)  
## Added  
- None.   

  
## Changed  
- Removed warning in `VariableFeaturePlot_scCustom` now fixed in Seurat release.  
   

## Fixes  
- Fixed error in `Add_Mito_Ribo_Seurat` causing failure due to error message when `overwrite = TRUE`.  
- Fixed error in `Add_Top_Gene_Pct_Seurat` to avoid issue that accidentally could call function on normalized data.  
- Fixed error in `Add_Top_Gene_Pct_Seurat` that caused error if more than one counts layer was present.
- Fixed error in `QC_Histogram` that prevented plotting or titling of plots.  



# scCustomize 2.0.0 (2023-11-13)  
## Added  
- Added support for metrics produced by Cell Ranger `multi` pipeline to `Read10X_Metrics` via new parameter `cellranger_multi`.
- Added `dot_size` parameter to `Seq_QC_Plot_*` family of functions.  
- Added two new sequencing QC functions to create and iterate barcode rank plots: `Barcode_Plot` and `Iterate_Barcode_Rank_Plot`.  
- Added `ident_legend` parameter to `QC_Plot_UMIvsGene` to control show/hide of the identity legend ([#121](https://github.com/samuel-marsh/scCustomize/issues/121)).  
- Added support for sparse matrix input in `CellBender_Feature_Diff`.  
- Added `min_count_label` in `CellBender_Diff_Plot` to better control feature labeling.  
- Allow specification of meta data column containing sample names/IDs in `Iterate_DimPlot_bySample` using new `sample_column` parameter.  
- Added new function `MAD_Stats` to calculate to the median absolute deviation of meta.data columns by grouping variable and across entire object.  
- Added new function `Add_Top_Gene_Pct_Seurat` to add another QC measure of cell complexity to object meta.data.  Returns percentage of counts occupied by top XX genes in each cell.  
- Added ability to provide set of custom features to `VariableFeaturePlot_scCustom` using `custom_features` parameter.  
- Added new overall cell QC metric function `Add_Cell_QC_Metrics` to simplify adding cell QC metrics.  Single function call to add Mito/Ribo Percentages, Cell Complexity, Top Gene Percentages, MSigDB Percentages, IEG Percentages, and/or Cell Cycle Scoring (human only).  
- Added 2 new gene lists to package data for use in `Add_Cell_QC_Metrics` function: "msigdb_qc_gene_list" and "ieg_gene_list".
- Added several internal functions to support new MsigDB and IEG capabilities of `Add_Cell_QC_Metrics`.  
- Added new parameters `plot_median` and `plot_boxplot` to `VlnPlot_scCustom` (and `VlnPlot_scCustom`-based plots; e.g., `QC_Plot_*` family) for added visualization.  
- Added `QC_Histogram` to plot QC features (or any feature) using simple histogram.  
- Added `FeatureScatter_scCustom` function to customize Seurat's `FeatureScatter` plots.  
- Added `figure_plot` parameter to all 2D DR (t-SNE, UMAP, etc) based plots ([#127](https://github.com/samuel-marsh/scCustomize/issues/127)).  

  
## Changed  
- Large scale under the hood code adjustments to ensure compatibility with Seurat V5 object structure.  
- Internal code syntax updates independent of Seurat functionality.  
- **HARD DEPRECATION** `Split_FeatureScatter` function has been completely deprecated and it's functionality has been moved to new `FeatureScatter_scCustom`.  
- **SOFT DEPRECATION** The parameter `gene_list` in `Iterate_FeaturePlot_scCustom` and `Iterate_VlnPlot_scCustom` has been soft-deprecated and replaced by `features` parameter.  Specifying `gene_list` will display deprecation warning but continue to function until next major update.  
- The above soft deprecation was to clarify that other features besides genes can be plotted and coincides with update to functions to allow for iterative plots of meta.data or reductions in addition to assay features ([#123](https://github.com/samuel-marsh/scCustomize/issues/123)).  
- Internal rewrite of `Read10X_Metrics` to use new internal helper functions.  
- Changed `Liger_to_Seurat` to transfer the liger_object@H slot in addition to H.norm slot already moved.  
- Replaced `length(x = colnames(x = obj)` with `length(x = Cells(x = obj)` for accurate plotting based on V5 object structure.  
- `Gene_Present` now accepts `assay` parameter.  
- Internal reorganization of some functions within `R/` for better organization.  
- Updated default scCustomize color palettes (`scCustomize_Palette`).  Now if number of colors is greater than 2 but less than 8 the default palette will be `ColorBlind_Pal` (previously it was "polychrome").  Polychrome remains the default when number of colors is between 9-36.  
- Updated parameter default within `scCustomize_Palette` to `ggplot_default_colors = FALSE` to avoid uncessary error when no value supplied.  
- Minimum version of scattermore package updated to v1.2.  
- `DimPlot_scCustom` will now set `label = TRUE` if `label.box` is set to TRUE but `label` is not changed from default.  
- Removed loading of full tidyverse in vignettes to remove from package suggests (lessen dependency installs when not completely needed).  
- Replace Seurat `PackageCheck` (now deprecated), with `rlang::is_installed()` for non-dependency checks.  
- Update vignettes with new features and bug fixes from old code.  
   

## Fixes  
- Fixed issue in `Read10X_Metrics` that caused errors when reading files on windows operating system ([#115](https://github.com/samuel-marsh/scCustomize/issues/115)).  
- Fixed issue in `Create_CellBender_Merged_Seurat` when feature names are changed (underscore to dash) during object creation ([#118](https://github.com/samuel-marsh/scCustomize/issues/118)).  
- Fixed error in `Read10X_h5_Mutli_Directory` when reading Cell Ranger `multi` directories.  
- Added new checks to `VlnPlot_scCustom`, `DimPlot_scCustom`, and `DotPlot_scCustom` to avoid otherwise ambiguous error messages ([#120](https://github.com/samuel-marsh/scCustomize/issues/120)).  
- Fixed internal check message accidentally user facing in `VlnPlot_scCustom` ([#122](https://github.com/samuel-marsh/scCustomize/issues/122)).  
- Fixed cli warning in `Cell_Highlight_Plot` that could cause function to error without proper error message.  
- Fixed handling of file names in `Read_*` functions to avoid unnecessary errors.  
- Replace superseded dplyr syntax/functionality `drop_na(.data[[var]]`, with current dplyr syntax.  
- Internal code fixes to accelerate plotting functions.  
- Fixed default plot colors in `VlnPlot`-based plots when `split.by` is not NULL.  
- Fixed error when trying to plot more than two variables with `group.by` when using `DimPlot_scCustom` ([#128](https://github.com/samuel-marsh/scCustomize/issues/128)).  
- Fixed errors in parameter description for `Add_Mito_Ribo_Seurat` and `Add_Mito_Ribo_LIGER` which incorrectly stated the names of new meta.data/cell.data columns to be added.  
- Fixed bug in `DotPlot_scCustom` that prevented it from working unless `group.by` parameter was explicitly added.  
- Fixed bug in `Case_Check` caused by typo.  
- Fixed color warning messages in `Cluster_Highlight_Plot` and `Meta_Highlight_Plot` that were too verbose.  
- Fixed bug in `Add_Mito_Ribo_Seurat` and `Add_Mito_Ribo_LIGER` which caused error when supplying custom list of features for non-default organism ([#133](https://github.com/samuel-marsh/scCustomize/issues/133)).  
- Fixed bug in `DimPlot_scCustom` preventing that errored when trying to split plot and use `figure_plot` at same time.  



# scCustomize 1.1.3 (2023-07-19)  
## Added  
- None.   
  
  
## Changed  
- None.   
   

## Fixes  
- Fixed manual entry for `Merge_Seurat_List`. 
   


# scCustomize 1.1.2 (2023-07-18)  
## Added  
- Added `aspect_ratio` parameter to all dimensionality reduction plots to control axes ratio of output plot.  
- Added `plot_median` and `median_size` parameters to `QC_Plots_*` functions.  
- Added `split_collect` parameter to `FeaturePlot_scCustom` to collect all guides when using `split.by` for a single feature ([#94](https://github.com/samuel-marsh/scCustomize/issues/94)).  
- Added new parameters to `Clustered_DotPlot` to allow modification of sizes of column text labels, legend text labels, and legend title labels ([#96](https://github.com/samuel-marsh/scCustomize/issues/96)).  
- Added new function `Merge_Sparse_Multimodal_All` for merging multi-modal data (1 matrix per modality) ([#104](https://github.com/samuel-marsh/scCustomize/issues/104)).  
- Added new parameter to `Clustered_DotPlot` named `row_label_fontface` to allow control of fontface used for row labels ([#103](https://github.com/samuel-marsh/scCustomize/issues/103)).  
- Added helper utility `Reduction_Loading_Present`, in part to fix issue with `FeaturePlot_scCustom` and internal feature checking.  
- Added ability to turn off feature/ident clustering in `Clustered_DotPlot` using new parameters: `cluster_feature`, `cluster_ident` ([#106](https://github.com/samuel-marsh/scCustomize/issues/106)).  
- Added `dot_size` parameter to statistics plotting functions `Plot_Cells_per_Sample` and `Plot_Median_*` family.  
- Added new parameter `no_legend` to `Iterate_Meta_Highlight_Plot` to allow for plotting with a plot title instead of plot legend ([#108](https://github.com/samuel-marsh/scCustomize/issues/108)).  
  
  
## Changed  
- Moved `QC_Plots_Feature` to use `VlnPlot_scCustom` under the hood like rest of `QC_Plots_*` functions.  
- Renamed parameter `abort` in `Meta_Present` to `return_none` to align with `Gene_Present` and `Reduction_Loading_Present`.  
- Replace superseded dplyr syntax/functionality `summarise_at`, `select(.data[[var]])`, and `rename(.data[[var]])` with current dplyr syntax.  
- Internal rewrite of plotting sections within `Iterate_Cluster_Highlight_Plot` and `Iterate_Meta_Highlight_Plot` to align with recent updates to base `Cluster_Highlight_Plot` and `Meta_Highlight_Plot` functions.  
   

## Fixes  
- Fixed `QC_Plots_Feature` to respect parameters when passing to `VlnPlot` ([#91](https://github.com/samuel-marsh/scCustomize/issues/91)).  
- Fixed `Read_CellBender_h5_*` functions to support CellBender outputs from STARsolo- or Cell Ranger (pre-V3)-processed data ([#99](https://github.com/samuel-marsh/scCustomize/issues/99)).  
- Fixed `FeaturePlot_scCustom` to allow for plotting of dimensionality reduction loadings ([#97](https://github.com/samuel-marsh/scCustomize/issues/97)).  
- Fixed `Read10X_Multi_Directory` and `Read10X_h5_Multi_Directory` to support files processed with Cell Ranger `multi` pipeline.  
- Fixed bug in `Merge_Seurat_List` that prevented `add.cell.id` from adding correct cell name prefixes ([#113](https://github.com/samuel-marsh/scCustomize/issues/113)).  
   


# scCustomize 1.1.1 (2023-01-13)  
## Added  
- Added `label_color_num` parameter to `PalettePlot` allow control of color labeling.  
- Added ability to rotate x-axis of `Stacked_VlnPlot` 90 degrees or 45 (previously possible) ([#84](https://github.com/samuel-marsh/scCustomize/issues/84)).  
- Added error checks to `Merge_Seurat_List` to avoid ambiguous error messages on failure.  
- Added `Case_Check` checks/messages to all feature-based plotting functions.  
  
## Changed  
- **BREAKING CHANGE** Parameter in `PalettePlot` has been changed from `palette` to `pal`.  
- Updated `PalettePlot` to support `pal` of class "colors".  
- Moved viridis package to Suggests and use paletteer package for viridis palette shortcut functions.  
- Fixed color palette continuity in `Cluster_Highlight_Plot` and `Meta_Highlight_Plot`.  
- `Fetch_Meta` is now S3 generic function that can handle either Seurat or LIGER objects.  
- Rearrange base R code within `R/` scripts for better organization.  
- Completed move of all scCustomize error/warning messages from base R to cli/rlang framework.  
- Move feature checking to internal function.  

## Fixes  
- Fixed potential for column name collision error in `Add_Mito_Ribo_Seurat` and `Add_Mito_Ribo_LIGER`.  
- Fixed `Add_Mito_Ribo_Seurat` to respect provided `mito_name`, `ribo_name` and `mito_ribo_name` values.  
- Updated out-dated documentation for number of package functions.  
- Typo/styling fixes.  
 
 
# scCustomize 1.1.0 (2022-12-22)  
## Added  
- Added `merge` parameter to `Read10X_GEO`, `Read10X_h5_GEO`, `Read_GEO_Delim` and `Read_CellBender_h5_Multi_File`.  
- Added `raster.dpi` parameter to `DimPlot_LIGER`.  
- Added `label` parameter to `FeaturePlot_scCustom` to avoid error collision ([#80](https://github.com/samuel-marsh/scCustomize/issues/80)).  
- Added `vln_linewidth` parameter to control violin outline line width ([#32](https://github.com/samuel-marsh/scCustomize/issues/32)).  
- Added quick meta data getter function `Fetch_Meta` for returning data.frame of object meta data.  
- Added `Extract_Sample_Meta` to extract sample-level meta data from object.  
- Added `Cell_Highlight_Plot` for highlight plots of custom cells not in active ident or meta data.  
- Added `flip` parameter to `Clustered_DotPlot` to enable axes flipping ([#69](https://github.com/samuel-marsh/scCustomize/issues/69)).  

## Changed  
- Updated Imports/Suggests for CRAN compatibility.  
- Under the hood code updates for CRAN compatibility.  
- Rearrange base R code within `R/` scripts for better organization.  

## Fixes  
- Fixed missing documentation for number of package functions.  
- Typo/styling fixes.  

 
# scCustomize 1.0.2 (2022-11-22)  
## Added  
- None. 

## Changed  
- Updated required Seurat version (v4.3.0) to avoid bug in `FindMarkers`.  

## Fixes  
- None.  

 
# scCustomize 1.0.1 (2022-11-10)  
## Added  
- Added `CellBender_Feature_Diff` to return data.frame with count sums and differences between raw and CellBender assays.  
- Added `CellBender_Diff_Plot` to plot differences between raw and CellBender assays using data from `CellBender_Feature_Diff`.  

## Changed  
- **BREAKING CHANGE** Function name changed, `Add_CellBender_Diff` is new name for `Add_Cell_Bender_Diff` in order to unify function names for CellBender related functions.  
- Updated CellBender vignette with new functions.

## Fixes  
- Fixed for automatic color palette selection when only plotting one group.

 
# scCustomize 1.0.0 (2022-10-25)  
## Added
- Added `mito_name` parameter to `QC_Plots_Mito` to allow for custom specification of meta data column name that contains mitochondrial information.
- Added `QC_Plots_Combined_Vln()` function to return patchwork layout of 3 QC plots.
- Added Rhesus Macaque (macaca mulatta) to the accepted species list for `Add_Mito_Ribo_Seurat()` and `Add_Mito_Ribo_LIGER()` ([#28](https://github.com/samuel-marsh/scCustomize/issues/28)).
- Added `alpha_exp` and `alpha_na_exp` parameters to `FeaturePlot_scCustom` to allow for control of color scale transparency ([#21](https://github.com/samuel-marsh/scCustomize/issues/21)).
- `*_Highlight_Plot` functions can now plot multiple variables simultaneously using either one color for all variables or one color per variable ([#34](https://github.com/samuel-marsh/scCustomize/issues/34)).
- Added parameter `figure_plot` to `DimPlot_scCustom()`.  This removes axes and axes labels and adds axis legend on left bottom corner of plot ([#40](https://github.com/samuel-marsh/scCustomize/issues/40)).
- Added parameter `plot_legend` to `Stacked_VlnPlot`.  This solves issue with returning only one shared legend across all features being plotted ([#48](https://github.com/samuel-marsh/scCustomize/issues/48)).
- Added `Add_Cell_Complexity_Seurat` and `Add_Cell_Complexity_LIGER` functions to add cell QC complexity/novelty metric (log10(Genes) / log10(UMIs)).  
- Added `QC_Plots_Complexity` plot for quick plotting of cell complexity score.
- Added 3 new CellBender functions `Read_CellBender_h5_Mat`, `Read_CellBender_h5_Multi_Directory`, `Read_CellBender_h5_Multi_File` to enable easy reading of new CellBender output files.
- Added `raster.dpi` parameter from Seurat to all `DimPlot` `FeaturePlot` or `FeatureScatter` based functions.  
- Added `add.noise` parameter from Seurat to `VlnPlot_scCustom` `Stacked_VlnPlot` functions.  
- Added `group.by` as default listed parameter to added to all`VlnPlot` based `QC_Plot_*`.  
- Added `ensembl_ids` parameter for `Add_Mito_Ribo_*` functions.  If `ensembl_ids = TRUE` functions will retrieve stored ensembl IDs representing mitochondrial and ribosomal genes for accepted default species.  
- Added parameter `label_feature_yaxis` to `FeaturePlot_scCustom`.  Allows for plotting of feature names on secondary y-axis when using `split.by` ([#60](https://github.com/samuel-marsh/scCustomize/issues/60)).  
- Added `Add_Sample_Meta` function for addition of sample-level meta data to cell-level `@meta.data` slot of Seurat objects.  
- Added a matrix check in `Read_GEO_Delim` to check for issues with imported matrices.  Check is modified version of `SeuratObject::CheckMatrix` called `CheckMatrix_scCustom()`.  Will warn if infinite, logical, non-integer (whole), or NA/NaN values are detected in input matrix.  
- `QC_Plot_UMIvsGene` will now returned filtered correlation value that takes into account `meta_gradient_name` if provided in addition to nFeature_RNA and nCount_RNA.
- Added new function `Variable_Features_ALL_LIGER` which allows for detection/selection of variable genes from entire LIGER object instead of iterating by dataset.
- Vignettes/Website updated with new function examples.  

## Changed
- **DEPENDENCY CHANGE** The required version of Seurat has been changed due to errors caused by updates to Matrix package and handling of sparse matrices.  To avoid errors version requirement for Seurat has been updated to 4.2.0. 
- **DEPENDENCY CHANGE** The dittoSeq package has been moved to Suggests to aid package installation. To catch errors a `PackageCheck` warning has been added where needed.  
- **BREAKING CHANGE** Function name for iterative `VlnPlot` has been changed to `Iterate_VlnPlot_scCustom` to reflect that it now uses `VlnPlot_scCustom` to generate plots. 
- `QC_Plot_*` functions now use `VlnPlot_scCustom` internally to unify color scheme and rasterization parameters.
- `*_Highlight_Plot` functions no longer display "Unselected" in plot legend and uses `DimPlot_scCustom` to generate plots ([#34](https://github.com/samuel-marsh/scCustomize/issues/34)).
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
- Fixed bug in `Plot_Density_Custom` when supplying `custom_palette` and multiple features. ([#51](https://github.com/samuel-marsh/scCustomize/issues/53)).
- Fixed bug in `Clustered_DotPlot` so that legend with identities is displayed by factor level of Seurat object idents ([#55](https://github.com/samuel-marsh/scCustomize/issues/55)).
- Fixed bug in `Split_FeatureScatter` to remove test code that prevented function from working properly ([#57](https://github.com/samuel-marsh/scCustomize/issues/57)).
- Fixed bug in `DimPlot_All_Samples`, `Split_FeatureScatter`, and `DimPlot_scCustom` that ignored factor order when plotting groups.
- Fixed error due to deprecation of functions in Matrix package v1.5-0+ ([#61](https://github.com/samuel-marsh/scCustomize/issues/61)).  
- Fixed error that prevent returning `FeaturePlot_scCustom` when setting `split.by` and one or more of features provided was not present in object ([#64](https://github.com/samuel-marsh/scCustomize/issues/64)).  
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
- Fixed `Read_Metrics_10X` errors that occurred due to differing outputs depending on Cell Ranger version or type of assay.
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
