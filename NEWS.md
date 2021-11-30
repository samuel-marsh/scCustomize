## scCustomize 0.6.X (2021-XX-XX)
### Added
- Addition 1

### Changes
- Change 1

### Fixes
- Fixed barcode name duplication checks in `Merge_Sparse_Data_All`. ([#8](https://github.com/samuel-marsh/scCustomize/issues/8))
- Fixed package imports in DESCRIPTION to avoid installation errors.
- Fixed NULL check in `Read_Metrics_10X`.
 

## scCustomize 0.6.1 (2021-11-19)
### Added
- Added plot spacing control to `StackedVlnPlot` with parameters `plot_spacing` and `spacing_unit`. ([#6](https://github.com/samuel-marsh/scCustomize/issues/6))
- Added `scCustomize_Palette` function select palette to use (simplify internal code).

### Changes
- Changed citation info to reflect global DOI and not version DOI.

### Fixes
- Restore package color palette defaults to `Iterate_VlnPlot`.  
- Fix `Iterate_...` function checks for file path parameter if `file_path = NULL`.
  
## scCustomize 0.6.0 (2021-11-16)
### Added
- scCustomize is public!!  Version 0.6.0 is released!

### Changes
- Many function names have changed since private release see reference page/manual for updated function names.
