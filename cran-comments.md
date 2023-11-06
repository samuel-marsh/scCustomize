## Minor Version Update 
This is a major version update to v2.0.0. In this version I have:  

- Added a number of new functions, added new function parameters, and fixed bugs (see News.md).  
- Ensured compatibility with major version of Seurat package.  
- Fixed in example code causing current CRAN check errors with current package version (v1.1.3).


## R CMD check results

0 errors | 0 warnings | 1 notes

### Test environments  
- Run locally, R4.3.1, Platform: x86_64-apple-darwin17.0 (64-bit) with `devtools:check()`.  
- Also run via GitHub Actions via `usethis::use_github_action_check_standard`
    - macos-latest (release), windows-latest (release), ubuntu-latest (devel), ubuntu-latest (release), ubuntu-latest (oldrel-1).  

## NOTES
1. Imports includes 28 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.  
    - I have worked to move/remove as many IMPORTS to SUGGESTS as possible.  This package aims to simplify a number of different
    visualizations/code tasks in scRNA-seq analysis and as such does have diverse array of dependencies.  I will monitor
    to ensure package functionality.  
2. Suggests or Enhances not in mainstream repositories:
     rliger  
     - I hope this issue to be fixed soon as the maintainers have assigned team member to fix issue
     (https://github.com/welch-lab/liger/issues/293).  However, in the interim this loss will not effect package global package
     functionality.  Those looking for rliger functionality can still easily download versions from CRAN archive and GitHub.  
