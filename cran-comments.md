## Resubmission  
This is a resubmission. In this version I have:  

- Single-quoted software names in DESCRIPTION.  
- Omitted redundant "Tools" in package title.  
- Added doi and citation links to description text.
- Fixed instances where TRUE/FALSE were abbreviated.  
- Added return values to documentation of exported functions where missing.  
- Where possible I have eliminated `\dontrun{}` in documentation examples.  
    - In remaining instances there is not lightweight example that can be used in documentation.
    These are documented in more detail in accompanying pkgdown website which utilizes heavier weight example data.  


## R CMD check results

0 errors | 0 warnings | 2 notes

### Test environments  
- Run locally, R4.1.2, Platform: x86_64-apple-darwin17.0 (64-bit) with `devtools:check()`.  
- Also run via GitHub Actions via `usethis::use_github_action_check_standard`
    - macos-latest (release), windows-latest (release), ubuntu-latest (devel), ubuntu-latest (release), ubuntu-latest (oldrel-1).  

## NOTES
1. Maintainer: ‘Samuel Marsh <samuel.marsh@childrens.harvard.edu>’
New submission


2. Imports includes 29 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.  
    - I have worked to move/remove as many IMPORTS to SUGGESTS as possible.  This package aims to simplify a number of different
    visualizations/code tasks in scRNA-seq analysis and as such does have diverse array of dependencies.  I will monitor
    to ensure package functionality.
