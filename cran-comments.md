## R CMD check results

0 errors | 0 warnings | 2 notes

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
