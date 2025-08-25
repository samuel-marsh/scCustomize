## Hotfix to v3.1.0 and v3.1.1  
Fixed fuerther errors of doi/url formatting in manual docs.  Sorry for missing them in earlier submission!


## Update 
This is major update to version 3.1.0.  Contains bug fixes and new functions. This version also updates maintainer email address and adds one new author to package.  See News.md for full changelog.  



## R CMD check results

0 errors | 0 warnings | 1 notes

### Test environments  
- Run locally, R4.5.0, Platform: aarch64-apple-darwin20 with `devtools:check()`.  
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

## Other Notes
1. GitHub Actions check returning strange errors only on macos (release) and only sometimes).  
    - The errors are from failures running package examples.  This includes functions that have been part of prior CRAN releases.
    NO errors are found when checking locally on macos platform using R 4.4.0 and none are found in GitHub Actions check on linux
    or windows platforms.  I believe to be error in GitHub Actions workflow and I have therefore refrained from adding `dontrun`
    to examples that run fine on other platforms.  This was also the case with v2.0.0, v2.0.1, 2.1.0, 2.1.2 which passed macos checks on CRAN,
    furthering it is likely a GitHub Actions issue.
