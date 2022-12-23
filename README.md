
# scCustomize <img src="man/figures/scCustomize_Logo.svg" align="right" width="150"/>

[![CRAN
Version](https://img.shields.io/cran/v/scCustomize?color=green&label=CRAN)](https://cran.r-project.org/package=scCustomize)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/scCustomize)](https://cran.r-project.org/package=scCustomize)
[![R-CMD-check](https://github.com/samuel-marsh/scCustomize/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/samuel-marsh/scCustomize/actions/workflows/R-CMD-check.yaml)
[![license](https://img.shields.io/github/license/samuel-marsh/scCustomize)](https://github.com/samuel-marsh/scCustomize/blob/master/LICENSE.md)
[![issues](https://img.shields.io/github/issues/samuel-marsh/scCustomize)](https://github.com/samuel-marsh/scCustomize/issues)
[![DOI](https://zenodo.org/badge/411807769.svg)](https://zenodo.org/badge/latestdoi/411807769)

scCustomize is a collection of functions created and/or curated to aid
in the visualization and analysis of single-cell data using R.

## Vignettes/Tutorials

See accompanying [scCustomize
website](https://samuel-marsh.github.io/scCustomize/) for detailed
tutorials of all aspects of scCustomize functionality.

## Installing scCustomize

scCustomize can be installed from
[CRAN](https://cran.r-project.org/package=scCustomize) on all platforms.
For more detailed instructions see
[Installation](https://samuel-marsh.github.io/scCustomize/articles/Installation.html).

    # Base R
    install.packages("scCustomize")

    # Using pak
    pak::pkg_install("scCustomize")

## Release Notes

A full copy of the changes in each version can be found in the
[NEWS/ChangeLog](https://samuel-marsh.github.io/scCustomize/news/index.html).

**Develop branch**  
I also maintain a separate development branch<sup>\*</sup> that can be
installed by supplying `ref = "develop"` in the devtools or remotes
installation command. Version scheme vX.X.X.9yyy.  
<sup>\*</sup>*Note: While this branch is typically mostly stable it may
contain breaking issues/bugs.*  
I do try and keep [development
ChangeLog](https://github.com/samuel-marsh/scCustomize/blob/develop/NEWS.md)
up to date so it’s easier to follow changes than reading commit history.

## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/samuel-marsh/scCustomize/issues) with details of the issue.

- If possible please include a reproducible example (suggest using
  [SeuratData package](https://github.com/satijalab/seurat-data) pbmc
  dataset for lightweight examples.)

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/samuel-marsh/scCustomize/issues).

- Even if you don’t know how to implement/incorporate with current
  package go ahead a submit!

#### [Pull Requests](https://github.com/samuel-marsh/scCustomize/pulls) are welcome for bug fixes, new features, or enhancements.

- Please set PR to merge with “develop” branch and provide description
  of what the PR contains (referencing existing issue(s) if
  appropriate).
