# scCustomize <img src="https://github.com/samuel-marsh/scCustomize/blob/master/data/scCustomize_Logo.png?raw=true" alt="drawing" width="150" align="right"/>  

![GitHub all releases](https://img.shields.io/github/downloads/samuel-marsh/scCustomize/total?style=flat-square)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/samuel-marsh/scCustomize?style=flat-square)
![version](https://img.shields.io/badge/version-v0.5.0-success?style=flat-square)
[![license](https://img.shields.io/github/license/samuel-marsh/scCustomize?style=flat-square)](https://github.com/samuel-marsh/scCustomize/blob/master/LICENSE)
[![last_commit](https://img.shields.io/github/last-commit/samuel-marsh/scCustomize?style=flat-square)](https://github.com/samuel-marsh/scCustomize/commits) [![issues](https://img.shields.io/github/issues/samuel-marsh/scCustomize?style=flat-square)](https://github.com/samuel-marsh/scCustomize/issues)
[![commit_freq](https://img.shields.io/github/commit-activity/m/samuel-marsh/scCustomize?style=flat-square)](https://github.com/samuel-marsh/scCustomize/commits)  


scCustomize is a collection of functions created and/or curated to aid in the analysis of single-cell data using R.

### About/Goals
The goals of scCustomize are to make things "easier", faster, and more accessible for many common tasks in scRNA-seq analysis.  
- Customized versions of many commonly used plotting functions (and some custom ones) to allow greater flexibility in visualization and more aesthetic visuals.
- Easy iterative functionality.  Many plotting and read/write functions can be easily automated with loops, apply, purrr etc.  However, these can be intimidating to novice user and often can be made easier through wrapping into a function.
- Replace multiple lines of code with single function.  There are often aspects of analysis that are run everytime and sometimes multiple times within an analysis that take multiple lines of code.  scCustomize takes some of these and provides single function to streamline the process.

Currently the package is primarily centered around interactivity with Seurat Objects.  If users are interested in adapting functions (or creating separate functions) to provide comparable use with SCE or other object formats I would be happy to add them.  See below for more info on PRs.


## Installing scCustomize
scCustomize can be installed from GitHub using either devtools or remotes package:
```
devtools::install_github(repo = "samuel-marsh/scCustomize")

remotes::install_github(repo = "samuel-marsh/scCustomize")
```
If you had previously installed [colorway](https://github.com/hypercompetent/colorway) package prior to installing scCustomize please make sure to update to >= v0.2.0.


## Branches  
### Master branch
I will push full releases with version scheme vX.X.X.  
See [News.md](https://github.com/samuel-marsh/scCustomize/blob/master/News.md) file for ChangeLog with additions, changes, and fixes contained in each release.


### develop branch
I also maintain a separate development branch<sup>\*</sup> that can be installed by supplying `ref = "develop"` in the devtools or remotes installation command.  Version scheme vX.X.X.yyy.  

  - <sup>\*</sup>*Note: While this branch is typically mostly stable it may contain breaking issues/bugs.*  
  - I do try and keep [development ChangeLog](https://github.com/samuel-marsh/scCustomize/blob/develop/News.md) up to date so it's easier to follow changes than reading commit history.
  

## Vignettes/Tutorials  
See accompanying pkgdown website for detailed tutorials of all aspects of scCustomize functionality.  [Website Link](coming_soon).

## Bug Reports/New Features
#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/samuel-marsh/scCustomize/issues) with details of the issue.
- If possible please include a reproducible example (suggest using [SeuratData package](https://github.com/satijalab/seurat-data) pbmc dataset for lightweight examples.)

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/samuel-marsh/scCustomize/issues).
- Even if you don't know how to implement/incorporate with current package go ahead a submit!
  
#### [Pull Requests](https://github.com/samuel-marsh/scCustomize/pulls) are welcome for bug fixes, new features, or enhancements.
- Please set PR to merge with "develop" branch and provide description of what the PR contains (referencing existing issue(s) if appropriate).
