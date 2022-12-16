## R CMD check results

0 errors | 0 warnings | 1 notes

## NOTE
- Imports includes 28 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.  
    - I have worked to move/remove as many as possible.  This package aims to simplify a number of different
    code tasks in scRNA-seq analysis and as such does have diverse array of dependecies.  I will monitor
    to ensure package functionality.
