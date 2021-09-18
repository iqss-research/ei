# ei
EI: A(n R) Program for Ecological Inference

Authors: Gary King, Margaret Roberts

This program provides a method of inferring individual behavior from aggregate data. It implements the statistical procedures, diagnostics, and graphics from the book [A Solution to the Ecological Inference Problem: Reconstructing Individual Behavior from Aggregate Data](https://gking.harvard.edu/eicamera/kinroot.html) (Princeton: Princeton University Press, 1997), by Gary King.

For more information, see here: https://gking.harvard.edu/eir


# EI Package TODO:
Last Updated: 2021.09.04

## General (Big) Comments:

- [ ] Use `cli` for progress reporting
- [ ] Replace workhorse functions in `R/zzz.R` with Rcpp implementations
- [ ] Use GitHub actions to auto-check the package on commit
- [ ] Add testthat tests
- [ ] Remove imports for single functions where possible

## Specific Comments

### Documentation

- [x] Roxygenize all functions and data
- [ ] Add better data descriptions
- [ ] Ensure return values exist for all functions

### Functions

- [ ] .samp rewrite in Rcpp
- [ ] Add `print.ei` for `"ei"` class objects
- [ ] Allowing users to extract values used in the plot for customization

### Input check

- [ ] ei(): check the `total` argument

## Plotting
Replace base R plots with ggplots

- [ ] .tomog -> plot_tomog
- [ ] .tomogl -> plot_tomogl
- [ ] .tomog80CI -> plot_tomog80CI
- [ ] .tomog95CI -> plot_tomog95CI
- [ ] .tomogE -> plot_tomogE
- [ ] .tomogP2 -> plot_tomogP2
- [ ] .betabd -> plot_betabd
- [ ] .betawd -> plot_betawd
- [ ] .xt -> plot_xt
- [ ] .xtc -> plot_xtc
- [ ] .xtfit -> plot_xtfit
- [ ] .xtfitg -> plot_xtfitg
- [ ] .estsims -> plot_estsims
- [ ] .boundXb  -> plot_boundXb
- [ ] .boundXw -> plot_boundXw
- [ ] .truthfn -> plot_truthfn
- [ ] .bndplot -> plot_bndplot
- [ ] .movieD -> plot_movieD
- [ ] .movie -> plot_movie

### Testing

- [ ] bounds1
- [x] ei
- [ ] ei.sim
- [ ] eiread
- [ ] tomogRxC
- [ ] tomogRxC3d
- [ ] plot


### Data
- [ ] rename sample (avoid overlaps with base:: functions)


### Minor points
- [ ] `message()` instead of `print()` (right now both are used)
