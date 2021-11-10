# ei
EI: A(n R) Program for Ecological Inference

Authors: Gary King, Margaret Roberts

This program provides a method of inferring individual behavior from aggregate data. It implements the statistical procedures, diagnostics, and graphics from the book [A Solution to the Ecological Inference Problem: Reconstructing Individual Behavior from Aggregate Data](https://gking.harvard.edu/eicamera/kinroot.html) (Princeton: Princeton University Press, 1997), by Gary King.

For more information, see here: https://gking.harvard.edu/eir


# EI Package TODO:
Last Updated: 2021.10.19

## General (Big) Comments:

- [x] Use `cli` for progress reporting
- [ ] Replace workhorse functions in `R/zzz.R` with Rcpp implementations
- [x] Use GitHub actions to auto-check the package on commit
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
- [ ] Rename functions that start with `.` and use `@noRd` tag instead

### Input check

- [ ] ei(): check the `total` argument

## Plotting
Replace base R plots with ggplots

- [x] .tomog -> plot_tomog
- [x] .tomogl -> plot_tomog + option (contour_ML)
- [x] .tomog80CI -> plot_tomog + option (CI)
- [x] .tomog95CI -> plot_tomog + option (CI)
- [x] .tomogE -> plot_tomog + option (points)
- [x] .tomogP2 -> plot_tomog + option (contour_posterior)
- [x] .betabd -> plot_density + option (betab)
- [x] .betawd -> plot_density + option (betaw)
- [x] .xt -> plot_xt
- [x] .xtc -> plot_xtc
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

- [x] bounds1
- [x] ei
- [x] ei.sim
- [x] eiread
- [ ] tomogRxC
- [ ] tomogRxC3d
- [ ] plot


### Data
- [x] rename sample (avoid overlaps with `base::` functions)


### Minor points
- [ ] `message()` instead of `print()` (right now both are used)
