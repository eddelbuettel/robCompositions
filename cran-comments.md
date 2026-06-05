## Submission of robCompositions 2.6.0

This is a feature release. Version 2.5.0 was prepared but never submitted to
CRAN; its changes are rolled into this 2.6.0 submission, so the list below
covers all changes since 2.4.2 (the current CRAN version).

### Main changes since 2.4.2 (current CRAN version)

* New function `cellNetCoDa()` for cellwise-robust network (Gaussian graphical
  model) estimation of compositional data, with `print`/`summary`/`plot`
  methods. Cellwise outliers are detected and down-weighted before a sparse
  graph is fitted by neighbourhood selection (`"mb"`) or the graphical lasso on
  the centred-logratio covariance (`"glasso"`); the penalty is chosen by StARS.
  The two back-end packages `glmnet` and `glasso` are in `Suggests` and used
  conditionally via `requireNamespace()` (in the function, examples and tests).
* New function `cellPcaCoDa()` for cellwise-robust principal component analysis
  of compositional data, with `print`/`summary`/`plot` methods, implementing
  Templ (2026, *Advances in Data Analysis and Classification*).
* New function `contaminate_simplex()` to generate cellwise-contaminated
  compositions for stress-testing robust estimators.
* `imputeBDLs()` has been renamed to `imputeBDL()`. `imputeBDLs()` and
  `impRZilr()` are kept as deprecated aliases (via `.Deprecated()`) that forward
  to it, so existing code keeps working. A latent bug was fixed so that values
  below the detection limit supplied as small positive numbers (rather than
  exact zeros) are now imputed; results for the documented usage (rounded zeros
  pre-coded as 0) are unchanged.
* Dependency cleanup: `ggplot2`, `pls` and `data.table` moved from `Depends` to
  `Imports`.
* See NEWS for the full list.

## Test environments

* local: macOS 26 (aarch64-apple-darwin20), R 4.5.2 -- 0 errors | 0 warnings | 1 note
* win-builder R-release and R-devel (checked separately)

## R CMD check results

`R CMD check --as-cran` gives **0 ERRORs, 0 WARNINGs**. The single remaining NOTE
is local-only and does not occur on the CRAN check infrastructure:

* **HTML manual / HTML Tidy** -- "Skipping checking HTML validation: 'tidy'
  doesn't look like recent enough HTML Tidy." This reflects the local macOS
  system `tidy` version only; CRAN's machines have a current HTML Tidy.

(On some local runs an additional transient "unable to verify current time" NOTE
appears when the check machine cannot reach its time server; it is unrelated to
the package.)

Note on URLs: documentation and vignette links use the canonical https form or
DOIs.

## Downstream dependencies

No reverse dependencies are broken by this release. The public API is only
extended (new functions) or extended compatibly: `imputeBDLs()` and `impRZilr()`
remain available as deprecated aliases with unchanged results for their
documented usage. Confirmed with a reverse-dependency check at submission time.
