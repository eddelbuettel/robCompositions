## Submission of robCompositions 2.5.0

This is a feature release.

### Main changes since 2.4.2 (current CRAN version)

* New function `cellPcaCoDa()` for cellwise-robust principal component analysis
  of compositional data, with `print`/`summary`/`plot` methods, implementing
  Templ (2026, *Advances in Data Analysis and Classification*).
* New function `contaminate_simplex()` to generate cellwise-contaminated
  compositions for stress-testing robust estimators.
* Dependency cleanup: `ggplot2`, `pls` and `data.table` moved from `Depends` to
  `Imports`; `GGally`, `ggfortify` and `kernlab` moved to `Suggests` and used
  conditionally via `requireNamespace()`. `daFisher()` no longer imports all of
  `rrcov` (only `stats::biplot` is now imported as the generic), removing a
  namespace-replacement warning.
* See NEWS for the full list.

## Test environments

* local: macOS 26 (aarch64-apple-darwin20), R 4.5.2 -- 0 errors | 0 warnings | 1 note
* (win-builder R-release and R-devel checked separately)

## R CMD check results

`R CMD check --as-cran` gives **0 ERRORs, 0 WARNINGs**. CRAN incoming feasibility
is OK. The single remaining NOTE is local-only and does not occur on the CRAN
check infrastructure:

* **HTML manual / HTML Tidy** -- "Skipping checking HTML validation: 'tidy'
  doesn't look like recent enough HTML Tidy." This reflects the local macOS
  system `tidy` version only; CRAN's machines have a current HTML Tidy.

(On some local runs an additional transient "unable to verify current time" NOTE
appears when the check machine cannot reach its time server; it is unrelated to
the package.)

Note on URLs: documentation and vignette links to publisher/portal pages that
block automated requests or no longer resolve have been removed or replaced with
DOIs; remaining CRAN links use the canonical https form.

## Downstream dependencies

No reverse dependencies are broken by this release (the public API is only
extended, not changed). Confirmed with a reverse-dependency check at submission
time.
