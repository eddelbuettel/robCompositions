## Submission of robCompositions 2.5.0

This is a feature release.

### Main changes since 2.4.2 (current CRAN version)

* New function `cellPcaCoDa()` for cellwise-robust principal component analysis
  of compositional data, with `print`/`summary`/`plot` methods, implementing
  Templ (2026, *Advances in Data Analysis and Classification*).
* New function `contaminate_simplex()` to generate cellwise-contaminated
  compositions for stress-testing robust estimators.
* See NEWS for the full list.

## Test environments

* local: macOS 26 (aarch64-apple-darwin20), R 4.5.2 -- 0 errors | 0 warnings | 0 notes
  (other than the two purely local notes described below)
* (please also run win-builder release/devel and R-hub before final submission)

## R CMD check results

`R CMD check --as-cran` gives **0 ERRORs, 0 WARNINGs**. CRAN incoming feasibility
is OK. The only two remaining NOTEs are artefacts of the local check machine and
do **not** occur on the CRAN check infrastructure:

1. **"unable to verify current time"** -- the local checker could not reach its
   online time server; this is transient/network and not related to the package
   (the same check returned OK on other local runs).

2. **HTML manual / HTML Tidy** -- "Skipping checking HTML validation: 'tidy'
   doesn't look like recent enough HTML Tidy." This reflects the local macOS
   system `tidy` version only; CRAN's machines have a current HTML Tidy.

Note on URLs: documentation links to publisher/portal pages that block automated
requests have been removed in favour of DOIs (`\doi{}`) where available; the OECD
data-portal references (which return HTTP 403 to any non-browser client) are now
given as plain-text attribution rather than checked links.

## Downstream dependencies

No reverse dependencies are broken by this release (the public API is only
extended, not changed). Please confirm with a reverse-dependency check at
submission time.
