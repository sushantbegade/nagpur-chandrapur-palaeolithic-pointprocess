# Generative point-process modelling of Palaeolithic site-pattern formation in Nagpur and Chandrapur districts, Central India

[![Code DOI](https://zenodo.org/badge/DOI/<CODE_DOI>.svg)](https://doi.org/<CODE_DOI>)
[![Data DOI](https://zenodo.org/badge/DOI/<DATA_DOI>.svg)](https://doi.org/<DATA_DOI>)

This repository contains the R code for the study

Manuscript Title: Beyond Environmental Suitability: A Statistical-Mechanical Framework for Generative Point-Process
Modelling of Palaeolithic Site-Pattern Formation in Nagpur and Chandrapur Districts, Central India

Author: Sushant Begade

## What the study does
The recorded Palaeolithic site inventory (197 open-air sites) is treated as a spatial point pattern.
Competing generative explanations are compared in a fixed hierarchy: environmental inhomogeneous Poisson
model (M1), M1 plus assumed observation offsets (M2), environment plus a latent log-Gaussian Cox field
(M3), and environment plus a Geyer saturation interaction (M4). The study is a mechanism-identification
study, not a suitability-mapping study. Statistical mechanics is used as a formalism only; no thermodynamic
claim about past populations is made.

## Principal findings (statistical level only)
- The environmental Poisson model (E3 trend) does not reproduce the recorded pattern: observed site pairs at 2-5 km
  are about three times the simulated mean (global envelope test p = 0.031, simplified window; 0.018, exact window).
- Fixed observation offsets reduce the deficiency only slightly (O2 p = 0.0387 at nsim 10000; O3 p = 0.024).
- An environment-plus-latent-field model (LGCP) is not rejected (p = 0.192), so the hierarchy ends at M3.
- A simulation-based recovery study (1,900 synthetic replicates) shows that Gibbs interaction cannot be identified
  at this sample size with the frozen pipeline (power 0.00-0.02 for gamma up to 2.5; recovery of the LGCP mechanism
  0.04-0.37). Gate 6 therefore fails and no interaction claim is made.
- The result is robust to covariate set, grid resolution, quadrature, edge correction, window definition,
  coordinate jitter and redeposited-site exclusion, and not retained in a few variations (Middle Palaeolithic subset,
  a free coefficient on a proxy containing a pattern-informed term, and removal of GIS-map-coordinate sites).
- Close site pairs are concentrated among sites recorded by the same lead author (exploratory analyses). The data
  cannot separate survey and discovery structure, real local clustering and repeated records of one locality.
These are statistical findings about the recorded pattern. They are not evidence for a behavioural mechanism.

## What is not in this repository
Large simulation outputs, derived rasters and site tables are in the Zenodo data record (DOI above).
Third-party source layers (geology, lithology, geomorphology, rivers, DEM, district boundary) are not redistributed;
their sources, dates and checksums are listed in docs/external_download_log.csv, docs/input_manifest.csv and
docs/freeze_v1_manifest.csv.

## How to reproduce
Four levels, in increasing cost:
A. Tables and figures from saved outputs (minutes): download the Zenodo record, place it under the repository root,
   run scripts 64, 65, 66, 67, 68.
B. Empirical tests from the deposited covariate caches (hours): scripts 42-63 (see SCRIPT_INDEX.md).
C. Recovery study (about 14 hours on four workers): scripts 36-41.
D. Rebuild covariates from raw layers: scripts 01-15 (requires the third-party layers listed in the manifests).
Setup:
1. Install R (version in package_versions.csv) and the packages listed there.
2. Edit R/config.R (RAW_DIR, REPO_DIR) and the single CONFIG_PATH line at the top of each script. Paths were replaced
   by /ABSOLUTE/PATH/TO/... in this public copy; no analysis code was changed.
3. Restart R before every long run. Scripts open all inputs before any sink() and write sessionInfo at the end.
4. All random seeds derive from MASTER_SEED (R/freeze_v1_params.R); seed offsets are documented in docs/freeze_v1.md
   and in the deviation log. Results can differ in the last digits across spatstat versions.

## Data, licences and attribution
Code: MIT (LICENSE). Compiled tables, results and figures: CC BY 4.0. Third-party material keeps its own terms:
ESA WorldCover 2021 v200 (CC BY 4.0; (c) ESA WorldCover project / Contains modified Copernicus Sentinel data (2021)
processed by ESA WorldCover consortium; DOI 10.5281/zenodo.7254221), OpenStreetMap roads ((c) OpenStreetMap
contributors, ODbL 1.0), NGDR layers [LICENCE AND SOURCE TO BE CONFIRMED], DEM [PRODUCT AND LICENCE TO BE CONFIRMED],
district boundary [SOURCE TO BE CONFIRMED].

## Known limitations
GIS-map coordinates (60 sites) have an unknown datum; the DEM product is unknown; the covariates are present-day
proxies; observation offsets are assumed strengths, not estimates; the unwinsorised-TPI sensitivity was not run
[DELETE IF RUN]; the recovery study tested Gibbs interaction only up to gamma = 2.5 at one saturation (3).