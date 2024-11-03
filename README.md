# Gammapy field-of-view (FOV) analysis using 3D acceptance map

These codes were developed for VERITAS analysis of 1LHAASO J0622+3754, an extended (r86 = 1 deg) very-high-energy gamma-ray source. These codes allow one to generate a 3D acceptance map and perform mimic data analysis (and bias correction) using Gammapy. Part of the code accesses the VERITAS database to fetch run information. Errors/suggestions should be directed to **Jooyun Woo (jw3855@columbia.edu)**.

# Environment variables to be set

- VTS_GP: path to these codes
- DB_HOST: VERITAS database hostname
- DB_USER: VERITAS database username
