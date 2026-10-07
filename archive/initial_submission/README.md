# Initial submission (archived)

Code and figures of the initial submission, kept for reference. They are superseded by the
revised analysis in the repository root and are not used for any result of the revised
manuscript.

- `base_model/`: base model with the CS risk fixed at 15% (10% and 20% in sensitivity analyses)
- `alternative_model/`: alternative model with a reporting fraction for each of three calendar periods
- `figures/`: Figures 1-3 and Supplementary Figures (`figure2bSup.tiff`, `figure3aSup.tiff`)

The scripts expect `Case.csv` and `CS.csv` (in `data/` of this repository) in the working
directory. The file names passed to `cmdstan_model()` (`syp_BaseModel_ver1.stan`,
`syp_AltModel.stan`) correspond to `BaseModel_final.stan` and `AltModel_final.stan`.
