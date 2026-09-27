# Experimental Data and Model Validation Section Plan

> For a later section (§3 or §4) of the JDSMC paper — exact placement TBD.

This file collects content deferred from the Modelling section that pertains to experimental setup, data characteristics, and model validation.

---

## Source Material (Thesis)

| Thesis Source | Content | Notes |
|---|---|---|
| `Thesis/secs/2-data/1-res_time.tex` | SCR-ASC chamber dimensions table (from Cummins), RTD figure | Chamber geometry |
| `Thesis/secs/2-data/2-data_units.tex` | ppm ↔ mol/m³ conversion, measurement unit table | Unit conversions |
| `Thesis/secs/2-data/0-subs/0-1_test_cell_data.tex` | Test-cell data description (RMC, Hot-FTP cycles) | Data source |
| `Thesis/secs/2-data/0-subs/0-2_truck_data.tex` | Truck data description (4 long-haul trucks) | Data source |
| `Thesis/secs/2-data/3-cross.tex` | NOₓ sensor cross-sensitivity to ammonia | Sensor limitations |
| `Thesis/secs/2-data/4-data_preproc.tex` | Data preprocessing pipeline | Filtering, decimation |
| `Thesis/secs/4-NOx_mdl/3_validation.tex` | Model validation with test-cell data | Validation results |

---

## Planned Content

### SCR-ASC Chamber Dimensions
*(Source: `Thesis/secs/2-data/1-res_time.tex`)*

- Chamber geometry table (from Cummins):
  | Parameter | Value |
  |---|---|
  | SCR length | 9.5 in (24.13 cm) |
  | ASC length | 2 in (5.08 cm) |
  | Chamber diameter | 13 in (33.02 cm) |
  | Total SCR-ASC volume ($V_{scr}$) | 25,014 cm³ |
- These dimensions enter the residence time and the surface-to-volume conversion factor $k_{s2v} = A_{scr}/V$.

---

### Concentration Units: ppm to mol/m³
*(Source: `Thesis/secs/2-data/2-data_units.tex`)*

- Sensors report concentrations in ppm (mole fraction). Conversion to mol/m³ uses ideal gas volume:
  $$x \,[\text{mol/m}^3] = \frac{x \,[\text{ppm}]}{22.4 \times 10^3 \times \frac{273.15 + T}{273.15}} \times 10^{-3}$$
- This temperature-dependent conversion is applied during data preprocessing; within the model, all concentrations are in mol/m³.
- Measurement unit table:
  | Property | Unit |
  |---|---|
  | Concentration | $10^{-3}$ mol/m³ |
  | Temperature | $\times 10 + 200\,°C$ |
  | Flow Rate | $\times 10$ g/s |
  | Urea Injection | $\times 10^{-1}$ ml/s |
  | Time | s |

---

### Additional Content (TBD)

- Test-cell and truck data descriptions
- NOₓ sensor cross-sensitivity discussion
- Data preprocessing and filtering
- Model validation results (FTIR vs production sensor comparison)
- Aging detection results (Wald-test)
