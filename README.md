# Prostate Cancer Agent-Based Model (PhysiCell)

Agent-based model of prostate tumor evolution under androgen deprivation therapy (ADT), simulating competition between PTEN-normal (androgen-sensitive, S) and PTEN-deleted (androgen-resistant, R) cells. Built using PhysiCell v1.14.0.

---

## Cell Types

| Cell Type | Description |
|---|---|
| PTEN_normal (S) | Androgen-sensitive cells; growth depends on testosterone via Michaelis-Menten kinetics |
| PTEN_deleted (R) | Androgen-resistant cells; reduced testosterone dependence |

Both cell types share identical mechanical parameters (adhesion, repulsion, motility) — differentiation arises through cohort-specific growth rate parameters.

---

## Microenvironment Substrates

| Substrate | Diffusion coefficient | Decay rate | Boundary condition |
|---|---|---|---|
| Oxygen | 100 µm²/min | 0 1/min | 38 mmHg (Dirichlet, all faces) |
| Testosterone | 1 µm²/min | 0 1/min | 8 ng/ml (Dirichlet, all faces) |

---

## Clinical Cohorts

Simulations are parameterized by three clinical cohorts, each with fitted growth rate parameters for S and R cells:

| Cohort | Description | cohort value in XML |
|---|---|---|
| BR | Biopsy Resistant | 0 |
| TR | Treatment Resistant | 1 |
| CTR | Control | 2 |

Growth rates, apoptosis rates, and testosterone sensitivity (Michaelis-Menten parameters m, n, p) differ across cohorts. See `runs/ABMruns_masterlist_prostatecancer.csv` for all parameter values.

---

## Two Simulation Conditions

The two directories differ **only** in mechanical spatial boundary, allowing isolation of the effect of spatial confinement on tumor evolution:

| Condition | Directory | Mechanical boundary | Diffusion domain |
|---|---|---|---|
| Unconstrained growth | `ABM_unconstrained/` | ±250 µm (mesh bounding box) | ±250 µm |
| Spatially confined | `ABM_densepacking/` | ±125 µm (hardcoded in `core/PhysiCell_standard_models.cpp`) | ±250 µm |

In the confined condition, cells are mechanically restricted to a 250×250 µm region while substrate diffusion continues over the full 500×500 µm domain, mimicking tumor growth within a bounded ductal/acinar compartment.

---

## How to Reproduce a Simulation

### 1. Get PhysiCell v1.14.0
```bash
git clone https://github.com/MathCancer/PhysiCell.git
cd PhysiCell
git checkout v1.14.0
```

### 2. Copy files from this repo into your PhysiCell directory
For unconstrained runs:
```bash
cp ABM_unconstrained/config/PhysiCell_settings.xml PhysiCell/config/
cp ABM_unconstrained/config/cell_rules.csv PhysiCell/config/
cp ABM_unconstrained/custom_modules/custom.cpp PhysiCell/custom_modules/
cp ABM_unconstrained/custom_modules/custom.h PhysiCell/custom_modules/
cp ABM_unconstrained/core/PhysiCell_standard_models.cpp PhysiCell/core/
```
For spatially confined runs, use `ABM_densepacking/` instead.

### 3. Compile
```bash
cd PhysiCell
make
```

### 4. Configure a run from the masterlist
```bash
python3 runs/ABMruns_updatexml.py PhysiCell/config/PhysiCell_settings.xml runs/ABMruns_masterlist_prostatecancer.csv RUN_ID
```
Replace `RUN_ID` with any integer from `ABMruns_masterlist_prostatecancer.csv`.

### 5. Run
```bash
cd PhysiCell
./project
```
Or submit via SLURM using the scripts in `runs/`.

---

## Repository Structure

```
ABM_unconstrained/                        # Unconstrained simulation files
  config/
    PhysiCell_settings.xml
    cell_rules.csv
  custom_modules/
    custom.cpp
    custom.h
  core/
    PhysiCell_standard_models.cpp

ABM_densepacking/                         # Spatially confined simulation files
  config/
    PhysiCell_settings.xml
    cell_rules.csv
  custom_modules/
    custom.cpp
    custom.h
  core/
    PhysiCell_standard_models.cpp         # double boundary = 125

runs/
  ABMruns_masterlist_prostatecancer.csv   # Full parameter table for all runs

analysis/
  ABMruns_PCa_dataanalysis.py             # Main data analysis script

README.md
```
