# Ternary_lipid_FH_parameter

# GROMACS Simulation Setup -- all atom molecular dynamcis (AA-MD)

This directory contains all the input files necessary to run molecular dynamics (MD) simulations of lipid bilayer systems using **GROMACS**.  
The workflow follows a standard protocol: topology preparation, energy minimization, multi-step equilibration, and final production runs.

---

## File Organization

### Topology and Force Field Files
- **forcefield.itp** – General force field parameters.
- **DPPC.itp** – Parameters for DPPC lipid.
- **DUPC.itp** – Parameters for DUPC lipid.
- **CLA.itp** – Parameters for chloride ion.
- **SOD.itp** – Parameters for sodium ion.
- **TIP3.itp** – Parameters for TIP3P water model.

These `.itp` files define the bonded and non-bonded parameters for each molecule used in the simulation.

---

### MD Parameter Files (`.mdp`)
- **mdout.mdp** – Template file with general MD settings.
- **step6.0_minimization.mdp** – Energy minimization.
- **step6.1_equilibration.mdp** – Initial equilibration (short, with strong restraints).
- **step6.2_equilibration.mdp** – Continued equilibration with adjusted restraints.
- **step6.3_equilibration.mdp**
- **step6.4_equilibration.mdp**
- **step6.5_equilibration.mdp**
- **step6.6_equilibration.mdp**

  The equilibration steps gradually relax positional restraints, temperature, and pressure couplings to bring the system to stable conditions.

- **step7_production.mdp** – Final production MD parameters (unrestrained trajectory for analysis).

---

## Typical Workflow

1. **System preparation**  
   Generate coordinates using CHARMM-GUI.

2. **Energy Minimization**  
   ```bash
   gmx grompp -f step6.0_minimization.mdp -c initial.gro -p topol.top -o em.tpr
   gmx mdrun -deffnm em

   
# Martini 3.0 MD Simulation Setup -- CG-MD

This directory contains the input files required to run **coarse-grained molecular dynamics (CG-MD)** simulations with the **Martini 3.0 force field** in **GROMACS**.  
The workflow consists of energy minimization, equilibration, and production runs.

---

## File Organization

### Force Field and Topology Files
- **martini_v3.0.0.itp** – Core Martini 3.0 parameters.
- **martini_v3.0.0_ffbonded_v2_openbeta.itp** – Bonded interaction definitions (open beta).
- **martini_v3.0.0_ions_v1.itp** – Ion parameters.
- **martini_v3.0.0_phospholipids_PC_v2_openbeta.itp** – Phosphatidylcholine (PC) lipid parameters (open beta version).
- **martini_v3.0.0_solvents_v1.itp** – Solvent parameters (e.g., water, alcohols).
- **martini_v3.0.0_sterols_v1.itp** – Sterol parameters (e.g., cholesterol).

These files should be included in the `topol.top` file to build the full system topology.

---

### MD Parameter Files (`.mdp`)
- **1_minimization.mdp** – Energy minimization to remove steric clashes and relax the system.
- **2_equilibration.mdp** – Equilibration phase under controlled temperature and pressure with possible restraints.
- **3_production.mdp** – Main production run (default simulation temperature).
- **3_production_300.mdp** – Production run specifically at 300 K.

- **mdout.mdp** – Generic MD parameter output file.

---

## Workflow

1. **Energy Minimization**
   ```bash
   gmx grompp -f 1_minimization.mdp -c system.gro -p topol.top -o em.tpr
   gmx mdrun -deffnm em

2. **Equilibration**
   ```bash
   gmx grompp -f 2_equilibration.mdp -c em.gro -p topol.top -o eq.tpr
   gmx mdrun -deffnm eq
3. **Production Run**
      ```bash
   gmx grompp -f 3_production.mdp -c eq.gro -p topol.top -o md.tpr
   gmx mdrun -deffnm md


# χ–Parameter Extraction (RDF → PMF → ω → χ)

This folder documents the **concept and steps** for extracting Flory–Huggins interaction parameters \(\chi_{ij}\) from molecular dynamics data.

**Note:** **Home-made Python scripts** (based on **MDAnalysis** and **mdtraj**) were used a posteriori to analyze trajectories and generate plots.  
**Data availability:** To maximize transparency and reproducibility (and to support the manuscript), we **share the precomputed RDF files** (`.dat`) for all lipid–lipid pairs analyzed (see **Provided RDF data** below).

---

## What you need (inputs)

- Preprocessed, imaged, and centered MD trajectories (AA or CG)
- Topology files compatible with MDAnalysis/mdtraj
- Temperature \(T\) (or \(k_BT\)) used for each dataset
- **Pair distribution functions \(g_{ij}(r)\) (shared as `.dat`)** — you can use these directly, or regenerate RDFs with your own tools following the steps here

---

## What you get (outputs)

- \(W_{ij}(r)\): Potential of mean force (PMF) curves shifted to 0 at large \(r\)
- \(\omega_{ij}\): **Effective contact energies** from the first PMF minimum (self-terms \(\omega_{ii}=0\))
- \(\chi'_{ij}\): **FH parameters** in \(k_BT\) units
- \(\chi_{ij}\): **Lattice-scaled** FH parameters (optional; \(\chi_{ij}=z_{ij}\chi'_{ij}\))
- Summary tables and heatmaps suitable for figures and downstream FH/MC calculations

---

## Provided RDF data (to enhance reproducibility)

- **Location:** `CG_MD/all_beadpairs/`
- **Format:** ASCII `.dat` files, one per pair/state point. Default columns:


---

# Order Parameter — `S_mix` (Leaflet Mixing Entropy)

Quantifies **leaflet demixing** between two lipid species (e.g., **DPPC** vs **DIPC**) using a **Voronoi neighbor graph**.  
Cholesterol (CHOL) shapes local geometry but is **not** counted in the like/unlike fractions.

---

## Definition

\[
S_{\text{mix}} = -\big[\,X_{\mathrm{SL}}\log_2 X_{\mathrm{SL}} + X_{\mathrm{DL}}\log_2 X_{\mathrm{DL}}\,\big] \in [0,1]\ \text{bits}
\]

- \(X_{\mathrm{SL}}\): fraction of **similar** edges (A–A + B–B)  
- \(X_{\mathrm{DL}}\): fraction of **dissimilar** edges (A–B)  
- **Interpretation:**  
  - \(S_{\text{mix}} \approx 1\) → random mixing  
  - \(S_{\text{mix}} \to 0\) → strong demixing

> Units are **bits** (base-2 logarithm). Run **upper** and **lower** leaflets separately.

---

## Run

```bash
python Smix_calc.py \
  -s system.pdb \
  -f md.xtc \
  -leaflet upper \
  -skip 10 \
  -nt 8 \
  -protein NO




