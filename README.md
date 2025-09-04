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
