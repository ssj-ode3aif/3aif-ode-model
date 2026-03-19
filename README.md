# 3-AIF: Three-Axis Integrative Framework for Difficult-to-Treat Rheumatic Disease

[![DOI](https://img.shields.io/badge/bioRxiv-preprint-blue)](https://doi.org/XXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository contains the ODE simulation code for the **3-Axis Integrative Framework (3-AIF)**, a computational systems immunology model of difficult-to-treat (D2T) rheumatic disease.

The 3-AIF integrates three biological axes into a six-variable ordinary differential equation (ODE) system:

- **Axis 1**: AID-IgA-Microbiota (Mucosal Tolerance, *T*)
- **Axis 2**: Nerve-Adipose-Immune unit (Danger Signal, *D*)
- **Axis 3**: ISR/ISRmt (Cellular Stress, *S*)

with additional state variables for Energy Reserve (*E*), Recovery Capacity (*R*), and Microbiota Diversity (*M*), gated by mTORC1-AMPK Hill-function metabolic switching.

## Citation

> Jung S. The 3-Axis Integrative Framework (3-AIF) for Difficult-to-Treat Rheumatic Disease: An ODE-Based Systems Model with Transcriptomic Support and an Illustrative 15-Patient Case Series. *bioRxiv* (2026). doi: [PENDING]

## Repository Contents

```
3AIF-ODE-model/
├── README.md                   # This file
├── LICENSE                     # MIT License
├── requirements.txt            # Python dependencies
├── 3aif_simulation.py          # Main ODE system and simulation code
├── figures/                    # Generated figures (after running)
└── supplementary/
    └── default_parameters.csv  # Default parameter values (Supplementary Table S4)
```

## Requirements

- Python ≥ 3.8
- NumPy
- SciPy
- Matplotlib

Install dependencies:

```bash
pip install -r requirements.txt
```

## Usage

### Run all simulations and generate figures

```bash
python 3aif_simulation.py
```

This produces:

| Output | Corresponds to |
|--------|---------------|
| `figures/fig2_trajectories.png` | Figure 2: Three disease trajectory scenarios |
| `figures/fig3a_phase_portrait.png` | Figure 3A: T-D phase portrait |
| `figures/fig3b_metabolic_switch.png` | Figure 3B: mTORC1-AMPK crossover |
| `figures/fig3c_bifurcation.png` | Figure 3C: Bifurcation diagram |
| `figures/fig4a_sensitivity.png` | Figure 4A: Sensitivity analysis tornado plot |

### Use as a module

```python
from 3aif_simulation import run_simulation, DEFAULT_PARAMS

# Run with default parameters
sol = run_simulation(scenario='B')  # Chronic/D2T

# Run with custom parameters
custom = DEFAULT_PARAMS.copy()
custom['alpha2'] *= 1.5  # Increase NAM activation
sol = run_simulation(scenario='B', params=custom)

# Access results
T, D, S, E, R, M = sol.y
time = sol.t
```

## ODE System

The six coupled ODEs (all state variables normalised to [0, 1]):

```
dT/dt = α₁·M·(1−T) − β₁·D·T − γ₁·S·T − ε·env(t)·T
dD/dt = α₂·(1−T)·E + β₂ₑₙᵥ·env(t) − β₂·D + γ₂·S·(1−R)
dS/dt = α₃·D·(1−S) − β₃·R·S + γ₃·mTORC1(E)·(1−S)
dE/dt = −δ₁·D·E + δ₂·R·(1−E) − δ₃·S·E + δ₄·(1−E)
dR/dt = ε₁·AMPK(E)·(1−R) − ε₂·S·R − ε₃·D·R
dM/dt = μ₁·T·(1−M) − μ₂·D·M − μ₃·(1−E)·M
```

where:
- `mTORC1(E) = E² / (Km² + E²)`
- `AMPK(E) = (1−E)² / (Ka² + (1−E)²)`

## Parameter Definitions

See `supplementary/default_parameters.csv` and Supplementary Table S4 in the manuscript for full parameter definitions, biological interpretations, and clinical biomarker proxies.

## Numerical Methods

- Integration: Explicit Runge-Kutta (RK45) with adaptive step-size control
- Solver: `scipy.integrate.solve_ivp`
- Tolerances: rtol = 1e-8, atol = 1e-10

## Notes on Parameter Values

Simulation parameters were iteratively tuned to produce clinically plausible trajectories consistent with observed disease dynamics. They are **not derived from experimental data**. This is a proof-of-concept computational model; the parameter values serve to demonstrate the qualitative dynamics (bistability, hysteresis, attractor trapping) described in the manuscript. See the Limitations section of the paper for full discussion.

## Author

**Sungsoo Jung**, MD, PhD  
Division of Rheumatology, College of Medicine  
SoonChunHyang University Bucheon Hospital  
Republic of Korea

ORCID: [0009-0006-6885-192X](https://orcid.org/0009-0006-6885-192X)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
