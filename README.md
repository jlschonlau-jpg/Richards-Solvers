========================================
1D Richards Equation Solver (Python)
Author: Jacob Schonlau
Institution: Brigham Young University (Research)
Date: April 2026
========================================

Overview
--------
This repository contains a 1D numerical solver for
the mixed-head form of the Richards Equation. It is designed to simulate
unsaturated water flow through various soil textures.

The solver handles the highly non-linear Van Genuchten soil functions using
adaptive time-stepping and Picard iteration.

Governing Equation:
    C(h) * ∂h/∂t = ∂/∂z [ K(h) * ( ∂h/∂z - 1 ) ]

Project Structure
-----------------
richards_solver.py       -Core execution script and numerical engine.
params.json               -Configuration file for soil physics, grid setup,
                           solver settings, and boundary conditions.
requirements.txt          -Minimal dependencies (NumPy, SciPy, Matplotlib).
simulation_report.txt     -Auto-generated performance log summarizing runtime,
                           iteration counts, and Mass Balance Error (MBE).


Getting Started
---------------
1. Install Dependencies
       pip install -r requirements.txt

2. Running a Simulation
       python richards_solver.py --config params.json


Configuration Guide (params.json)
---------------------------------
The simulation is controlled entirely through the JSON configuration file.

Soil Physics:
    Uses the Van Genuchten model.
    Provide: alpha, n, theta_r, theta_s, Ks.

Grid Setup:
    total_depth     Column length (cm)
    node_spacing    Vertical discretization (dz)

Solver Settings:
    averaging_mode  "arithmetic", "geometric", "harmonic", or "upwind"
                    "all" → runs all four modes sequentially and compares
                    their performance in a single plot.
                    (Note in simple cases these will be so close you may not see the difference)

                    Upwind is usually the most stable of all of these.

    initial_dt      Starting time step.
                    Note: Time units are determined by the flux and Ks rates
                    (e.g., if Ks is cm/hr, max_time is in hours)


Numerical Details
-----------------
Spatial Discretization:
    Central finite difference on a staggered grid.
    z is defined as positive downwards (depth).

Temporal Discretization:
    Fully implicit (Backward Euler).

Linear Algebra:
    O(N) tridiagonal matrix inversion using scipy.linalg.solve_banded.

Adaptive Time-Stepping:
    Automatically scales dt based on Picard iteration counts to maintain
    stability during sharp wetting fronts.

Boundary Conditions:
    Top:    Neumann (constant flux / rainfall)
    Bottom: Dirichlet (constant pressure / water table)


Analyzing Results
-----------------
Upon completion, the solver generates a Matplotlib plot showing the final
moisture profile (theta).

Technical Audit:
    Check simulation_report.txt after each run.
    A successful simulation should have a Mass Balance Error (MBE)
    near zero (< 1e-7 cm).

If the MBE is high or the solver crashes:
    • Reduce node_spacing
    • Reduce initial_dt

The visualization engine employs Auto-Zooming Axes. Rather than using static bounds,
the plot window wraps tightly around the current theta values. This ensures that even the most subtle
capillary wicking or evaporation-driven drying is visible to the researcher,
even when those changes occur at the 4th or 5th decimal place.


Handover & Maintenance Notes
----------------------------
This codebase was designed to be modular for future researchers.

To Add New Soil Models:
    Extend the VanGenuchtenSoil class.

To Add Sink Terms:
    Modify the residual calculation in perform_timestep to include
    root water uptake or evaporation.

Stability Notes:
    The solver is optimized for sandy loam.
    For heavy clay scenarios, use a very small initial_dt (≈ 1e-5)
    to handle steep suction gradients.

Note:
    This solver was developed with collaborative AI assistance to ensure
    code readability and robust error handling. Users should verify
    research outputs against analytical or benchmark solutions.
  
