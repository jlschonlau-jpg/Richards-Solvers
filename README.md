## 1D Richards Equation Solver (Python)

### Overview
This solver provides a numerical engine for the **mixed-head form** of the Richards Equation.

The solver uses **Picard Iteration** and **adaptive time-stepping** to maintain stability during extreme events (like surface evaporation from wet clay or rapid infiltration into dry sand)

#### Governing Equation
The solver works with the vertical moisture movement defined by:

C(h) * ∂h/∂t = ∂/∂z [ K(h) * ( ∂h/∂z - 1 ) ]

Where:
* **h**: Pressure head
* **C(h)**: Specific moisture capacity (dθ/dh)
* **K(h)**: Hydraulic conductivity

### Getting Started

    python richards_solver.py --config params.json

### Configuration Guide (`params.json`)

The simulation parameters isn't in the "richards_solver.py", it is controlled via the JSON configuration file. If you don't use --config, the program will by default use "params.json". If you do, you can specify any scenario.

| Category | Parameter | Description |
| :--- | :--- | :--- |
| **Soil Properties** | `alpha`, `n` | Van Genuchten shape parameters. |
| | `theta_r`, `theta_s` | Residual and Saturated volumetric water content. |
| | `Ks` | Saturated hydraulic conductivity. |
| **Grid Setup** | `total_depth` | Total column length (e.g., cm). |
| | `node_spacing` | Vertical discretization (dz). |
| **Solver Settings**| `averaging_mode` | `arithmetic`, `geometric`, `harmonic`, `upwind`, or `all`. |
| | `initial_dt` | Starting time step (automatically adjusted during run). |
| | `max_time` | Total simulation duration. |

> **Note on Averages:** The flux is calculated in between each node, which requires the hydraulic conductivity (Ks). This must me some kind of average between the Ks value of each node. I implimented four different possible averaging modes, with the additional setting "all" which runs them all. `upwind` is generally the most stable averaging mode, but not always. I don't understand the nuance of when different modes might be better so I coded them all.
---

### Numerical Implementation

* **Spatial Discretization:** Central finite difference on a staggered grid. Depth (z) is defined as positive downwards.
* **Temporal Discretization:** Fully implicit (Backward Euler) to ensure stability at larger time steps.
* **Linear Algebra:** Employs an O(N) tridiagonal matrix inversion using `scipy.linalg.solve_banded`.
* **Adaptive Time-Stepping:** The solver monitors Picard convergence. It "chops" dt upon failure and expands dt when convergence is rapid (≤ 5 iterations).

---

### Analyzing Results

#### Visualization
The solver generates a profile plot of θ (Volumetric Water Content) vs. Depth. The axes utilize **Auto-Zooming**, wrapping tightly around current values to reveal subtle capillary wicking or evaporation-driven drying that might otherwise be missed. Because it zooms in, you sometimes it will appear as if the theta levels vary quite a bit but it turns out that the scale is just super zoomed in.

#### Technical Audit
After every run, check `simulation_report.txt`. 
* **Success Criteria:** A robust simulation should maintain a **Mass Balance Error (MBE)** near zero (< 1e-7 cm). 
* **Troubleshooting:** If the MBE is high or the solver fails, decrease `node_spacing` or the `initial_dt` in the config file.

---

### Citation & Acknowledgements
This solver was developed at **Brigham Young University** with collaborative AI assistance. This program has NOT been verified against benchmark solutions (e.g., HYDRUS-1D) 
