import os
import time
import json
import argparse
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
from typing import Tuple, Dict, Any, List

# ==========================================
# 1. SOIL PHYSICS MODULE
# ==========================================
class VanGenuchtenSoil:
    def __init__(self, theta_r: float, theta_s: float, alpha: float, n: float, Ks: float):
        self.theta_r = theta_r
        self.theta_s = theta_s
        self.alpha = alpha
        self.n = n
        self.m = 1.0 - (1.0 / n)
        self.Ks = Ks

    def calc_theta(self, h: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Calculates volumetric water content."""
        h_clipped = np.minimum(h, 0.0) 
        Se = (1.0 + (self.alpha * np.abs(h_clipped))**self.n)**(-self.m)
        return self.theta_r + (self.theta_s - self.theta_r) * Se

    def calc_K(self, h: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Calculates hydraulic conductivity."""
        h_clipped = np.minimum(h, 0.0)
        Se = (1.0 + (self.alpha * np.abs(h_clipped))**self.n)**(-self.m)
        return self.Ks * (Se**0.5) * (1.0 - (1.0 - Se**(1.0 / self.m))**self.m)**2

    def calc_C_analytical(self, h: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Calculates analytical Specific Moisture Capacity (C = dθ/dh)."""
        C = np.full_like(h, 1e-8) 
        mask_neg = h < 0.0
        if np.any(mask_neg):
            h_abs = np.abs(h[mask_neg])
            num = self.alpha * self.n * self.m * (self.theta_s - self.theta_r) * (self.alpha * h_abs)**(self.n - 1.0)
            den = (1.0 + (self.alpha * h_abs)**self.n)**(self.m + 1.0)
            C[mask_neg] = np.maximum(num / den, 1e-12)
        return C

    def calc_C_chord(self, h_iter: npt.NDArray[np.float64], h_old: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Calculates C using chord slope to improve Picard stability in dry soils."""
        dh = h_iter - h_old
        C = np.full_like(h_iter, 1e-8)
        
        # Use chord slope where dh is large enough to avoid division by zero
        mask_diff = np.abs(dh) > 1e-4
        if np.any(mask_diff):
            C[mask_diff] = (self.calc_theta(h_iter[mask_diff]) - self.calc_theta(h_old[mask_diff])) / dh[mask_diff]
            
        # Fallback to analytical derivative where head hasn't changed much
        mask_small = ~mask_diff
        if np.any(mask_small):
            C[mask_small] = self.calc_C_analytical(h_iter[mask_small])
            
        return np.maximum(C, 1e-12)

# ==========================================
# 2. SOLVER CORE (MATH & MATRICES)
# ==========================================
def get_k_face(K: npt.NDArray[np.float64], h: npt.NDArray[np.float64], dz: float, mean_type: str) -> npt.NDArray[np.float64]:
    """Calculates inter-node hydraulic conductivity."""
    K_top, K_bottom = K[:-1], K[1:]
    
    if mean_type == 'upwind':
        h_top, h_bottom = h[:-1], h[1:]
        return np.where((h_top - h_bottom + dz) >= 0, K_top, K_bottom)
    elif mean_type == 'arithmetic':
        return 0.5 * (K_top + K_bottom)
    elif mean_type == 'geometric':
        return np.sqrt(K_top * K_bottom + 1e-20)
    else: # harmonic
        return (2 * K_top * K_bottom) / (K_top + K_bottom + 1e-20)

def get_fluxes(h: npt.NDArray[np.float64], soil: VanGenuchtenSoil, dz: float, N: int, 
               bc_top_flux: float, bc_bot_head: float, mean_type: str) -> npt.NDArray[np.float64]:
    """Calculates fluxes across all cell faces."""
    flux = np.zeros(N + 1)
    flux[0] = bc_top_flux 
    K = soil.calc_K(h)
    K_face = get_k_face(K, h, dz, mean_type)
    
    flux[1:N] = -K_face * ((h[1:] - h[:-1]) / dz - 1.0)
    flux[N] = -K[-1] * ((bc_bot_head - h[-1]) / (dz/2.0) - 1.0)
    return flux

def perform_timestep(h_old: npt.NDArray[np.float64], theta_old: npt.NDArray[np.float64], 
                     soil: VanGenuchtenSoil, dt: float, dz: float, N: int, 
                     bc_top: float, bc_bot: float, mean_type: str, 
                     max_iters: int = 20, tol_h: float = 1e-4, tol_res: float = 1e-5) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], int, bool]:
    """Executes one Picard Iteration timestep."""
    h_iter = np.copy(h_old)
    
    for iteration in range(max_iters):
        theta_iter = soil.calc_theta(h_iter)
        K = soil.calc_K(h_iter)
        C = soil.calc_C_chord(h_iter, h_old)
        K_face = get_k_face(K, h_iter, dz, mean_type)
        
        flux = get_fluxes(h_iter, soil, dz, N, bc_top, bc_bot, mean_type)
        residual = (theta_iter - theta_old)/dt + (flux[1:] - flux[:-1])/dz
        
        A_diag = C / dt
        A_upper = -K_face / dz**2
        A_lower = -K_face / dz**2
        
        A_diag[:-1] += K_face / dz**2
        A_diag[1:]  += K_face / dz**2
        A_diag[-1] += K[-1] / (dz/2.0 * dz)

        ab = np.zeros((3, N))
        ab[0, 1:] = A_upper
        ab[1, :]  = A_diag
        ab[2, :-1] = A_lower
        
        try:
            dh = solve_banded((1, 1), ab, -residual)
        except Exception:
            return h_old, theta_old, iteration, False

        h_iter += dh
        
        # Dual Convergence criteria (Head change AND mass residual)
        if np.max(np.abs(dh)) < tol_h and np.max(np.abs(residual)) < tol_res:
            return h_iter, soil.calc_theta(h_iter), iteration + 1, True
            
    return h_old, theta_old, max_iters, False
# ==========================================
# 3. SIMULATION ORCHESTRATION
# ==========================================
def simulate(cfg: Dict[str, Any], mode: str) -> Dict[str, Any]:
    print(f"Running {mode} over {cfg['solver_settings']['max_time']} hours")
    """Runs the simulation for a specific averaging mode."""
    soil = VanGenuchtenSoil(**cfg['soil_properties'])
    L, dz = cfg['grid_setup']['total_depth'], cfg['grid_setup']['node_spacing']
    N = int(L / dz)
    
    t, total_iters, step_count, total_mbe = 0.0, 0, 0, 0.0
    dt = cfg['solver_settings'].get('initial_dt', 1e-4)
    dt_min, dt_max = 1e-6, cfg['solver_settings'].get('max_dt', 1.0)
    
    #Define a maximum step limit (default to 50,000 if not in JSON)
    max_steps = cfg['solver_settings'].get('max_steps', 50000)
    
    h_state = np.full(N, cfg['conditions']['initial_head'])
    theta_state = soil.calc_theta(h_state)
    
    start_wall = time.time()
    simulation_finished = True

    while t < cfg['solver_settings']['max_time']:
        if step_count >= max_steps:
            #Catch stagnation if the solver is spinning its wheels
            simulation_finished = False
            break

        storage_before = np.sum(theta_state) * dz
        
        h_new, th_new, iters, success = perform_timestep(
            h_state, theta_state, soil, dt, dz, N, 
            cfg['conditions']['top_flux'], cfg['conditions']['bottom_head'], mode
        )
        
        if success:
            # Calculate Mass Balance Error
            fluxes = get_fluxes(h_new, soil, dz, N, cfg['conditions']['top_flux'], cfg['conditions']['bottom_head'], mode)
            inflow, outflow = cfg['conditions']['top_flux'] * dt, fluxes[-1] * dt
            storage_after = np.sum(th_new) * dz
            total_mbe += (storage_after - storage_before) - (inflow - outflow)
            
            # Advance State
            h_state, theta_state = np.copy(h_new), np.copy(th_new)
            t += dt
            total_iters += iters
            step_count += 1
            
            if step_count % 500 == 0:
                print(f"  -> Sim Time: {t:>6.3f} / {cfg['solver_settings']['max_time']} h | Step: {step_count:<7} | dt: {dt:.2e}", end='\r', flush=True)
            # live progress checker
            
            # Adaptive Time Stepping (Success)
            if iters <= 5:
                dt = min(dt * 1.25, dt_max)
            elif iters > 12:
                dt = max(dt * 0.8, dt_min)
            
            # Adaptive Time Stepping (Success)
            if iters <= 5:
                dt = min(dt * 1.25, dt_max)
            elif iters > 12:
                dt = max(dt * 0.8, dt_min)
        else:
            # Timestep Retry Mechanism (Failure)
            dt *= 0.5
            if dt < dt_min: 
                simulation_finished = False 
                break

    return {
        "mode": mode,
        "success": simulation_finished,
        "runtime": time.time() - start_wall,
        "avg_iters": total_iters / step_count if step_count > 0 else 0,
        "mbe": total_mbe,
        "final_theta": theta_state,
        "final_time": t
    }
# ==========================================
# 4. I/O & EXECUTION
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="1D Mixed-Form Richards Equation Solver")
    parser.add_argument("--config", type=str, default="params.json", help="Path to configuration JSON")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Config file '{args.config}' not found.")
        return

    with open(args.config, 'r') as f:
        cfg = json.load(f)

    target_setting = cfg['solver_settings'].get('averaging_mode', 'harmonic').lower()
    modes_to_run = ['arithmetic', 'geometric', 'harmonic', 'upwind'] if target_setting == "all" else [target_setting]

    L, dz = cfg['grid_setup']['total_depth'], cfg['grid_setup']['node_spacing']
    z_nodes = np.linspace(dz/2, L - dz/2, int(L/dz))
    
    plt.figure(figsize=(9, 7))
    initial_theta = VanGenuchtenSoil(**cfg['soil_properties']).calc_theta(np.full(int(L/dz), cfg['conditions']['initial_head']))
    plt.plot(initial_theta, -z_nodes, 'k--', label='Initial Condition', alpha=0.6)

    print("="*60)
    print(f"{'Mode':<12} | {'Status':<10} | {'Time(s)':<8} | {'Avg Iter':<8} | {'MBE (cm)':<10}")
    print("-" * 60)

    # List to hold results for the text file
    report_lines = []

    for mode in modes_to_run:
        result = simulate(cfg, mode)
        status_str = "SUCCESS" if result['success'] else f"FAIL @ {result['final_time']:.2f}"
        plot_label = f"{mode} ({status_str})"
        
        # Print to terminal
        print(f"{result['mode']:<12} | {status_str:<10} | {result['runtime']:<8.2f} | {result['avg_iters']:<8.2f} | {result['mbe']:<10.2e}")
        
        # Save for text file
        report_lines.append(f"{result['mode']:<12} | {status_str:<18} | {result['runtime']:<10.2f} | {result['avg_iters']:<10.2f} | {result['mbe']:<12.2e}\n")
        
        plt.plot(result['final_theta'], -z_nodes, label=plot_label, linewidth=2)

    print("="*60)

    # WRITE THE REPORT FILE
    with open("simulation_report.txt", "w") as f:
        f.write("MIXED-FORM RICHARDS SOLVER REPORT\n")
        f.write("="*80 + "\n")
        f.write(f"{'Mode':<12} | {'Status':<18} | {'Time (s)':<10} | {'Avg Iter':<10} | {'MBE (cm)':<12}\n")
        f.write("-" * 80 + "\n")
        f.writelines(report_lines)
        f.write("="*80 + "\n")
        f.write(f"Timestamp: {time.ctime()}\n")

    plt.xlabel('Volumetric Water Content (θ)')
    plt.ylabel('Depth (cm)')
    plt.title(f'1D Mixed-Form Moisture Profiles (t = {cfg["solver_settings"]["max_time"]}h)')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.show()

if __name__ == "__main__":
    main()