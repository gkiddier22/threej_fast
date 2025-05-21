"""
Example computation of coupling matrix for EE Power Spectrum

This script:
1. Reads HEALPix masks and computes the mask window function.
2. Calls a compiled C program to compute the coupling matrix.
3. Loads and visualizes the resulting coupling matrix.
4. Saves the plot 

Author: Georgia Kiddier  
Date: 19-02-2025  
License: [Specify License]  
"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
import sys
from matplotlib.colors import LogNorm

# ---- Configuration ----
L_MAX = 1996  # Maximum multipole moment
MASK_FILE_1 = "./COM_Mask_Likelihood-polarization-143-hm1_2048_R3.00.fits"
MASK_FILE_2 = "./COM_Mask_Likelihood-polarization-143-hm2_2048_R3.00.fits"
C_PROGRAM = "./src/thrj_220_opt"  # Path to the compiled C program
OUTPUT_MATRIX_FILE = "./results/coupling_matrix_EE_opt.dat"  # Output binary file for the matrix
PLOT_SAVE_PATH = "./results/coupling_matrix_EE_opt.png"  # Output path for the plot
SPECTRUM_TYPE = "EE" 

# ---- Step 1: Compute Mask Window Function ----
def compute_mask_window_function(lmax, mask1_path, mask2_path):
    """
    Computes the window function from two HEALPix mask files.

    Parameters:
    -----------
    lmax : int
        Maximum multipole moment.
    mask1_path : str
        Path to the first mask file.
    mask2_path : str
        Path to the second mask file.

    Returns:
    --------
    np.ndarray
        The computed window function.
    str
        Filename where the window function is saved.
    """
    print("Reading HEALPix masks...")
    mask_1 = hp.read_map(mask1_path)
    mask_2 = hp.read_map(mask2_path)

    print("Computing the window function...")
    mask_cl = hp.anafast(mask_1, map2=mask_2, lmax=2*lmax, use_pixel_weights=True)
    window_filename = "./results" + f"w_l_{2*lmax}_{SPECTRUM_TYPE}.dat"
    
    mask_cl.tofile(window_filename)
    print(f"Window function saved to {window_filename}")

    return mask_cl, window_filename

# ---- Step 2: Run C Program for Coupling Matrix ----
def run_c_program(c_program, wlfname, output_matrix_file, lmax, spectrum_type):
    """
    Calls the external C program to compute the coupling matrix.

    Parameters:
    -----------
    c_program : str
        Path to the compiled C executable.
    wlfname : str
        Path to the computed window function file.
    output_matrix_file : str
        Path to save the computed matrix.
    lmax : int
        Maximum multipole moment.
    spectrum_type : str
        Type of power spectrum (e.g., "EE").
    """
    args = [c_program, wlfname, output_matrix_file, str(lmax), spectrum_type, "yes"]
    print(f"Running C program: {' '.join(args)}")

    try:
        subprocess.run(args, check=True)
        print(f"Coupling matrix saved to {output_matrix_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running the C program: {e}")
        sys.exit(1)

# ---- Step 3: Load and Plot the Coupling Matrix ----
def plot_coupling_matrix(matrix_file, lmax, save_path):
    """
    Loads the computed coupling matrix and visualizes it.

    Parameters:
    -----------
    matrix_file : str
        Path to the binary file containing the computed matrix.
    lmax : int
        Maximum multipole moment for reshaping the matrix.
    save_path : str
        Path to save the output plot.
    """
    print("Loading coupling matrix...")
    matrix = np.fromfile(matrix_file, dtype="float64").reshape((lmax + 1, lmax + 1))

    print("Plotting the coupling matrix...")
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    norm = LogNorm()  
    im = ax.imshow(matrix, cmap="viridis", norm=norm)

    ax.set_title(f"Coupling Matrix ({SPECTRUM_TYPE})")
    ax.set_xlabel(r"$\ell_1$")
    ax.set_ylabel(r"$\ell_2'$")

    cbar = fig.colorbar(im, ax=ax, orientation="vertical", fraction=0.05, pad=0.02)
    cbar.set_label("Coupling Strength")

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved to: {save_path}")

    plt.show()

# ---- Main Execution ----
if __name__ == "__main__":
    # Compute the window function
    wl, wl_filename = compute_mask_window_function(L_MAX, MASK_FILE_1, MASK_FILE_2)

    # Run the C program to compute the coupling matrix
    run_c_program(C_PROGRAM, wl_filename, OUTPUT_MATRIX_FILE, L_MAX, SPECTRUM_TYPE)

    # Plot and save the coupling matrix
    plot_coupling_matrix(OUTPUT_MATRIX_FILE, L_MAX, PLOT_SAVE_PATH)
