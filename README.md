# threej_fast

Efficient OpenMP-parallelized computation of Wigner 3j symbol-based mode-coupling matrices for CMB analysis.

This repository provides a lightweight, high-performance C implementation for generating mode-coupling matrices \( K^{EE} \) and \( K^{TT} \), using Wigner 3j symbols. It is designed for speed and portability, and includes example usage through a Python script.

## 🔧 Contents

- `thrj_220_opt.c`: Optimized C code using OpenMP (recommended for Linux)
- `thrj_000_opt.c`: Compatible version without OpenMP (for macOS/Clang)
- `tests/coupling_mat_example.py`: Example script for generating and plotting coupling matrices
- `data/`: Directory for input window functions and output matrices
- `3j_env.yml`: Conda environment file for Python dependencies
- `citation.cff`: Citation metadata for referencing this work

## 🚀 Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/gkiddier22/threej_fast.git
cd threej_fast
```

### 2. Compilation
Use the provided Makefile to compile both optimized executables:
```bash
make
```
This will build:
- thrj_220_opt: EE/EB coupling matrix calculator 
- thrj_000_opt: TT coupling matrix calculator 

 ### 3. Python environment (optional for plotting and testing)
 ```bash
conda env create -f 3j_env.yml
conda activate threej_env
```

### 4. Run example code 
```bash
python tests/coupling_mat_example.py
```
This script will:
- Load input window function data from the data/ directory
- Call the compiled thrj_220_opt or thrj_000_opt binary
- Plot the resulting coupling matrix for visual inspection

## 📁 Input Format

The C binaries expect input files containing one column of window function values for each \ell, from \ell = 0 to \ell_{\max}.

