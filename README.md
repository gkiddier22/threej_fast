# threej_fast

Efficient OpenMP-parallelized computation of Wigner 3j symbol-based mode-coupling matrices for CMB analysis.

This repository provides a lightweight, high-performance C implementation for generating mode-coupling matrices $K^{EE}$ and $K^{TT}$, using Wigner 3j symbols. It is designed for speed and portability, and includes example usage through a Python script.

## 🔧 Contents
```bash
threej_fast/
├── src/
│   ├── thrj_000_opt.c         # Computes scalar 3j symbols
│   └── thrj_220_opt.c         # Computes spin-2 3j symbols (EE-like)
│
├── tests/
│   └── coupling_mat_example.py  # Example script using thrj_XXX_opt
│
├── data/                      # Folder for any example input data files
│
├── Makefile                   # Builds both C programs with OpenMP
├── 3j_env.yml                 # Conda environment file
├── citation.cff              # Citation metadata
└── README.md                  # Project overview and usage
```

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
conda activate 3j_env
```

### 4. Run example code 
Change the pathnames for the variables `MASK_FILE_1` and `MASK_FILE_2` in the file `coupling_mat_example_EE.py` to .fits files containing the masks. Then run the script:

```bash
python tests/coupling_mat_example_EE.py
```
This script will:
- Load input window function data from the data/ directory
- Call the compiled thrj_220_opt or thrj_000_opt binary
- Plot the resulting coupling matrix for visual inspection

## 📁 Input Format

The C binaries expect input files containing one column of window function values for each \ell, from \ell = 0 to \ell_{\max}.

