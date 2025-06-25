======================================================================
Distributed DC Optimal Power Flow (DCOPF) Solver Using ADMM
======================================================================

Author      : Milad Hasanzadeh  
Affiliation : Department of Electrical and Computer Engineering,  
              Louisiana State University  
Email       : e.mhasanzadeh1377@yahoo.com  
Release     : June 2025  
License     : Academic/Research Use Only  
======================================================================

Overview
--------
This MATLAB script implements a **distributed DC Optimal Power Flow (DCOPF)** solver based on the **Alternating Direction Method of Multipliers (ADMM)**. It is designed to solve large-scale OPF problems by decomposing the system into multiple regions, each solving its subproblem locally, while coordinating across tie-lines to achieve consensus.

The solver is compatible with the output of the **Partitioning Code**, which prepares the multi-region MATPOWER data format.

======================================================================

Requirements
------------
- MATLAB R2020a or newer  
- [YALMIP toolbox](https://yalmip.github.io/download/) (added to MATLAB path)  
- [MATPOWER toolbox](https://matpower.org) (added to MATLAB path)

Ensure that both toolboxes are correctly installed and added to your MATLAB path before running the code.

======================================================================

Files Included
--------------
- `DistDCOPF.m`            : Main script to run the distributed DCOPF using ADMM  
- `DCOPF_SP.m`             : Helper function that solves the regional DCOPF subproblem  

======================================================================

How to Use
----------
1. Launch MATLAB and make sure the current directory includes `DistDCOPF.m` and the `partitioned_case.mat` file.
2. Run the script:

   >> DistDCOPF

3. When prompted, provide the following inputs:
   - Name of the partitioned `.mat` file (e.g., `pglib_opf_case300_ieee_4regions.mat`)
   - Penalty parameter (rho) for ADMM
   - Convergence threshold (worst primal residual)
   - Maximum number of iterations
   - Centralized optimal cost (if known, for calculating optimality gap; otherwise enter 0)

4. The script will:
   - Load the partitioned data
   - Augment each region with tie-line data
   - Initialize interface and dual variables
   - Enter the ADMM loop to solve each regional subproblem
   - Update Lagrange multipliers (dual variables)
   - Track convergence metrics and plot the results

======================================================================

Code Structure
--------------
The code proceeds through the following steps:

1. **Load Partitioned Data**  
   - Verifies and loads the `.mat` file that includes region-wise bus, branch, and generator data.

2. **Tie-Line Augmentation**  
   - Appends tie-line data to each region’s branch matrix.

3. **Variable Initialization**  
   - Initializes variables such as region sizes, dual multipliers, phase angles, and counters.

4. **ADMM Parameter Input**  
   - Reads user-defined penalty parameter, stopping threshold, and maximum iteration count.

5. **Subproblem Solver (DCOPF_SP)**  
   - Solves each regional DCOPF using YALMIP.  
   - Includes generator constraints, power balance, and line flow limits.  
   - Adds consistency terms to enforce agreement across tie-lines.

6. **Dual Variable Update**  
   - Updates Lagrange multipliers based on tie-line phase angle differences.

7. **Convergence Check**  
   - Tracks worst primal residual and computes the optimality gap (if applicable).  
   - Stops if the residual is below threshold or if the max iterations are reached.

8. **Plot Results**  
   - Error progression (log scale)  
   - Optional: Optimality gap over iterations

======================================================================

Outputs
-------
1. **Plots**
   - `Error Progression`: shows the maximum primal residual over iterations (log scale)
   - `Optimality Gap`: (if centralized cost is given) shows percentage gap between distributed and centralized solutions

2. **Printed Information in MATLAB Command Window**
   - Iteration-wise primal residuals
   - Convergence status and stopping condition
   - (Optional) Optimality gap at each iteration

======================================================================

Notes
-----
- The helper script `DCOPF_SP.m` automatically fetches each region’s parameters from the base workspace using `evalin` and `assignin`.  
- Each region’s bus angles are stored in symbolic form and updated at every iteration.  
- The power balance and line constraints are modeled using a standard DC power flow approximation.

======================================================================

Citation
--------
If you use this solver in your research or publication, please cite the following:

M. Hasanzadeh, "Toward a Standard Repository for Distributed Power System Optimization and Control," 2025.

======================================================================

Contact
-------
For questions, suggestions, or bug reports, please contact:  
Milad Hasanzadeh  
e.mhasanzadeh1377@yahoo.com
