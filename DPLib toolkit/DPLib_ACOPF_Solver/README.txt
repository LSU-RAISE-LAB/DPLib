======================================================================
Distributed AC Optimal Power Flow (ACOPF) Solver Using ADMM
======================================================================

Author      : Milad Hasanzadeh  
Affiliation : Department of Electrical and Computer Engineering  
              Louisiana State University  
Email       : e.mhasanzadeh1377@yahoo.com  
Release     : June 2025  
License     : Academic and Research Use Only  
======================================================================

Overview
--------
This MATLAB code implements a **distributed AC Optimal Power Flow (ACOPF)** solver based on the **Alternating Direction Method of Multipliers (ADMM)**. The solver is designed for use on large-scale power networks that have been partitioned into multiple regions using a graph-based method.

Each region solves its own ACOPF problem locally while coordinating with neighboring regions through consistency constraints at shared boundary buses (tie-lines). This framework supports privacy, modularity, and scalability in power system operation.

The solver works directly with the output of the **Partitioning Code** and uses **IPOPT** solver.

======================================================================

Requirements
------------
- MATLAB R2020a or newer  
- [MATPOWER toolbox](https://matpower.org) (must be added to MATLAB path)  
- [IPOPT solver](https://coin-or.github.io/Ipopt/) installed and configured with YALMIP  
- The `.mat` file created by the Partitioning Code

Ensure all toolboxes and solvers are properly installed and accessible before running the code.

======================================================================
Files Included
--------------

- DistACOPF.m            : Main driver script that initializes and runs the distributed ADMM-based ACOPF solver. It loads the partitioned MATPOWER case, sets ADMM parameters, calls regional solvers, updates consistency variables and dual multipliers, and checks convergence.

- acopf.m                : Core helper function that solves the local ACOPF problem for a given region using IPOPT. It receives the region data and optimization structures, and returns the optimal state variables and solver status.

- bounds.m               : Defines lower and upper bounds for the optimization variables in each region, including generator limits, voltage bounds, and power limits.

- constraintbounds.m     : Constructs the lower and upper bounds for the equality and inequality constraints, including power balance and generator output restrictions.

- constraints.m          : Evaluates all nonlinear constraints for a given region based on its decision variables. This includes nodal power balances and operating limits.

- doanalysis.m           : Performs post-processing of the ACOPF solution, including reporting mismatch values, verifying convergence, and computing residuals.

- gradient.m             : Computes the gradient of the regional objective function with respect to the optimization variables, needed by IPOPT.

- hessian.m              : Computes the Hessian of the Lagrangian (second-order derivatives), combining both objective and constraint contributions for each region.

- hessianstructure.m     : Returns the sparsity pattern (structure) of the Hessian matrix for use by IPOPT, allowing more efficient interior-point optimization.

- initialx0.m            : Generates an initial guess for the decision variables in each region, including voltage magnitudes, angles, and power outputs.

- jacobian.m             : Computes the Jacobian matrix of constraint functions (i.e., first derivatives of all constraints) required by IPOPT.

- jacobianstructure.m    : Defines the sparsity pattern of the Jacobian matrix, enabling IPOPT to exploit structural sparsity for faster computation.

- objective.m            : Evaluates the regional cost function (e.g., generator quadratic cost) for the current values of decision variables.

- preprocess.m           : Prepares and extracts necessary data for each region from the global MATPOWER case, including region-specific buses, branches, and tie-lines. It also defines the indexing for shared variables and constraints.

- whichstatus.m          : Interprets the IPOPT solver exit flags and provides user-friendly messages to indicate success, failure, or other solver outcomes.



======================================================================

How to Use
----------
1. Open MATLAB and navigate to the folder containing `DistACOPF.m`.  
2. Ensure your partitioned system `.mat` file is also in the same folder.  
3. Run the script:

   >> DistACOPF

4. When prompted, provide the following:
   - Name of the partitioned `.mat` file (e.g., `pglib_opf_case300_ieee_4regions.mat`)
   - ADMM penalty parameter (ρ)
   - Stopping threshold for the worst primal residual
   - Maximum number of iterations
   - Centralized OPF cost (optional; enter 0 to skip optimality gap calculation)

The solver will:
- Load the partitioned system
- Augment each region with interregional tie-line data
- Initialize required variables and parameters
- Enter an ADMM loop that solves local ACOPF problems, updates multipliers, and checks convergence

======================================================================

Code Structure
--------------
1. **Loading the Partitioned Data**  
   The script loads bus, branch, and generator data for all regions. Tie-line information is also processed to connect the regions for distributed coordination.

2. **Tie-Line Augmentation**  
   Each region's branch data is expanded to include the tie-lines it shares with neighbors. Additional buses are added with adjusted indices to maintain topology integrity.

3. **Initialization**  
   Important variables like the number of buses, tie-lines, and interface counters are initialized. Dual variables and phase angle/magnitude tracking arrays are also created.

4. **User Input for ADMM Parameters**  
   The user defines the penalty parameter (ρ), residual threshold, and max iterations. Optionally, the centralized cost is used to compute the optimality gap.

5. **Main ADMM Iteration Loop**  
   At each iteration, all regions solve their local ACOPF problems using the `ACOPF_SP.m` helper function. The Lagrange multipliers (dual variables) are then updated based on discrepancies in boundary bus values.

6. **Convergence Tracking**  
   The code calculates the maximum residual error and logs it across iterations. If the maximum error falls below the threshold, convergence is declared.

7. **Visualization**  
   After completion, the script generates:
   - A log-scale plot of primal residual errors over iterations
   - (If applicable) a plot of optimality gap percentages

======================================================================

Helper Function (ACOPF.m)
----------------------------
This file defines and solves the regional ACOPF subproblem using YALMIP and IPOPT. Here’s what it does:

- Retrieves region-specific data and interface mappings
- Defines decision variables:  
  * Real power \( P_g \)  
  * Reactive power \( Q_g \)  
  * Voltage magnitude \( V \)  
  * Voltage angle \( \theta \)

- Constructs the local AC power flow constraints:
  * Generator limits
  * Voltage magnitude limits
  * Real and reactive power balances at each bus
  * Line flow thermal limits

======================================================================

Outputs
-------
1. **Plots**
   - `Error Progression`: maximum primal residual per iteration
   - `Optimality Gap`: shows the percentage gap to centralized cost if provided

2. **Command Window Output**
   - Iteration number
   - Current maximum residual error
   - Whether convergence was achieved

======================================================================

Notes
-----
- The distributed ACOPF formulation follows a standard ADMM structure where each region is decoupled except at tie-line interfaces.
- Users can change the stopping criterion or penalty parameter to control the speed and accuracy of convergence.
- The code is modular and can be adapted for future improvements or extended to handle uncertainty, time-varying loads, or contingencies.

======================================================================

Citation
--------
If you use this code in your work, please cite the following:

M. Hasanzadeh and A. Kargarian, "DPLib: A Standard Benchmark Library for Distributed Power System Analysis and Optimization," http://arxiv.org/abs/2506.20819

======================================================================

Contact
-------
For any questions, bug reports, or suggestions, contact:  
Milad Hasanzadeh  
e.mhasanzadeh1377@yahoo.com
