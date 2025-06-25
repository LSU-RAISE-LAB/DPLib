======================================================================
 Partitioning Code for Distributed Power System Applications
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
This MATLAB-based code partitions a centralized MATPOWER case (e.g., case118, case1354, etc.) into multiple electrically coherent regions using a graph-theoretic approach. The resulting partitioned system is saved in formats that are directly usable in distributed DC and AC OPF simulations. It supports reproducible distributed studies and is compatible with the distributed solvers released alongside this work.

======================================================================

Requirements
------------
- MATLAB R2020a or newer  
- The following MATLAB toolboxes must be installed and added to the MATLAB path:
  * MATPOWER (https://matpower.org)
  * Statistics and Machine Learning Toolbox (for k-means clustering)

======================================================================

Files Included
--------------
- `partitioning_code.m`       : Main script to run the partitioning process.

======================================================================

How to Use
----------
1. Open MATLAB and set the current folder to where `partitioning_code.m` is located.
2. Ensure MATPOWER is installed and its path is added to MATLAB using `addpath(genpath('matpower_directory'))`.
3. Run the script by typing in the MATLAB command window:
   
   >> partitioning_code

4. The script will prompt you to:
   a) Enter the MATPOWER case name (e.g., case118)
   b) Enter the number of regions (e.g., 3)

5. The code will:
   - Load the original case data
   - Build the network graph using bus and branch connectivity
   - Compute the Laplacian matrix
   - Perform spectral clustering followed by k-means
   - Assign buses and branches to regions
   - Detect and store inter-regional tie-lines
   - Save results in `.mat`, `.m`, and `.png` formats

======================================================================

Outputs
-------
After successful execution, the following are generated:

1. **MAT File (.mat)**  
   Contains:
   - Each region's bus, branch, and generator data in MATPOWER format
   - Tie-line table indicating connections between regions
   - Number of regions and slack bus region ID

2. **M File (.m)**  
   - Saves regional assignments and parameters as MATLAB structures
   - Can be reused for visualization or debugging

3. **PNG Figure**  
   - Shows the entire system graph with colored nodes indicating regions
   - Tie-lines are highlighted between regions

4. **Command Window Output**  
   - Prints detailed information such as:
     * If bus renumbering was needed
     * Graph connectivity verification
     * Total number of regions and tie-lines
     * Slack bus location
     * Confirmation that all buses, branches, and generators are included
     * Output filenames and locations

======================================================================

Citation
--------
If you use this code in your academic work, please cite the following:

M. Hasanzadeh, "Toward a Standard Repository for Distributed Power System Optimization and Control," 2025.

======================================================================

Contact
-------
For questions, suggestions, or bug reports, please contact:  
Milad Hasanzadeh  
e.mhasanzadeh1377@yahoo.com

