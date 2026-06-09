# DPLib: A Standard Benchmark Library for Distributed Power System Analysis and Optimization

**DPLib** is an open-source MATLAB-based benchmark library for **distributed power system analysis, optimization, planning, and control**. It provides a standard and reproducible framework for converting centralized MATPOWER/PGLib test systems into **multi-region distributed benchmark cases**.

The main purpose of DPLib is to provide the **data foundation** needed for distributed power-system research. Instead of requiring researchers to manually divide centralized test systems into areas, DPLib provides a graph-based partitioning tool and ready-to-use distributed datasets with regional network data, tie-line information, and boundary connections.

DPLib can be used for a wide range of distributed power-system studies, including:

* distributed DC and AC optimal power flow,
* distributed generation expansion planning,
* distributed transmission expansion planning,
* distributed voltage-stability analysis,
* distributed state estimation,
* contingency analysis,
* decomposition-based optimization,
* privacy-preserving grid operation,
* distributed and decentralized control.

The included distributed DC OPF and AC OPF solvers are provided mainly as **validation tools** and example workflows to demonstrate that the generated distributed datasets can be solved and benchmarked using distributed optimization methods.

---

## Key Features

DPLib includes the following main components.

### 1. Standard Distributed Benchmark Datasets

DPLib provides ready-to-use multi-region benchmark datasets generated from standard MATPOWER/PGLib cases. Each distributed case represents the original centralized power network as a set of connected regions.

Each region contains its own:

* buses,
* generators,
* branches,
* local network data,
* boundary buses,
* tie-line information,
* neighboring-region connections.

These datasets allow researchers to test distributed methods on standard and reproducible power-system cases without manually creating multi-area models.

---

### 2. Graph-Based Partitioning Tool

This tool converts any centralized MATPOWER/PGLib case into a distributed multi-region benchmark case.

The partitioning tool uses:

* graph representation of the power network,
* graph Laplacian construction,
* spectral clustering,
* k-means clustering,
* optional weighted network partitioning,
* optional user-defined bus-to-region mapping.

The tool generates:

* distributed `.mat` files,
* distributed `.m` files,
* regional MATPOWER structures,
* tie-line and boundary information,
* region-to-bus assignments,
* topology visualization figures.

The partitioning tool can be used in two main modes:

1. **Automatic partitioning**, where the tool partitions the system using spectral clustering.
2. **User-defined partitioning**, where the user provides a custom `bus_region_map`.

This makes DPLib flexible for both benchmark generation and customized distributed power-system studies.

---

### 3. Automatic Benchmark-Generation Runner

DPLib also includes an automatic MATLAB runner for generating and testing distributed benchmark cases:

```matlab
run_DPLib_partition_OPF
```

This function can automatically:

1. take centralized MATPOWER/PGLib case names,
2. partition each case into a specified number of regions,
3. load and verify the generated distributed case,
4. summarize the number of buses, branches, generators, regions, and tie-lines,
5. optionally run distributed OPF validation studies,
6. save convergence figures,
7. save summary tables in `.mat`, `.csv`, and `.xlsx` formats,
8. print compact LaTeX tables for papers, appendices, and response letters.

The runner supports different modes:

```matlab
'partition'   % only generate distributed benchmark cases
'dc'          % generate case and run distributed DC OPF validation
'ac'          % generate case and run distributed AC OPF validation
'both'        % generate case and run both DC and AC OPF validation
'none'        % same as partition-only mode
```

Example:

```matlab
case_list = {
    'pglib_opf_case4917_goc', 10
};

[DPLibMainSummary, all_results] = run_DPLib_partition_OPF(case_list, 'partition');
```

Example with OPF validation:

```matlab
[DPLibMainSummary, all_results] = run_DPLib_partition_OPF(case_list, 'both', ...
    'force_partition', true, ...
    'use_weighted', false, ...
    'bus_region_map', [], ...
    'rho_dc', 5000, ...
    'rho_ac', 5000, ...
    'thresholdDC', 1e-4, ...
    'thresholdAC', 1e-4, ...
    'maxIterDC', 3000, ...
    'maxIterAC', 3000, ...
    'holds', 1e8, ...
    'saveFigures', true, ...
    'saveSummary', true, ...
    'summaryName', 'DPLib_OPF_summary');
```

---

## Repository Structure

```text
DPLib/
│
├── DPLib toolkit/
│   │
│   ├── DPLib_partitioning_tool/
│   │   └── Graph-based spectral partitioning tool
│   │
│   ├── DPLib_main_runner/
│   │   └── Automatic distributed benchmark-generation runner
│   │
│   ├── DPLib_DCOPF_Solver/
│   │   └── Distributed DC OPF validation solver based on ADMM and YALMIP
│   │
│   └── DPLib_ACOPF_Solver/
│       └── Distributed AC OPF validation solver based on ADMM and IPOPT
│
├── pglib_opf_case5_pjm/
├── pglib_opf_case24_ieee_rts/
├── pglib_opf_case30_ieee/
├── pglib_opf_case39_epri/
├── pglib_opf_case57_ieee/
├── pglib_opf_case60_c/
├── pglib_opf_case89_pegase/
├── pglib_opf_case118_ieee/
├── pglib_opf_case162_ieee_dtc/
├── pglib_opf_case179_goc/
├── pglib_opf_case197_snem/
├── pglib_opf_case200_tamu/
├── pglib_opf_case240_pserc/
├── pglib_opf_case300_ieee/
├── pglib_opf_case793_goc/
├── pglib_opf_case1354_pegase/
├── pglib_opf_case1803_snem/
├── pglib_opf_case1888_rte/
├── pglib_opf_case1951_rte/
├── pglib_opf_case2000_tamu/
├── pglib_opf_case2316_sdet/
├── pglib_opf_case2383wp_k/
├── pglib_opf_case2736sp_k/
├── pglib_opf_case2742_goc/
├── pglib_opf_case2746wp_k/
├── pglib_opf_case2869_pegase/
├── pglib_opf_case4020_goc/
├── pglib_opf_case4917_goc/
├── pglib_opf_case7336_epigrids/
├── pglib_opf_case10192_epigrids/
├── pglib_opf_case10480_goc/
├── pglib_opf_case13659_pegase/
├── pglib_opf_case20758_epigrids/
└── ...
```

Each `pglib_opf_case*` folder contains one or more distributed benchmark cases generated from a centralized MATPOWER/PGLib system. 

---

## What Is a DPLib Case?

A standard centralized MATPOWER/PGLib case represents the full power network as one integrated system.

A DPLib case represents the same network as a set of connected regions. Each region contains local power-system data, while inter-regional connections are represented through tie-lines and boundary information.

For example, a centralized case such as:

```text
pglib_opf_case4917_goc
```

can be partitioned into 10 regions and saved as:

```text
pglib_opf_case4917_goc_10regions.mat
pglib_opf_case4917_goc_10regions.m
pglib_opf_case4917_goc_10regions_topology.png
```

This distributed case can then be used as input data for distributed optimization, planning, control, and analysis algorithms.

---

## Requirements

DPLib requires the following software packages:

* MATLAB R2020a or newer
* MATPOWER
* PGLib-OPF cases
* Statistics and Machine Learning Toolbox

Additional packages are required only when using the OPF validation solvers:

* YALMIP, required for the distributed DC OPF validation solver
* IPOPT, required for the distributed AC OPF validation solver

Before using DPLib, make sure all required folders are added to the MATLAB path. For AC OPF validation, IPOPT must also be installed and properly linked to MATLAB.

---

## Basic Usage

### 1. Generate a Distributed Benchmark Case

```matlab
case_name = 'pglib_opf_case4917_goc';
num_regions = 10;
use_weighted = false;
bus_region_map = [];

result = partitioning_code(case_name, num_regions, use_weighted, bus_region_map);
```

This generates:

```text
pglib_opf_case4917_goc_10regions.mat
pglib_opf_case4917_goc_10regions.m
pglib_opf_case4917_goc_10regions_topology.png
```

---

### 2. Generate Multiple Distributed Benchmark Cases

```matlab
case_list = {
    'pglib_opf_case1803_snem',    6
    'pglib_opf_case2736sp_k',     6
    'pglib_opf_case2742_goc',     7
    'pglib_opf_case57_ieee',      2
    'pglib_opf_case60_c',         2
    'pglib_opf_case4020_goc',     8
    'pglib_opf_case5_pjm',        2
    'pglib_opf_case39_epri',      2
};

[DPLibMainSummary, all_results] = run_DPLib_partition_OPF(case_list, 'partition');
```

This generates distributed benchmark datasets and summary files without running OPF validation.

---

### 3. Run the Full Workflow with OPF Validation

```matlab
case_list = {
    'pglib_opf_case4917_goc', 10
};

[DPLibMainSummary, all_results] = run_DPLib_partition_OPF(case_list, 'both');
```

This generates the distributed case, runs centralized and distributed OPF validation, saves convergence plots, and creates summary files.

---

## Output Files

Depending on the selected workflow, DPLib can generate the following outputs.

### Partitioning Outputs

```text
caseName_Nregions.mat
caseName_Nregions.m
caseName_Nregions_topology.png
```

### Summary Outputs

```text
DPLib_main_runner_summary.mat
DPLib_main_runner_summary.csv
DPLib_main_runner_summary.xlsx
```

The summary table includes:

* centralized case name,
* partitioned case name,
* number of regions,
* number of buses,
* number of branches,
* number of generators,
* number of tie-lines,
* partitioning time,
* solver validation status, if OPF validation is selected,
* error messages, if any.

### Optional OPF Validation Figure Outputs

If the OPF validation solvers are used, DPLib can also generate:

```text
caseName_Nregions_DC_error.png
caseName_Nregions_DC_gap.png
caseName_Nregions_DC_rho.png

caseName_Nregions_AC_error.png
caseName_Nregions_AC_gap.png
caseName_Nregions_AC_rho.png
```

---

## Python Interface

DPLib cases can also be loaded in Python using the provided helper function. This allows users to access regional MATPOWER data structures, including bus, generator, branch, and cost matrices, for use in Python-based optimization, planning, voltage-stability analysis, or machine-learning workflows.

---

## Applications

DPLib is designed as a distributed benchmark-data library. The generated multi-region datasets can be used as base data for many distributed power-system research problems.

Possible applications include:

* distributed optimal power flow,
* distributed generation expansion planning,
* distributed transmission expansion planning,
* multi-area unit commitment,
* distributed voltage-stability analysis,
* distributed state estimation,
* distributed contingency analysis,
* distributed reliability assessment,
* distributed restoration and resilience studies,
* distributed market-clearing studies,
* decomposition-based optimization,
* privacy-preserving grid operation,
* distributed and decentralized control.

DPLib provides the network and regional decomposition structure required for these studies. Depending on the application, users may add extra problem-specific data on top of the generated DPLib cases.

---

## Notes on Planning, Unit Commitment, and Voltage-Stability Studies

DPLib provides:

* network topology,
* regional decomposition,
* bus data,
* generator data,
* branch data,
* tie-line data,
* boundary-bus information,
* region-level MATPOWER structures.

These data are useful base inputs for many distributed power-system problems.

For **distributed generation expansion planning**, users may add:

* candidate generator locations,
* investment costs,
* technology types,
* planning horizon,
* demand growth scenarios,
* renewable generation scenarios,
* budget constraints.

For **distributed transmission expansion planning**, users may add:

* candidate transmission lines,
* candidate corridor data,
* line construction costs,
* planning horizon,
* contingency scenarios,
* reliability criteria,
* budget constraints.

For **multi-area unit commitment**, users may add:

* time horizon,
* hourly load profiles,
* generator minimum up/down times,
* startup and shutdown costs,
* ramping limits,
* initial generator status,
* reserve requirements,
* renewable generation profiles.

For **distributed voltage-stability analysis**, DPLib cases are more directly usable because they already contain AC network data, voltage magnitudes, reactive-power data, generator limits, branch parameters, and regional tie-line information.

---

## OPF Solvers as Validation Tools

DPLib includes distributed DC OPF and AC OPF solvers as example validation workflows.

These solvers are not the only purpose of DPLib. They are included to show that the generated distributed benchmark cases can be solved, tested, and compared against centralized OPF results.

### Distributed DC OPF Validation Solver

The distributed DC OPF solver is based on ADMM and implemented using MATLAB and YALMIP. Each region solves a local DC OPF problem and coordinates with neighboring regions through boundary variables and tie-line information.

### Distributed AC OPF Validation Solver

The distributed AC OPF solver is based on ADMM and nonlinear programming. It uses IPOPT to solve regional AC OPF subproblems and validates the generated distributed cases under nonlinear AC network constraints.

These validation solvers help users check:

* correctness of generated distributed datasets,
* consistency between centralized and distributed formulations,
* convergence behavior,
* scalability with respect to network size and number of regions.

---

## Citation

If you use DPLib, its generated datasets, partitioning tool, or validation solvers in your research, please cite:

> M. Hasanzadeh and A. Kargarian, “DPLib: A Standard Benchmark Library for Distributed Power System Analysis and Optimization,” [arXiv:2506.20819, 2025](https://arxiv.org/abs/2506.20819).

BibTeX:

```bibtex
@article{hasanzadeh2025dplib,
  title={DPLib: A Standard Benchmark Library for Distributed Power System Analysis and Optimization},
  author={Hasanzadeh, Milad and Kargarian, Amin},
  journal={arXiv preprint arXiv:2506.20819},
  year={2025}
}
```

---

## Author

**Milad Hasanzadeh**
Department of Electrical and Computer Engineering
Louisiana State University
Email: [e.mhasanzadeh1377@yahoo.com](mailto:e.mhasanzadeh1377@yahoo.com)

---

## Release Information

**Initial release:** June 2025
**Current version:** Updated DPLib benchmark and automatic runner version
**License:** Academic and Research Use Only

