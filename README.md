# DPLib: A Standard Benchmark Library for Distributed Power System Analysis and Optimization

**DPLib** is a comprehensive MATLAB-based framework for **Distributed Power System Analysis, Optimization, and Control**. It establishes the **first standard distributed dataset repository** for power systems by partitioning centralized MATPOWER cases into electrically coherent regions and benchmarking them using modular solvers.

---

## 🔧 What’s Included

This repository includes the following three core tools:

1. **Graph-Based Partitioning Tool**  
   Partitions any MATPOWER case into multiple regions using graph Laplacian construction, **spectral clustering**, and **k-means**. It outputs `.m` and `.mat` files and visualizes regions and tie-lines.

2. **Distributed DC OPF Solver**  
   A fully modular **ADMM-based DCOPF solver** implemented in **YALMIP**. Each region solves its local OPF and coordinates across tie-lines using dual variables.

3. **Distributed AC OPF Solver**  
   A distributed **ACOPF solver using ADMM and IPOPT**. It supports scalable, privacy-preserving operation over large-scale power systems partitioned into regions.

The DCOPF and ACOPF solvers are designed to **validate and benchmark** the generated distributed datasets under realistic power system operating conditions.

---

## 📁 Repository Structure

DPLib/
│
├── DPLib toolkit/                   # Graph partitioning + solvers
│   ├── DPLib_partitioning_tool/    # Graph-based spectral clustering tool
│   ├── DPLib_DCOPF_Solver/         # Distributed DC OPF (ADMM, YALMIP)
│   └── DPLib_ACOPF_Solver/         # Distributed AC OPF (ADMM, IPOPT)
│
├── pglib_opf_case118_ieee/         
├── pglib_opf_case1354_pegase/
├── pglib_opf_case162_ieee_dtc/
├── pglib_opf_case179_goc/
├── pglib_opf_case1888_rte/
├── pglib_opf_case1951_rte/
├── pglib_opf_case2000_tamu/
├── pglib_opf_case200_tamu/
├── pglib_opf_case2383wp_k/
├── pglib_opf_case240_pserc/
├── pglib_opf_case2746wp_k/
├── pglib_opf_case2869_pegase/
└── pglib_opf_case300_ieee/
└── ....


Each `pglib_opf_case*` folder contains a distributed benchmark version of a MATPOWER test case partitioned into multiple regions using the DPLib partitioning tool. These datasets are ready for validation using the DCOPF and ACOPF solvers.

---

## 📦 Requirements

- MATLAB R2020a or newer  
- [MATPOWER toolbox](https://matpower.org)  
- [YALMIP toolbox](https://yalmip.github.io/download/)  
- [IPOPT solver](https://coin-or.github.io/Ipopt/) (required for ACOPF)  
- Statistics and Machine Learning Toolbox (for `kmeans` clustering)

Ensure all toolboxes and solvers are installed and added to your MATLAB path.

---

## 📌 Usage

Each tool includes a standalone README inside its folder:

- `DPLib_partitioning_tool/README.md`
- `DPLib_DCOPF_Solver/README.md`
- `DPLib_ACOPF_Solver/README.md`

---

## 📚 Citation

If you use this framework or its solvers in your research, please cite:

> M. Hasanzadeh, "Toward a Standard Repository for Distributed Power System Optimization and Control," 2025.

---

## 👤 Author

**Milad Hasanzadeh**  
Department of Electrical and Computer Engineering  
Louisiana State University  
📧 e.mhasanzadeh1377@yahoo.com  

📅 **Release Date:** June 2025  
📄 **License:** Academic and Research Use Only  

---

## 🔗 Project URL

[https://github.com/LSU-RAISE-LAB/DPLib.git](https://github.com/LSU-RAISE-LAB/DPLib.git)
