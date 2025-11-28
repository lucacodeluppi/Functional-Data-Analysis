# 📄 Functional and Topological Analysis of Australian Weekly Mortality Data  
*A Case Study using FDA and TDA Techniques*  

This repository contains the full R code and materials used for the analysis presented in:

> **Functional and Topological Analysis of Australian Weekly Mortality Data**  
> *Luca Codeluppi – Department of Economics, Management and Quantitative Methods (DEMM)*  
> *(See paper included in the repo.)*

---

## 🔍 Overview

This project applies **Functional Data Analysis (FDA)** and **Topological Data Analysis (TDA)** to the Australian STMF weekly mortality dataset (2015–2025).  
The goal is to study structural changes across time, compare pre/post-COVID mortality patterns, and explore hidden geometric features in the mortality time series.

The analysis includes:  
✔ smoothing of weekly mortality curves  
✔ age-specific functional PCA  
✔ functional regression (mortality ~ population shares)  
✔ functional ANOVA with permutation test  
✔ sliding-window embeddings  
✔ persistence diagrams, barcodes, and bottleneck distances  

---

## 📊 Methods Used

### **1. Functional Data Analysis (FDA)**
- B-spline smoothing of weekly total deaths  
- Creation of functional data objects (fd)  
- **Functional PCA** on age-mortality profiles  
- **Functional regression** over weeks  
- **Functional ANOVA (FANOVA)** with permutation F-test (global significance test)

### **2. Topological Data Analysis (TDA)**
- Sliding-window embedding of detrended mortality signals  
- Vietoris–Rips complexes  
- Persistent homology (H0 & H1)  
- Persistence diagrams and barcodes  
- Bottleneck distance between early-season vs late-season weeks in 2020  

---

## 📁 Repository Structure
📂 /src
- data_preprocessing.R
- functional_smoothing.R
- fpca.R
- functional_regression.R
- fanova.R
- tda_analysis.R
