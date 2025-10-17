**Bayesian Community Detection in Assortative Stochastic Block Models**
This project implements a Bayesian approach to community detection in networks using Gibbs sampling applied to the Stochastic Block Model (SBM) framework. It introduces an assortative constraint that enforces stronger within-community connections, improving clustering accuracy and robustness.


📘**Overview**

Community detection aims to identify groups of nodes in a network that are more densely connected internally than with other groups.
This project compares:
Standard SBM – no structural constraints.
Assortative SBM – enforces 
𝑃𝑎𝑎>𝑃𝑎𝑏

**A Bayesian inference framework with Gibbs sampling is used to estimate:**

.Community assignments (z)

.Edge probabilities (P)

.Community proportions (θ)

.Assortativity threshold (τ)


🧠**Key Features**

.Bayesian formulation with Dirichlet and Beta priors

.Gibbs sampling algorithm for posterior inference

.Monte Carlo simulations to compare Standard vs. Assortative SBM

.Evaluation with Adjusted Rand Index (ARI)

.Visualizations: boxplots, histograms, τ distributions, and performance scatterplots

📂**Files**

report.pdf – Full project report, including theory, results, and discussion

code.R – Complete R implementation (network generation, Gibbs sampling, evaluation, and visualization)

⚙️**How to Run**

Install R dependencies:

install.packages(c("igraph", "parallel", "ggplot2"))

Run the main R script:

source("code.R")

The script automatically:

Generates synthetic networks

Runs Gibbs sampling for both models

Computes ARI metrics

Produces visual comparisons

**Result Summary:** 

| Model           | k | Mean ARI  | SD(ARI)   | Mean τ | SD(τ) |
| --------------- | - | --------- | --------- | ------ | ----- |
| Standard SBM    | 3 | 0.760     | 0.150     | —      | —     |
| Assortative SBM | 3 | **0.880** | **0.075** | 0.208  | 0.042 |
| Standard SBM    | 4 | 0.350     | 0.190     | —      | —     |
| Assortative SBM | 4 | 0.620     | 0.160     | 0.225  | 0.035 |
