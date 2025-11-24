# SUMC (Spatially Unified Meta-Clustering)

### Chen Di

A computational framework for robust integration of multiple spatial transcriptomics datasets across platforms and batches.

![png](images/workflow.png)

## Overview

SUMC (Spatially Unified Meta-Clustering) is a deep learning-based framework that addresses the critical challenge of batch effects in cross-sample spatial transcriptomics analysis. It employs a three-stage clustering strategy to align spatial architectures and identify conserved spatial patterns across diverse datasets.

## Repository Contents

SUMC.ipynb, provides the core SUMC algorithm implementation, includes a runnable example demonstrating model training

prepareForSUMC.R, Prepares input files required for running SUMC. Processes spatial data to generate feature and similarity matrices.

data/, example of two input files.
