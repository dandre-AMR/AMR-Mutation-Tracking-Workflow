# AMR-Mutation-Tracking-Workflow
A reproducible computational genomics workflow for large-scale extraction, annotation, and temporal analysis of antimicrobial resistance (AMR) genes in Klebsiella pneumoniae and Acinetobacter baumannii, integrating AMRFinderPlus outputs with custom Python scripts for mutation tracking and amino acid analysis.

# Integrated Computational Genomics Workflow for AMR Gene Extraction and Temporal Analysis of ARG Mutations

**Author:** Damani Andre  
**Affiliations:** Morehouse School of Medicine, Bioinformatics Core, and collaborating institutions  
**Year:** 2025  

---

## Overview

This repository contains the complete set of Python scripts, example data, and documentation supporting the paper:

> **Andre D.**, Aslam A., Udoh U., Nwajiobi-Princewill P., Patel A., Ponce R., Weiss D.S., Idris M., Onyedibe K.I.  
> *Integrated Computational Genomics Workflow for Antimicrobial Resistance Gene (ARG) Extraction and Temporal Analysis of ARG Mutations in* *Klebsiella pneumoniae* *and* *Acinetobacter baumannii* (2025).

This workflow automates large-scale antimicrobial resistance (AMR) gene analysis using public genome assemblies. It performs genome retrieval, AMRFinderPlus annotation, sequence extraction, mutation frequency tracking, and amino acid composition analysis to characterize gene-level variability across time.

---

## Key Features

- **Automated genome retrieval:** Downloads *K. pneumoniae* and *A. baumannii* assemblies from NCBI RefSeq and GenBank.
- **Resistance gene annotation:** Batch executes NCBI AMRFinderPlus with full metadata logging.
- **Sequence extraction:** Retrieves annotated gene sequences with strand orientation and coordinates.
- **Mutation frequency tracking:** Identifies unique sequence variants and generates year-by-year binary matrices.
- **Amino acid composition:** Converts nucleotide sequences to codons and tracks amino acid profile changes across time.
- **Reproducible outputs:** Produces structured, machine learning–ready tables for mutation and compositional analysis.

---

## Repository Structure

