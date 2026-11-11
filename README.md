# AMR-Mutation-Tracking-Workflow

Code supporting:

> Andre D., Aslam A., Udoh U., Nwajiobi-Princewill P., Patel A., Tyson G., Idris M.Y., Onyedibe K.I.
> *Quantifying Antimicrobial Resistance Gene Sequence Stability Across 170,000 Genomes of Klebsiella pneumoniae and Acinetobacter baumannii.*

The workflow retrieves *K. pneumoniae* and *A. baumannii* assemblies from NCBI, annotates resistance genes with AMRFinderPlus, and summarizes the sequence diversity of each gene and allele with the variant proportion (VP = S/N) and the Stability Index (SI).

## Repository structure

```
environment.yml                    Conda environment (AMRFinderPlus 3.12.6, Python packages)
data/stability_tables/             Per-gene stability tables used in the paper
  KP_Stability.xlsx                K. pneumoniae
  AB_Stability.xlsx                A. baumannii
scripts/pipeline/
  fetch_all_kp_ab.py               Download assemblies from NCBI RefSeq/GenBank (Section 2.1)
  adddate.py                       Add assembly dates from the NCBI assembly summaries (Section 2.1)
scripts/stability/
  kp_stability_by_total.py         Per-gene VP and SI tables (Sections 2.4.2 and 2.5)
scripts/aa_composition/
  build_aa_composition_per_assembly.py   Amino acid composition by assembly and year (Figures S1 to S13)
scripts/validation/
  si_validation.py                 SI validation and comparative figures (Section 2.6, Figures 2, 3, 5, S15 to S17, Tables S4 and S5)
```

## Setup

```
conda env create -f environment.yml
conda activate amr-stability
amrfinder -u   # then pin the database to release 2024-03-25 used in the paper
```

AMRFinderPlus was run with default criteria (98% identity and 98% coverage) and NCBI AMR database release 2024-03-25.

## Usage

1. Download genomes (assemblies used in the paper were downloaded on May 12, 2026):

   ```
   python scripts/pipeline/fetch_all_kp_ab.py <basedir> --species KP --threads 8
   python scripts/pipeline/fetch_all_kp_ab.py <basedir> --species AB --threads 8
   ```

   The accession lists are written to `<basedir>/kp_accessions.txt` and `<basedir>/ab_accessions.txt`.

2. Add assembly dates. `adddate.py` reads its paths from the `BASE` variable at the top of the script; edit it to point to your project directory.

3. Build the stability table for each species from the compiled sequence table:

   ```
   python scripts/stability/kp_stability_by_total.py --infile KP_Compiled_With_Sequence_And_Date_CodonsAA.tsv --outfile KP_Stability.xlsx
   ```

   Output columns include N_total, S_variants, p_max, p_tilde_max, Evenness_E, Redundancy_R, Simpson_D and `SI_star_2comp`, which is the SI reported in the paper:

   SI = p~max x (1 - E), where p~max = (x_max + 0.5) / (N + 1) and E is Shannon evenness (E = 0 when S = 1).

   `SI_star` is the three-component form (SI x R) shown for comparison in Tables S1 and S2.

4. Reproduce the validation analyses and comparative figures from the stability tables:

   ```
   python scripts/validation/si_validation.py \
       --kp data/stability_tables/KP_Stability.xlsx \
       --ab data/stability_tables/AB_Stability.xlsx \
       --outdir results/validation
   ```

   This runs in a few minutes and writes the figures, `shared_genes.csv`, `shared_alleles.csv`, `subsampling.csv` and `summary.json` (correlations, subsampling bias by sample size, and counts quoted in Section 3.5).

## Data

Genome sequences and derived data are available on Zenodo: https://zenodo.org/records/17806965

## Contact

Kenneth I. Onyedibe, Morehouse School of Medicine (konyedibe@msm.edu)
