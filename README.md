# Phylogeny Estimation Methods Comparison

This repository contains code and results for comparing different phylogeny estimation methods on simulated datasets. The analysis focuses on evaluating the performance of FastTree (with GTR and JC models) and Neighbor Joining (with JC, LogDet, and p-distance metrics) on the 1000M1 and 1000M4 model conditions from Liu et al. (2009).

## Overview

Phylogenetic tree estimation is a fundamental problem in computational biology. This project compares different methods for their accuracy in reconstructing evolutionary relationships from sequence data. We use simulated datasets with known true trees to precisely measure reconstruction accuracy.

## Datasets

We use the 1000M1 and 1000M4 model conditions from:

Liu, K., Raghavan, S., Nelesen, S., Linder, C. R., & Warnow, T. (2009). "Rapid and accurate large-scale coestimation of sequence alignments and phylogenetic trees." Science, 324(5934), 1561-1564.

Each dataset consists of:
- 1999 sequences per replicate
- Multiple replicates per model condition
- True trees (PIMT) for accuracy assessment
- Sequences evolved under GTR+Gamma model with different parameters for each condition

## Methods Evaluated

1. **FastTree with GTR model**
   - Command: `FastTree -gtr -nt -fastest alignment > output.tree`

2. **FastTree with JC model**
   - Command: `FastTree -nt -fastest alignment > output.tree`

3. **Neighbor Joining with JC distances**
   - Command: `fastme -i alignment.phy -o output.tree -m N -d J -T threads`

4. **Neighbor Joining with LogDet distances**
   - Command: `fastme -i alignment.phy -o output.tree -m N -d L -T threads`

5. **Neighbor Joining with p-distances**
   - Command: `fastme -i alignment.phy -o output.tree -m N -d p -T threads`

## Repository Structure

```
├── scripts/
│   ├── run_phylogeny_analysis.sh      # Main workflow script
│   ├── compare_trees.py               # Tree comparison script using DendroPy
│   └── summarize_results.py           # Result aggregation script
├── data/
│   ├── 1000M1/                        # 1000M1 model condition datasets
│   └── 1000M4/                        # 1000M4 model condition datasets
├── results/
│   ├── summary_results.csv            # Aggregated results
│   ├── 1000M1/                        # Results for 1000M1 condition
│   └── 1000M4/                        # Results for 1000M4 condition
├── figures/
│   ├── fp_rates_by_method.png         # False Positive rates visualization
│   └── fn_rates_by_method.png         # False Negative rates visualization
└── paper/
    └── phylogeny_comparison.tex       # LaTeX paper
```

## Key Findings

1. FastTree consistently outperforms Neighbor Joining methods by a large margin (7-10 times lower error rates)
2. Within FastTree, the GTR model improves reconstruction accuracy by 10-18% compared to JC
3. All Neighbor Joining implementations (with different distance metrics) produced identical results
4. The 1000M4 model condition proved easier for reconstruction than 1000M1, despite having higher evolutionary rates

## Running the Analysis

### Requirements

- FastTree 2.1.11
- FastME 2.1.6
- Python 3.x with the following packages:
  - DendroPy 4.x
  - BioPython
  - NumPy
  - Matplotlib
- Slurm workload manager (for HPC execution)

### Execution

1. Clone this repository:
```
git clone https://github.com/yourusername/phylogeny-estimation-comparison.git
cd phylogeny-estimation-comparison
```

2. Set up the directory structure:
```
mkdir -p data/{1000M1,1000M4} results/{1000M1,1000M4} figures
```

3. Download the datasets from the original source and place them in the appropriate directories.

4. Run the analysis:
```
sbatch scripts/run_phylogeny_analysis.sh
```

5. Generate visualizations:
```
python scripts/plot_results.py
```

## Citation

If you use this code or results in your research, please cite:

```
@article{liu2009rapid,
  title={Rapid and accurate large-scale coestimation of sequence alignments and phylogenetic trees},
  author={Liu, Kevin and Raghavan, Sindhu and Nelesen, Serita and Linder, C Randal and Warnow, Tandy},
  journal={Science},
  volume={324},
  number={5934},
  pages={1561--1564},
  year={2009},
  publisher={American Association for the Advancement of Science}
}

@article{price2010fasttree,
  title={FastTree 2--approximately maximum-likelihood trees for large alignments},
  author={Price, Morgan N and Dehal, Paramvir S and Arkin, Adam P},
  journal={PloS one},
  volume={5},
  number={3},
  pages={e9490},
  year={2010},
  publisher={Public Library of Science}
}

@article{lefort2015fastme,
  title={FastME 2.0: a comprehensive, accurate, and fast distance-based phylogeny inference program},
  author={Lefort, Vincent and Desper, Richard and Gascuel, Olivier},
  journal={Molecular Biology and Evolution},
  volume={32},
  number={10},
  pages={2798--2800},
  year={2015},
  publisher={Oxford University Press}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- The datasets used in this analysis were introduced by Liu et al. (2009)
- This analysis was conducted as part of a computational biology course project
