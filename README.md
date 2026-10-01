# HBB Variant Effect Analysis

## Question
How do known point mutations in the human beta-globin gene (HBB)
change the protein sequence?

## Method
- Fetched HBB mRNA (NM_000518.5) from NCBI using Biopython (Entrez)
- Extracted the coding sequence and applied known variants
  (HbS, HbC, HbE)
- Translated wild-type and mutant sequences, aligned the proteins
  (BLOSUM62), and classified each effect
- Visualized variant positions in R (ggplot2)

## Results
| Variant | Protein change | Effect | Protein identity |
|---|---|---|---|
| HbS (sickle) c.20A>T | p.E7V | missense | 99.32% |
| HbC c.19G>A | p.E7K | missense | 99.32% |
| HbE c.79G>A | p.E27K | missense | 99.32% |

All three variants are missense mutations that change one amino
acid of the 147-residue beta-globin protein.



![Variant map](results/hbb_variant_map.png)



## Notes
Variants use c. numbering (A of ATG = 1). The sickle mutation is
called Glu6Val in older papers and p.Glu7Val in current notation.

## Tools
Python 3, Biopython, R, ggplot2

## How to run
1. pip install biopython
2. python variant_analysis.py
3. Run plot_variants.R in RStudio