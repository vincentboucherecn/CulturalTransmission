# Cultural Transmission and Socialization Spillovers in Education
## Vincent Boucher, Carlo L. Del Bello, Fabrizio Panebianco, Thierry Verdier and Yves Zenou

## Replication Files:

### Main Files

+ *allfunction.R*: It's an R script collecting all of the functions called by the other scripts. This script does not have to be launch manually. It is launched by the other scripts.

+ *import.do*: This do-file has to be lauched FIRST. It keeps the relevant variables from the original datasets and merges them in a single file.

+ *children.R*: This Rscript has to be launched SECOND. It maximizes the joint likelihood P(G,s). Estimates are saved in the file "outestim_partial.RData". Produces Table 1--4.

+ *parents_clusters.R*: This Rscript has to be launched THIRD. It performs OLS estimates, with bootstrapped SE for the parents' model. It uses "outestim_partial.RData" as input. The overall results are saved in "outestim_final.RData". Produces Table 5 and Figure 2--3.

+ *policy.R*: This Rscript has to be launched FOURTH. It performs counterfactual policy simulations. It uses "outestim_final.RData" and the results are saved in "total_withpolicy.RData". Produces Figures 4--6.

### Robustness analysis

+ "additional_tables.R": Produces Figures C.1 and C.2

+ "children_alt.R: Comparison with probit model. Produces Table 7.

+ "modelfit.R": Homophily fit. Produces Table 6.

+ "parents_het.R": Finite mixture model. Produces Tables C.1 and C.2 
