# PatHapouf

This is a R script that finds out the most probable number of haplomes per subpopulations present in a sequenced genetic pool. It was initially developped to analyze the content of ant queens spermathecae, that may contain a number of paternal haplomes of various origins. 

## Principle

The general idea of the analysis is to evaluate the likelihood of all possible combinations of haplome numbers up to a maximum number per subpopulation, $m_{max}$. For any datased X where G subpopulations are present, that is all $p(X \vert m)$, with $m = \\\{m_{1}, m_{2}, ..., m_{G}\\\}$ and $m_{g} \in [0,m_{max}]$. Then, one can find out the most probable haplome number combination directly, or look at the most probable haplome number for any subpopulation $g$ independently by comparing all $p(X \vert m_{g})$. Any $p(X \vert m_{g} = n)$ is obtained as $\sum{p(X \vert m, m_{g} = n)}$.

The dataset X is assumed to be read counts for each present allele at a number L of polymorphic sites, as can typically be obtained from NGS reads alignments. 

The analysis also needs pre-acquired and independent estimates of allele frequencies in each subpopulation. This is critical as variance in allele frequency constitutes the signal that is used to tell haplome number combinations apart. 

Read counts and alleles frequencies for each site should be bundled into a single input file (see Input section and example input file).

## Underlying model

The model assumed to generate the observed read counts dataset is the following. 


## Setup
## Input
## Output
