# PatHapOuf

This is a R script that helps finding out the most probable number of haplomes per subpopulations present in a sequenced genetic pool. It was initially developped to analyze the content of ant queens spermathecae, that may contain a number of paternal haplomes of various origins. 

The general idea behind this tool is to evaluate the likelihood of all possible combinations of haplome numbers up to a maximum number per subpopulation. Then, one can find out the most probable haplome number combination directly, or look at the most probable haplome number for any subpopulation independently, by summing likelihoods over combinations that support each possible number.

The analysis is meant to be applied to per-allele read counts (i.e. alleles coverage) at a number of polymorphic sites, as can typically be obtained from NGS reads alignments. Importantly, the analysis also needs pre-acquired and independent estimates of allele frequencies in each subpopulation for each site.

# Underlying model
## Data
Assume that per-allele read counts were obtained for $L$ sites and arranged into a set $X$ of size $L$, with each element $X_{i}$ giving per-allele read counts at site $i$. With $A_{i}$ the number of possible alleles at site $i$, that is $X_{i} = \\\{c_{1}, c_{2}, ..., c_{A_{i}}\\\}$. Denote $C_{i} = \sum_{k=1} c_{k}$ the total coverage at site i.

* dataset $X$ is a set with L elements $X_{i}$ 
plop$m_{max}$. For any datased X where G subpopulations are present, that is all $p(X \vert m)$, with $m = \\\{m_{1}, m_{2}, ..., m_{G}\\\}$ and $m_{g} \in [0,m_{max}]$. Then, one can find out the most probable haplome number combination directly, or look at the most probable haplome number for any subpopulation $g$ independently by comparing all $p(X \vert m_{g})$. Any $p(X \vert m_{g} = n)$ is obtained as $\sum{p(X \vert m, m_{g} = n)}$.

The model assumed to generate the observed read counts dataset is the following. First, a true number of haplomes are sampled from the various subpopulation to constitute the genetic pool. At this point, no particular haplome number combination is preferred (i.e., uniform prior on $m$) 

## Error parameter
# Setup
# Input

Read counts and alleles frequencies for each site should be bundled into a single input file (see Input section and example input file).
# Output
