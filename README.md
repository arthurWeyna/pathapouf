# PatHapOuf

This is a R script that helps finding out the most probable number of haplomes per subpopulations present in a sequenced genetic pool. It was initially developped to analyze the content of ant queens spermathecae, that may contain a number of paternal haplomes of various origins. 

The general idea behind this tool is to evaluate the likelihood of all possible combinations of haplome numbers up to a (user-defined) maximum number per subpopulation. Then, one can find out the most probable haplome number combination directly, or look at the most probable haplome number for any subpopulation independently, by summing likelihoods over combinations that support each possible number.

The analysis is meant to be applied to per-allele read counts (i.e. alleles coverage) at a number of polymorphic sites, as can typically be obtained from NGS reads alignments. Importantly, the analysis also needs independent pre-acquired estimates of allele frequencies in each subpopulation.

## Underlying model
Assume that by sequencing a given genetic pool, per-allele read counts were obtained for $L$ sites and arranged into a set $X$, where each element $X_{i}$ gives per-allele read counts at site $i$. Letting $A_{i}$ be the number of possible alleles at site $i$, that is: 

$$X_{i} = \\\{c_{1}, c_{2}, ..., c_{A_{i}}\\\}, \textrm{ with total coverage } C_{i} = \sum_{k=1}^{A_{i}} c_{k}$$

Also assume that there exist $G$ subpopulations whose members are susceptible to be represented in the sequenced genetic pool. If allele frequencies at each site and in each subpopulation are known *a priori*, they can be arranged in a set $F$ where each element $F_{i}$ is a $G$ x $A_{i}$ matrix, whose $j^{th}$ row gives allele frequencies in subpopulation $j$ at site $i$.

We want to compute the likelihood of $X$ for any given combination $M$ of haplome number per subpopulations. Letting $m_{max}$ be the maximum number of haplomes allowed per subpopulation, that is:

$$ p(X \vert M) = \prod_{i=1}^{L} p(X_{i} \vert M) \textrm{  with  } M =\\\{m_{1}, m_{2}, ..., m_{G}\\\}; m_{j} \in [0,m_{max}] $$

At first, it is simplest to assume that $X_{i}$ follows a multinomial distribution:

$$ X_{i} \sim Multinomial(C_{i}, S_{i})$$

where $S_{i} = \\\{s_{1}, s_{2}, ..., s_{A_{i}}\\\}$ is the vector of realized allele frequencies in the sequenced genetic pool. Letting $M_{tot} = \sum_{j=1}^{G}$

## Error parameter
# Setup
# Input

Read counts and alleles frequencies for each site should be bundled into a single input file (see Input section and example input file).
# Output
