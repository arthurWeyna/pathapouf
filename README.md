# PatHapOuf

This is a R script that helps finding out the most probable number of haplomes per subpopulations present in a sequenced genetic pool. It was initially developped to analyze the content of ant queens spermathecae, that may contain a number of paternal haplomes of various origins. 

The general idea behind this tool is to evaluate the likelihood of all possible combinations of haplome numbers up to a (user-defined) maximum number per subpopulation. Then, one can find out the most probable haplome number combination directly, or look at the most probable haplome number for any subpopulation independently, by summing likelihoods over combinations that support each possible number.

The analysis is meant to be applied to per-allele read counts (i.e. alleles coverage) at a number of polymorphic sites, as can typically be obtained from NGS reads alignments. Importantly, the analysis also needs independent pre-acquired estimates of allele frequencies in each subpopulation.

## Underlying model

Assume that by sequencing a given genetic pool, per-allele read counts were obtained for $L$ sites and arranged into a set $X$, where each element $X_{i}$ gives per-allele read counts at site $i$. Letting $A_{i}$ be the number of possible alleles at site $i$, that is: 

$$ X_{i} = \\\{c_{i1}, c_{i2}, ..., c_{iA_{i}}\\\}, \textrm{ with total coverage } C_{i} = \sum_{k=1}^{A_{i}} c_{ik}$$

Also assume that there exist $G$ subpopulations whose members are susceptible to be represented in the sequenced genetic pool. If allele frequencies at each site and in each subpopulation are known *a priori*, they can be arranged in a set $F$ where each element $F_{i}$ is a $G$ x $A_{i}$ matrix, whose $j^{th}$ row gives allele frequencies in subpopulation $j$ at site $i$.

We want to compute the likelihood of $X$ for any given combination $M$ of haplome number per subpopulations. Letting $m_{max}$ be the maximum number of haplomes allowed per subpopulation, that is:

$$ p(X \vert M) = \prod_{i=1}^{L} p(X_{i} \vert M) \textrm{  with  } M =\\\{m_{1}, m_{2}, ..., m_{G}\\\}\textrm{, } m_{j} \in [0,m_{max}]\textrm{ and } M_{tot} = \sum_{j=1}^{G} m_{j} \in [1,G \times m_{max}] $$

At first, it is simplest to assume that $X_{i}$ follows a multinomial distribution:

$$ X_{i} \sim Multinomial(C_{i}, S_{i})$$

where $S_{i} = \\\{s_{i1}, s_{i2}, ..., s_{iA_{i}}\\\}$ is the vector of realized allele frequencies in the sequenced genetic pool at site $i$. Assuming equal contribution of each haplome to the pool, the vector $S_{i}$ can be written as $A_{i}/M_{tot}$, where $A_{i} = \\\{a_{i1}, a_{i2}, ..., a_{iA_{i}}\\\}$ is a vector whose $k^{th}$ element gives the realized number of haplome that carry allele $k$ within the sequenced genetic pool at site $i$. Because $A_{i}$ has a finite number $A_{i}^{T} = (M_{tot}+1)^{A_{i}}$ of possible values, we can write:

$$ p(X \vert M) = \prod_{i=1}^{L} p(X_{i} \vert M) = \prod_{i=1}^{L}\sum_{l=1}^{A_{i}^{T}}[p(X_{i} \vert A_{i})p(A_{i} \vert M)] $$ 

This is useful because $p(A_{i} \vert M)$ can be computed. Assuming that individuals from every subpopulation have an equal probability to end up in the sequenced genetic pool, $A_{i}$ follows a Poisson-Multinomial distribution whose parameter matrix can be expressed in terms of $M$ and $F_{i}$:

$$ A_{i} \sim PMD(P) \textrm{ with } P=IF_{i} $$

where $I$ is a $M_{tot}$ x $G$ indicator matrix whose element $[n,j]$ is $1$ if the pool's $n^{th}$ haplome originates from subpopulation $j$, and $0$ otherwise.

## Error parameter

Sadly, likelihood computations as presented above are very sensitive to error because even small sequencing or bioinformatic errors can render a dataset impossible under any $M$. To compensate, PatHapOuf assumes that instead of a binomial distribution with parameters $C_{i}$ and $S_{i}$, $X_{i}$ follows a Dirichlet-multinomial distribution with parameters $C_{i}$ and $D_{i}$ where the contrentration parameter $D_{i}$ depends on both $S_{i}$ and on a fixed user-defined error parameter $e$ in the following way:

$$ D_{i} = \frac{S_{i}}{e} + (1-S_{i})^e $$
 
With this setup, $e$ is to be set between 0 and 1. A small positive $e$ will make $D_{i}$ arbitrarily large and will make $D_{i}/sum(D_{i})$ approach $S_{i}$. In other words, a small $e$ will let the user approach the simpler multinomial case with parameter $S_{i}$ (*i.e.* no error). Conversely, a larger $e$ will make all elements of $D{i}$ approach one, which corresponds to a binomial distribution with uniform parameters (*i.e.* no information). 

Empirically, a good value of $e$ for real data seems to be 0.01, but this will ultimately depend on the error rate within one's specific dataset. Users should try different values as the relationship between $e$ and PatHapOuf's result can be quite insightful.

Setting $e=0$ when calling the script will in fact use the smallest respresentable value to avoid division by 0.

# Requirements

PatHapOuf should work on any up-to-date **R** install, but relies on three packages. **matrixStats**, **extraDistr** and **PoissonMultinomial**. Install these before running!

# Input

Read counts and alleles frequencies for each site should be bundled into a single input file (see example input file).

# Output
