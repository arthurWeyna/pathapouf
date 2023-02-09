# PatHapouf

This is a R script that finds out the most probable number of haplomes per subpopulations present in a sequenced genetic pool. It was initially developped to analyze the content of ant queens spermathecae, that may contain a number of paternal haplomes of various origins. 

## Principle

The general idea of the analysis is to evaluate the likelihood of all possible combinations of haplome numbers up to a maximum number per subpopulation, $m_{max}$. For any datased X where G subpopulations are present, that is all $p(X \vert m)$, with $m = \\\{m_{1}, m_{2}, ..., m_{G}\\\}$ and $m_{g} \in [0,m_{max}]$. Then, one can find out the most probable haplome number combination directly, or look at the most probable haplome number for any subpopulation $g$ independently by comparing all $p(X \vert m_{g})$. Any $p(X \vert m_{g} = n)$ is obatined as $\sum{p(X \vert m)}$

## Setup
## Input
## Output
