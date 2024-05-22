# Welcome to the Granger Causality in Extremes Package!

For any questions or error reports, please contact [juraj.bodik@unil.ch](mailto:juraj.bodik@unil.ch).

## Overview
This package contains two main functions: `Extreme_causality_test` and `Extreme_causality_full_graph_estimate`. Both functions are defined in the file 'Main functions', without any packages required.

## Extreme_causality_test

This function tests whether the tails/extremes of time series `X` cause those of time series `Y`, given confounders `Z`.

### Function Inputs:
- **x**: A numeric vector representing the first time series (potential cause).
- **y**: A numeric vector representing the second time series (potential effect).
- **z**: A data.frame of potential confounders. Set to NULL if there are no confounders.
- **lag_future**: The time delay for the effect from `x` to `y`. This is the coefficient $p$ in Appendix A of the manuscript.
- **p_value_computation**: If `p_value_computation=FALSE` it returns the output based on (fast) Algorithm 1.  If `p_value_computation=TRUE` it computes the p-value for the hypothesis $H_0: X \text{ does not cause } Y \text{ in extremes given } Z$. If $p\_value < 0.05$, we conclude that $X$ causes $Y$ given $Z$. 
- **bootstrap_repetitions**: The number of bootstrap repetitions for p-value computation. More repetitions yield more precise p-values but require longer computation time.
- **both_tails**: Set to TRUE to consider both large and extremely negative values. For example, in GARCH models, both tails are of interest, while in VAR models, only large values might be relevant.
- **nu_x**: The coefficient $\tau_X$ or $k_n$ in the manuscript, defined as $k = \lfloor n^{\nu_x} \rfloor$. If strong hidden confounding is expected, set $\nu_x$ to 0.4 or 0.5.
- **q_y**: The coefficient $\tau_y = q_y \times n$, describing the conditioning on $Y_t$. For large auto-correlation in $Y$, set $q_y$ to 0.1 or less. Note that in the manuscript, $q_y$ is defined as $1 - q_y$.
- **q_z**: The coefficient $\tau_z = q_z \times n$, describing the conditioning on $Z_t$. This is irrelevant if $Z$ is NULL. For strong confounding effects, set $q_z$ to 0.2 or 0.3.
- **lag_past**: The lag from $Z$ to $(X, Y)$. If the common cause has different lags to $X$ and $Y$, it may cause spurious causality between $X$ and $Y$. Ensure `lag_past` is larger than this lag.

### Function Outputs:
- **output**: Either 'Evidence of causality' or 'No causality' based on Algorithm 1 from the manuscript.
- **p_value_tail**: This is not shown if `p_value_computation` is FALSE. Rejection indicates evidence of causality in the tail. It corresponds to the p-value for the hypothesis $H_0: X \text{ does not cause } Y \text{ in the tail given } Z$, based on bootstrapping. Often `p-value = 1` which means that $\text{CTC} < \text{baseline}$.
- **p_value_extreme**: This is not shown if `p_value_computation` is FALSE. Rejection indicates evidence of causality in extremes. It corresponds to the p-value for the hypothesis $\hat{\Gamma}< (1 +3 \hat{\Gamma}^{baseline})/4$.
- **CTC**: The coefficient $\hat{\Gamma}_{X \rightarrow Y | Z}$.
- **baseline**: The baseline coefficient $\hat{\Gamma}^{\text{baseline}}_{X \rightarrow Y | Z}$.

## Extreme_causality_full_graph_estimate

This function estimates the causal graph (path diagram) between a set of time series.

### Function Inputs:
- **w**: A data.frame of all time series, which should be numeric and of the same length.
- **lag_future**: Same as in `Extreme_causality_test`.
- **both_tails**: Same as in `Extreme_causality_test`.
- **nu_x**: Same as in `Extreme_causality_test`.
- **q_y**: Same as in `Extreme_causality_test`.
- **q_z**: Same as in `Extreme_causality_test`.
- **lag_past**: Same as in `Extreme_causality_test`.

### Function Outputs:
- **G$G**: A graph defined by its edges. Each row corresponds to an edge from the first column pointing to the second column. Use `graph <- graph_from_edgelist(G$G)` from the `igraph` library to obtain the graph environment.
- **G$weights**: Weights corresponding to each edge, representing how close the coefficient $\hat{\Gamma}$ is to 1. If $\hat{\Gamma}_{X \rightarrow Y | Z} = 1$, the weight is 1. The weight is 0 if $\hat{\Gamma} = \left(1 + \hat{\Gamma}^{\text{baseline}}\right) / 2$.

## Example Usage

```r
library(igraph)   # For visualizing the final graph estimates
library(EnvStats) # Or any other package that can generate Pareto noise

# Example: Generating a 4-dimensional VAR time series with lag = 2
n <- 5000
epsilon_x <- rpareto(n, 1, 1)
epsilon_y <- rpareto(n, 1, 1)
epsilon_z1 <- rpareto(n, 1, 1)
epsilon_z2 <- rpareto(n, 1, 1)
x <- rep(0, n)
y <- rep(0, n)
z1 <- rep(0, n)
z2 <- rep(0, n)

for (i in 3:n) {
  z1[i] <- 0.5 * z1[i - 1] + epsilon_z1[i]
  z2[i] <- 0.5 * z2[i - 1] + epsilon_z2[i]
  x[i] <- 0.5 * x[i - 1] + 0.5 * z1[i - 2] + 0.5 * z2[i - 1] + epsilon_x[i]
  y[i] <- 0.5 * y[i - 1] + 0.5 * z1[i - 1] + 0.5 * z2[i - 2] + 0.2 * x[i - 1] + epsilon_y[i]
}

z <- data.frame(z1, z2)
w <- data.frame(z1, z2, x, y)

# Running the extreme causality tests
Extreme_causality_test(x, y, z, lag_future = 2, p_value_computation = FALSE)
Extreme_causality_test(y, x, z, lag_future = 2, p_value_computation = FALSE)

# Estimating the full causality graph
G <- Extreme_causality_full_graph_estimate(w, lag_future = 2)  # Try it out also with lag = 1. You will see that the lagged edges disappear

# Visualizing the graph using igraph package
graph <- graph_from_edgelist(G$G)
V(graph)$name <- names(w)
plot(graph, layout = layout_nicely(graph), vertex.label = V(graph)$name)

