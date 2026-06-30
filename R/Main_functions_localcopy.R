
## Main functions definitions

#' Extreme Granger-causal test
#' 
#' This function tests whether the tails/extremes of a time series `X` cause those of a time series `Y`, given potential confounders `Z`.
#' 
#' @param x A numeric vector representing the first time series (potential cause).
#' @param y A numeric vector representing the second time series (potential effect).
#' @param z A `data.frame` of potential confounders. Set to `NULL` if there are no confounders.
#' @param max_causal_lag The time delay for the effect from `x` to `y`. This is the coefficient 'p' in Appendix A of the manuscript.
#' @param max_confounder_lag The lag from \eqn{Z} to \eqn{(X, Y)}. If the common cause has different lags to \eqn{X} and \eqn{Y}, it may cause spurious causality between \eqn{X} and \eqn{Y}. 
#'   Ensure `max_confounder_lag` is larger than this lag.
#' @param nu_x The coefficient \eqn{\tau_X} or \eqn{k_n} in the manuscript, defined as \eqn{k_n = \lfloor n^{\nu_x} \rfloor}. If strong hidden confounding is expected, set `nu_x` to 0.4 or 0.5.
#' @param q_y The coefficient \eqn{\tau_y = q_y \times n}, describing the conditioning on \eqn{Y_t}. For large auto-correlation in \eqn{Y}, set `q_y` to 0.1 or less. Note that in the manuscript, \eqn{q_y} is defined as `1 - q_y`.
#' @param q_z The coefficient \eqn{\tau_z = q_z \times n}, describing the conditioning on \eqn{Z_t}. This is irrelevant if `z` is `NULL`. 
#' For strong confounding effects, set `q_z` to 0.2 or 0.3. Note that in the manuscript, \eqn{q_z} is defined as `1 - q_z`.
#' @param both_tails Set to `TRUE` to consider both large and extremely negative values. 
#'   For example, in GARCH models, both tails are of interest, while in VAR models, only large values might be relevant.
#' @param instant Whether instantaneous effects should be captured; defaults to `FALSE`.
#' @param p_value_computation If set to `FALSE` the faster "Algorithm 1" is used. 
#'   If set to `TRUE` the p-value for the hypothesis \eqn{H_0: X \text{ does not cause } Y \text{ in extremes given } Z} is computed. 
#'   If `p_value < 0.05`, we conclude that \eqn{X} causes \eqn{Y} given \eqn{Z}. 
#' @param bootstrap_repetitions The number of bootstrap repetitions for p-value computation. More repetitions yield more precise p-values but require longer computation time.
#' @param choice_of_F Choice of F in the coefficient. Leave default unless you want to reproduce the results from the manuscript
#' 
#' @returns A named list containing:
#' \item{is_causal}{A logical value indicating whether evidence of causality is detected;}
#' \item{output}{Either 'Evidence of causality' or 'No causality' based on Algorithm 1 from the manuscript;}
#' \item{p_value_tail}{This is not shown if `p_value_computation==FALSE`. Rejection indicates evidence of causality in tail. 
#'                     It corresponds to the p-value for the hypothesis H_0: "X does not cause Y in tail given Z", based on bootstrapping. Often `p_value==1` which means that `CTC<baseline`;}
#' \item{p_value_extreme}{This is not shown if `p_value_computation==FALSE`. Rejection indicates evidence of causality in extremes. 
#'                        It corresponds to the p-value for the hypothesis H_0: \eqn{\hat{\Gamma}_{X\rightarrow Y | Z} < (1 + 3 \cdot \hat{\Gamma}^{baseline}_{X\rightarrow Y | Z}) / 4};}
#' \item{CTC}{The coefficient \eqn{\hat{\Gamma}_{X\rightarrow Y | Z}};}
#' \item{baseline}{The baseline coefficient \eqn{\hat{\Gamma}^{baseline}_{X\rightarrow Y | Z}}.}
#' 
#' @importFrom stats median quantile ecdf
#' @export
Extreme_causality_test <- function(x, y, z = NULL, 
                                   max_causal_lag = 1, max_confounder_lag = 0, 
                                   nu_x = 0.3, q_y = 0.2, q_z = 0.1, 
                                   both_tails = TRUE, instant = FALSE, 
                                   p_value_computation = FALSE, bootstrap_repetitions = 50, choice_of_F = 0.5) {
  n <- length(x)
  z <- data.frame(z)
  d <- ncol(z)
  if (d >= 2) {
    q_z <- (2 / d) * q_z
  }
  tau_y <- q_y * n
  tau_z <- q_z * n
  if (both_tails) {
    if (!all(x >= 0)) x <- abs(x - stats::median(x))
    if (!all(y >= 0)) y <- abs(y - stats::median(y))
    if (d >= 1 & !all(z >= 0)) z <- abs(z - apply(z, 2, stats::median))
  }


  F_u <- function(imput, output) { # This is the F^{truc}_Y(t) function that we opt for
    ifelse(output < stats::quantile(imput, choice_of_F), 0, (stats::ecdf(imput)(output)))
  }


  CTC_baseline <- function(x, y, z = NULL, max_causal_lag = max_causal_lag, max_confounder_lag = max_confounder_lag, tau_y = tau_y, tau_z = tau_z, instant = instant) {
    x_to_y_masking <- c()
    for (i in 1:(max_confounder_lag + 1)) {
      yy <- c(rep(0, i - 1), y)
      x_to_y_masking <- c(x_to_y_masking, which(yy > sort(y)[(n - tau_y)]))

      if (d > 0) {
        for (j in 1:d) {
          zz <- c(rep(0, i - 1), z[, j])
          x_to_y_masking <- c(x_to_y_masking, which(zz > sort(z[, j])[(n - tau_z)]))
        }
      }
      x_to_y_masking <- unique(x_to_y_masking)
    }

    future_y <- run_future_max(y, max_causal_lag, instant = instant)
    # baseline1 =mean(F_u(y[-x_to_y_masking], future_y[-x_to_y_masking]))
    baseline2 <- mean(F_u(y, future_y[-x_to_y_masking]))

    return(baseline2)
  }


  CTC_masked <- function(x, y, z = NULL, max_causal_lag, max_confounder_lag, tau_y, tau_z, instant) {
    x_to_y_masking <- c()
    for (i in 1:(max_confounder_lag + 1)) {
      yy <- c(rep(0, i - 1), y)
      x_to_y_masking <- c(x_to_y_masking, which(yy > sort(y)[(n - tau_y)]))

      if (d > 0) {
        for (j in 1:d) {
          zz <- c(rep(0, i - 1), z[, j])
          x_to_y_masking <- c(x_to_y_masking, which(zz > sort(z[, j])[(n - tau_z)]))
        }
      }
      x_to_y_masking <- unique(x_to_y_masking)
    }

    n <- length(x)
    k <- round((n - length(x_to_y_masking))^(nu_x))

    new_x <- x[-x_to_y_masking]

    future_y <- run_future_max(y, max_causal_lag, instant = instant)

    future_y <- future_y[-x_to_y_masking]

    top_x <- which(new_x >= sort(new_x)[length(new_x) - k + 1])


    # return( mean(F_u(y[-x_to_y_masking], future_y[top_x]))  )
    return(mean(F_u(y, future_y[top_x])))
  }

  baseline <- CTC_baseline(x, y, z = z, max_causal_lag = max_causal_lag, max_confounder_lag = max_confounder_lag, tau_y = tau_y, tau_z = tau_z, instant = instant)
  CTC <- CTC_masked(x, y, z, max_causal_lag = max_causal_lag, max_confounder_lag = max_confounder_lag, tau_y = tau_y, tau_z = tau_z, instant = instant)


  if (!p_value_computation) {
    if (CTC <= (1 + baseline) / 2) {
      return(list(is_causal = FALSE, output = "No causality", CTC = CTC, baseline = baseline))
    }
    if (CTC > (1 + baseline) / 2) {
      return(list(is_causal = TRUE, output = "Evidence of causality", CTC = CTC, baseline = baseline))
    }
  }


  if (CTC <= baseline) {
    result <- 0
    result2 <- 1
  }
  if (CTC > baseline) {
    result <- c()
    result2 <- c()

    for (i in 1:bootstrap_repetitions) {
      if (d == 0) original_data <- data.frame(x, y) else original_data <- data.frame(x, y, z)
      new_time_series <- switcher_for_test(original_data, number_of_blocks = floor(sqrt(length(x))))
      tilde_x <- new_time_series[, 1]
      tilde_y <- new_time_series[, 2]
      if (d == 0) {
        tilde_z <- NULL
      }
      if (d > 0) {
        tilde_z <- data.frame(new_time_series[, 3:(d + 2)])
      }
      result <- c(result, CTC_masked(tilde_x, tilde_y, tilde_z, max_causal_lag = max_causal_lag, max_confounder_lag = max_confounder_lag, tau_y = tau_y, tau_z = tau_z, instant = instant))
      result2 <- c(result2, CTC_baseline(tilde_x, tilde_y, tilde_z, max_causal_lag = max_causal_lag, max_confounder_lag = max_confounder_lag, tau_y = tau_y, tau_z = tau_z, instant = instant))
    }
  }


  if (CTC > (1 + baseline) / 2) {
    is_causal <- TRUE
    output <- "Evidence of causality"
  } else {
    is_causal <- FALSE
    output <- "No causality"
  }
  return(list(
    is_causal = is_causal,
    output = output, 
    p_value_tail = mean(result <= result2), 
    p_value_extreme = mean(result <= (1 + 3 * result2) / 4), 
    CTC = CTC, 
    baseline = baseline
  ))
}



#' Extreme Granger-causal graph estimate
#' 
#' This function estimates the causal graph (path diagram) between a set of time series.
#' Use [Extreme_causality_graph_parallel()] instead, for faster parallel computations.
#'
#' @param w A `data.frame` of all time series, which should be numeric and of the same length.
#' @inheritParams Extreme_causality_test
#' @param p_value_based If `FALSE`, Algorithm 1 is used for inferring the edges. 
#'   If `TRUE`, the testing procedure with a cut-off p-value of `p_value_cutoff` is used for detecting the presence of an edge. 
#'   These procedures typically output similar results, but the testing procedure is significantly slower. 
#' @param p_value_cutoff P-value cut-off level to reject the absence of an edge in the estimated graph.
#'
#' @returns A named list containing:
#' \item{G}{A graph defined by its edges. Each row corresponds to an edge from the first column pointing to the second column. 
#'          Use `graph <- graph_from_edgelist(G$G)` from the [igraph] package to obtain the graph environment;}
#' \item{weights}{Weights corresponding to each edge, representing how close the coefficient \eqn{\hat{\Gamma}_{X\rightarrow Y | Z}} is to 1. 
#'                If \eqn{\hat{\Gamma}_{X\rightarrow Y | Z} = 1}, the weight is 1. 
#'                The weight is 0 if \eqn{\hat{\Gamma}_{X\rightarrow Y | Z} = (1 + \hat{\Gamma}^{baseline}_{X\rightarrow Y | Z}) / 2}.}
#' @export
Extreme_causality_graph <- function(w, max_causal_lag = 1, max_confounder_lag = 0, 
                                                  nu_x = 0.3, q_y = 0.2, q_z = 0.1, instant = FALSE, 
                                                  both_tails = TRUE, p_value_based = FALSE, p_value_cutoff = 0.05) {
  m <- ncol(w)

  # Step 1: Pairwise
  G <- c()
  for (i in 1:m) {
    for (j in (1:m)[-i]) {
      x <- w[, i]
      y <- w[, j]
      CTC <- Extreme_causality_test(x, y, z = NULL, 
                                    nu_x = nu_x, q_y = q_y, q_z = q_z, 
                                    max_causal_lag = max_causal_lag, instant = instant, max_confounder_lag = max_confounder_lag, 
                                    both_tails = both_tails, p_value_computation = p_value_based, bootstrap_repetitions = 5 / p_value_cutoff)
      if (p_value_based) {
        if (CTC$p_value_tail <= p_value_cutoff) G <- rbind(G, c(i, j))
      } else {
        # if (CTC$output == "Evidence of causality") G <- rbind(G, c(i, j))
        if (CTC$is_causal) G <- rbind(G, c(i, j))
      }
    }
  }

  # Step 2: Multivariate
  if (all(G == FALSE)) {
    return("Result: Empty graph")
  } else { # if G is non-empty
    
    indexes_to_erase <- .Machine$integer.max
    edges_weights <- c()
    for (i in 1:nrow(G)) {
      x <- w[, G[i, 1]]
      y <- w[, G[i, 2]]

      z_indexes <- intersect(find_parents(G, G[i, 1]), find_parents(G, G[i, 2]))
      if (all(z_indexes == FALSE)) {
        z <- NULL
      } else {
        z <- data.frame(w[, z_indexes])
      }

      CTC <- Extreme_causality_test(x, y, z = z, 
                                    nu_x = nu_x, q_y = q_y, q_z = q_z, 
                                    max_causal_lag = max_causal_lag, instant = instant, max_confounder_lag = max_confounder_lag, 
                                    both_tails = both_tails, p_value_computation = p_value_based, bootstrap_repetitions = 5 / p_value_cutoff)
      if (p_value_based) {
        if (CTC$p_value_tail > p_value_cutoff) {
          indexes_to_erase <- c(indexes_to_erase, i)
        } else {
          edges_weights <- c(
            edges_weights,
            compute_edge_weights_for_p_values(p_value = CTC$p_value_tail, p_value_cutoff = p_value_cutoff)
          )
        }
      } else {
        # if (CTC$output == "No causality") {
        if (!CTC$is_causal) {
          indexes_to_erase <- c(indexes_to_erase, i)
        } else {
          edges_weights <- c(
            edges_weights,
            compute_edge_weights(CTC$CTC, CTC$baseline)
          )
        }
      }
    }
  }

  return(list(G = G[-indexes_to_erase, , drop=FALSE], weights = edges_weights))
}


#' Extreme Granger-causal graph parallelized estimate
#' 
#' This function estimates the causal graph (path diagram) between a set of time series, 
#' using `doFuture` parallellization for computational efficiency.
#' It does the same as [Extreme_causality_graph()], but potentially faster if appropriate parallelization is used.
#'
#' @inheritParams Extreme_causality_graph
#' @param strategy One of `"sequential"` (default), `"multisession"`, `"multicore"`, or `"mixed"`.
#' @param n_workers A positive numeric scalar or a function specifying the maximum number of parallel futures
#' that can be active at the same time before blocking.
#' If a function, it is called without arguments when the future is created and its value is used to configure the workers.
#' The function should return a numeric scalar.
#' Defaults to [future::availableCores()]`-1` if `NULL` (default), with `"multicore"` constraint in the relevant case.
#' Ignored if `strategy=="sequential"`.
#'
#' @inherit Extreme_causality_graph return
#' @export
Extreme_causality_graph_parallel <- function(w, max_causal_lag = 1, max_confounder_lag = 0, 
                                                      nu_x = 0.3, q_y = 0.2, q_z = 0.1, instant = FALSE, 
                                                      both_tails = TRUE, p_value_based = FALSE, p_value_cutoff = 0.05,
                                                      strategy = c("sequential", "multisession", "multicore", "mixed"), n_workers = NULL) {
  m <- ncol(w)
  strategy <- match.arg(strategy)
  `%fun%` <- set_doFuture_strategy(strategy, n_workers)

  # Step 1: Pairwise
  index_pairs <- expand.grid('j' = 1:m, 'i' = 1:m)
  index_pairs <- index_pairs[index_pairs$i != index_pairs$j, c('i', 'j')]
  
  idx <- NULL
  G <- foreach::foreach(idx = 1:nrow(index_pairs), .combine = rbind, .options.future = list(seed = p_value_based)) %fun% {
    i <- index_pairs$i[idx]
    j <- index_pairs$j[idx]
    x <- w[, i]
    y <- w[, j]
    CTC <- Extreme_causality_test(x, y, z = NULL, 
                                  nu_x = nu_x, q_y = q_y, q_z = q_z, 
                                  max_causal_lag = max_causal_lag, instant = instant, max_confounder_lag = max_confounder_lag, 
                                  both_tails = both_tails, p_value_computation = p_value_based, bootstrap_repetitions = 5 / p_value_cutoff)
    if (p_value_based) {
      if (CTC$p_value_tail <= p_value_cutoff) {
        tmp_toadd <- c(i, j)
      } else {
        tmp_toadd <- NULL
      }
    } else {
      # if (CTC$output == "Evidence of causality") {
      if (CTC$is_causal) {
        tmp_toadd <- c(i, j)
      } else {
        tmp_toadd <- NULL
      }
    }
    tmp_toadd
  }

  # Step 2: Multivariate
  if (all(G == FALSE)) {
    end_doFuture_strategy()
    return("Result: Empty graph")
  } else { # if G is non-empty
    
    tmp_res_multiv <- foreach::foreach(i = 1:nrow(G), .options.future = list(seed = p_value_based)) %fun% {
      x <- w[, G[i, 1]]
      y <- w[, G[i, 2]]

      z_indexes <- intersect(find_parents(G, G[i, 1]), find_parents(G, G[i, 2]))
      if (all(z_indexes == FALSE)) {
        z <- NULL
      } else {
        z <- data.frame(w[, z_indexes])
      }

      CTC <- Extreme_causality_test(x, y, z = z, 
                                    nu_x = nu_x, q_y = q_y, q_z = q_z, 
                                    max_causal_lag = max_causal_lag, instant = instant, max_confounder_lag = max_confounder_lag, 
                                    both_tails = both_tails, p_value_computation = p_value_based, bootstrap_repetitions = 5 / p_value_cutoff)
      if (p_value_based) {
        if (CTC$p_value_tail > p_value_cutoff) {
          tmp_feres <- list(index = i, weight = NA)
        } else {
          tmp_feres <- list(index = NA, weight = compute_edge_weights_for_p_values(p_value = CTC$p_value_tail, p_value_cutoff = p_value_cutoff))
        }
      } else {
        # if (CTC$output == "No causality") {
        if (!CTC$is_causal) {
          tmp_feres <- list(index = i, weight = NA)
        } else {
          tmp_feres <- list(index = NA, weight = compute_edge_weights(CTC$CTC, CTC$baseline))
        }
      }
      tmp_feres
    }
    
    indexes_to_erase <- c(unlist(lapply(tmp_res_multiv, function(r) r$index)), .Machine$integer.max)
    indexes_to_erase <- indexes_to_erase[!is.na(indexes_to_erase)]
    edges_weights <- unlist(lapply(tmp_res_multiv, function(r) r$weight))
    edges_weights <- edges_weights[!is.na(edges_weights)]
    
  }
  end_doFuture_strategy()

  return(list(G = G[-indexes_to_erase, , drop=FALSE], weights = edges_weights))
}


## Helper functions

#' run_future_max
#'
#' @param x .
#' @param k .
#' @param instant .
#'
#' @returns .
#'
#' @keywords internal
run_future_max <- function(x, k = 3, instant = TRUE) { # instant= do we want to consider Y_0?
  if (instant) {
    q <- 0
  } else {
    q <- 1
  }
  n <- length(x)
  y <- c()
  for (i in 1:n) {
    y <- rbind(y, max(x[(min(n, i + q)):(min(n, i + k))]))
  }
  return(y)
}

#' switcher_for_test
#'
#' @param x .
#' @param number_of_blocks .
#'
#' @returns .
#'
#' @keywords internal
switcher_for_test <- function(x, number_of_blocks = 15) { # resampling block-wise
  # x must be a matrix (rows=time, columns=variables)

  n <- nrow(x)
  m <- n %/% number_of_blocks # length of one block

  y <- c()
  for (i in 1:number_of_blocks) { # choose one random block with the beginning uniformly chosen from 1:(n-m)
    kocka <- sample(1:(n - m), 1)
    for (j in 1:m) {
      y <- rbind(y, x[kocka + j, ]) # Add this block to our resampled series
    }
  }
  # what to do with the ending if it is not divisible? we just add one random block with a smaller length to obtain times eris with length $n$ again
  if (ncol(x) == 1) { # Code for one-dimensional time series
    if (n %% number_of_blocks != 0) {
      y <- c(y, x[((number_of_blocks * m + 1):n), ])
    }
  } else {
    if (n %% number_of_blocks != 0) {
      k <- n - number_of_blocks * m
      kocka <- sample(1:(n - k), 1)

      for (j in 1:k) {
        y <- rbind(y, x[kocka + j, ])
      }
    }
  }
  return(data.frame(y))
}


#' find_parents
#'
#' @param G .
#' @param vertex .
#'
#' @returns Parents
#'
#' @keywords internal
find_parents <- function(G, vertex) {
  result <- c()
  for (i in 1:nrow(G)) {
    if (G[i, 2] == vertex) {
      result <- c(result, G[i, 1])
    }
  }
  return(result)
}

#' compute_edge_weights
#'
#' @param CTC .
#' @param baseline .
#'
#' @returns Edge weights.
#'
#' @keywords internal
compute_edge_weights <- function(CTC, baseline) {
  (CTC - ((1 + baseline) / 2)) / (1 - ((1 + baseline) / 2))
}

#' compute_edge_weights_for_p_values
#'
#' @param p_value .
#' @param p_value_cutoff .
#'
#' @returns Edge weights.
#'
#' @keywords internal
compute_edge_weights_for_p_values <- function(p_value, p_value_cutoff) {
  (1 / p_value_cutoff) * (p_value_cutoff - p_value)
}

# TODO: un-nest other functions (F_u, CTC_baseline, CTC_masked)? (They are not self-contained.)
