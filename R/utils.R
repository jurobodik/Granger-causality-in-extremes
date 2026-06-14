## Opasche's utils file, partially shared across projects.

#' Check directory existence
#'
#' @description Checks if the desired directory exists. If not, the desired directory is created.
#'
#' @param dir_name Path to the desired directory, as a string.
#' @param recursive Should elements of the path other than the last be created?
#' If `TRUE`, behaves like the Unix command `mkdir -p`.
#' @param no_warning Whether to cancel the warning issued if a directory is created (bool).
#'
#' @export
#'
#' @examples \dontrun{check_directory("./results/my_new_folder")}
check_directory <- function(dir_name, recursive=TRUE, no_warning=FALSE){
  if (!dir.exists(dir_name)){
    dir.create(dir_name, recursive=recursive)
    if(!no_warning){
      warning(paste0("The following given directory did not exist and was created by 'check_directory': ", dir_name))
    }
  }
}


#' Safe RDS save
#'
#' @description Safe version of [saveRDS()].
#' If the given save path (i.e. `dirname(file_path)`) does not exist, it is created instead of raising an error.
#'
#' @param object R variable or object to save on disk.
#' @param file_path Path and name of the save file, as a string.
#' @param recursive Should elements of the path other than the last be created?
#' If `TRUE`, behaves like the Unix command `mkdir -p`.
#' @param no_warning Whether to cancel the warning issued if a directory is created (bool).
#'
#' @export
#'
#' @examples \dontrun{safe_save_rds(c(1, 2, 8), "./results/my_new_folder/my_vector.rds")}
safe_save_rds <- function(object, file_path, recursive=TRUE, no_warning=FALSE){
  dir_name <- dirname(file_path)
  check_directory(dir_name, recursive=recursive, no_warning=no_warning)
  
  saveRDS(object, file = file_path)
  
}


# ==== Parallel helpers ====

#' Get doFuture operator
#'
#' @param strategy One of `"sequential"` (default), `"multisession"`, `"multicore"`, or `"mixed"`.
#'
#' @return Returns the appropriate operator to use in a [foreach::foreach()] loop.
#' The \code{\link[foreach]{\%do\%}} operator is returned if `strategy=="sequential"`.
#' Otherwise, the \code{\link[foreach]{\%dopar\%}} operator is returned.
#' @export
#' @importFrom foreach %do% %dopar%
#'
#' @examples `%fun%` <- get_doFuture_operator("sequential")
get_doFuture_operator_deprecated <- function(strategy=c("sequential", "multisession", "multicore", "mixed")){
  
  strategy <- match.arg(strategy)
  
  if(strategy == "sequential"){
    return(foreach::`%do%`)
  } else {
    return(foreach::`%dopar%`)
  }
}


#' Set a doFuture execution strategy
#'
#' @param strategy One of `"sequential"` (default), `"multisession"`, `"multicore"`, or `"mixed"`.
#' @param n_workers A positive numeric scalar or a function specifying the maximum number of parallel futures
#' that can be active at the same time before blocking.
#' If a function, it is called without arguments when the future is created and its value is used to configure the workers.
#' The function should return a numeric scalar.
#' Defaults to [future::availableCores()]`-1` if `NULL` (default), with `"multicore"` constraint in the relevant case.
#' Ignored if `strategy=="sequential"`.
#'
#' @return The corresponding [get_doFuture_operator()] operator to use in a [foreach::foreach()] loop.
#' @export
#' @importFrom foreach %do% %dopar%
#' @importFrom future availableCores plan sequential multisession multicore tweak
#' @importFrom doFuture registerDoFuture
#'
#' @examples \dontrun{
#' `%fun%` <- set_doFuture_strategy("multisession", n_workers=3)
#' # perform foreach::foreach loop
#' end_doFuture_strategy()
#' }
set_doFuture_strategy <- function(strategy=c("sequential", "multisession", "multicore", "mixed"),
                                  n_workers=NULL){
  strategy <- match.arg(strategy)
  
  doFuture::registerDoFuture()
  if(strategy == "sequential"){
    future::plan(future::sequential)
    
  } else if (strategy == "multisession"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores() - 1, 1)
    }
    future::plan(future::multisession, workers = n_workers)
    
  } else if (strategy == "multicore"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores(constraints = "multicore") - 1, 1)
    }
    future::plan(future::multicore, workers = n_workers)
    
  } else if (strategy == "mixed"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores() - 1, 1)
    }
    strategy_1 <- future::tweak(future::sequential)
    strategy_2 <- future::tweak(future::multisession, workers = n_workers)
    future::plan(list(strategy_1, strategy_2))
  }
  # return(get_doFuture_operator(strategy))
  return(doFuture::`%dofuture%`)
}


#' End the currently set doFuture strategy
#'
#' @description Resets the default strategy using `future::plan("default")`.
#'
#' @export
#' @importFrom future plan
#'
#' @examples \dontrun{
#' `%fun%` <- set_doFuture_strategy("multisession", n_workers=3)
#' # perform foreach::foreach loop
#' end_doFuture_strategy()
#' }
end_doFuture_strategy <- function(){
  
  future::plan("default")
}


