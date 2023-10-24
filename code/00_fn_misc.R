# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Functions: miscellaneous 

#' Shortcut for dir(..., full.names=T)
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
dirf <- function(...) {
  dir(..., full.names=T)
}










#' Adapted from butcher::axe_env for bayesian(engine='brms')
#'
#' @param x 
#' @param verbose 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
axe_env_bayesian <- function(x, verbose = FALSE, ...) {
  old <- x
  for(i in seq_along(x$fit$actions$model$spec$args)) {
    attr(x$fit$actions$model$spec$args[[i]], ".Environment") <- rlang::base_env()
  }
  for(i in seq_along(x$fit$fit$spec$args)) {
    attr(x$fit$fit$spec$args[[i]], ".Environment") <- rlang::base_env()
  }
  for(i in seq_along(x$fit$fit$spec$method$fit$args)) {
    attr(x$fit$fit$spec$method$fit$args[[i]], ".Environment") <- rlang::base_env()
  }
  attr(x$fit$actions$model$formula, ".Environment") <- rlang::base_env()
  
  return(x)
}







get_intervals <- function(df, y, type="hdci") {
  library(tidybayes)
  
  ci_fun <- switch(type,
                   "hdci"=hdci,
                   "qi"=qi)
  df |>
    summarise(mn=mean({{y}}),
              med=median({{y}}),
              L025=ci_fun({{y}}, 0.95)[,1],
              L05=ci_fun({{y}}, 0.9)[,1],
              L10=ci_fun({{y}}, 0.8)[,1],
              L15=ci_fun({{y}}, 0.7)[,1],
              L20=ci_fun({{y}}, 0.6)[,1],
              L25=ci_fun({{y}}, 0.5)[,1],
              L30=ci_fun({{y}}, 0.4)[,1],
              L35=ci_fun({{y}}, 0.3)[,1],
              L40=ci_fun({{y}}, 0.2)[,1],
              L45=ci_fun({{y}}, 0.1)[,1],
              L55=ci_fun({{y}}, 0.1)[,2],
              L60=ci_fun({{y}}, 0.2)[,2],
              L65=ci_fun({{y}}, 0.3)[,2],
              L70=ci_fun({{y}}, 0.4)[,2],
              L75=ci_fun({{y}}, 0.5)[,2],
              L80=ci_fun({{y}}, 0.6)[,2],
              L85=ci_fun({{y}}, 0.7)[,2],
              L90=ci_fun({{y}}, 0.8)[,2],
              L95=ci_fun({{y}}, 0.9)[,2],
              L975=ci_fun({{y}}, 0.95)[,2])
}





#' Calculate pseudo-R2s
#' 
#' Currently supports McFadden and the Veall-Zimmermann correction of the Aldrich-Nelson pseudo-R2
#'
#' @param dat.df 
#' @param type 
#'
#' @return
#' @export
#'
#' @examples
calc_R2 <- function(dat.df, type="mf", ...) {
  if(type=="mf") {
    return(
      dat.df |>
        group_by(..., model, PCA, covSet) |>
        summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) |>
        group_by(...) |>
        arrange(y, model) |>
        mutate(R2=1 - LL/first(LL)) |>
        select(-LL)
    )
  }
  if(type=="vz") {
    return(
      dat.df |>
        group_by(..., model, PCA, covSet) |>
        summarise(LL=sum(dbinom(alert, 1, prA1, log=T)),
                  N=n(),
                  p=mean(alert)) |>
        group_by(...) |>
        arrange(y, model) |>
        mutate(R2=(-2*(first(LL)-LL))/(-2*(first(LL)-LL)+N) / 
                 ( (-2*(p*log(p)+(1-p)*log(1-p)))/(1-2*(p*log(p)+(1-p)*log(1-p))) )) |>
        select(-N, -p, -LL)
    )
  }
}






#' Omit rows with NAs in any of the columns specified
#' 
#' Similar to na.omit, but only filters based on the columns given.
#'
#' @param df Dataframe
#' @param ... Unquoted column names
#'
#' @return
#' @export
#'
#' @examples
na_omit_col <- function(df, ...) {
  na_rows <- df |>
    ungroup() |> 
    select(any_of(...))|>
    mutate(row_id=row_number()) |>
    na.omit()
  return(df[na_rows$row_id,])
}



