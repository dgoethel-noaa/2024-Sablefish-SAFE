#' Title Get OSA residuals
#'
#' @param obs observed matrix (n_years, n_bins)
#' @param pred predicted matrix (n_years, n_bins)
#' @param iss Input sample size
#' @param iter_wt Iterative weights (if any)
#' @param index Index (age bins or length bins)
#' @param drop_bin Which age or length bin to drop
#'
#' @return
#' @export
#'
#' @examples
get_osa_res = function(obs, pred, iss, iter_wt, index, drop_bin) {
  require(compResidual)
  # Do some munging
  years = rownames(pred) # get years
  ess_wt = round(iss * iter_wt, 0) # get ess
  obs1 = round((obs + 0.001) * (ess_wt)) # calculate observations in numbers
  pred = pred + 0.001 # add constant so we can divide finite
  # determine which bin to drop (software drops last bin automatically)
  obs1 <- cbind(obs1[,-drop_bin], obs1[,drop_bin])
  pred1 <- cbind(pred[,-drop_bin], pred[,drop_bin])
  
  # Get OSA here
  osa_res <- NULL
  while (is.null(osa_res) || !is.finite(sum(osa_res))) {
    osa_res <- resMulti(t(obs1), t(pred1)) 
  } # caluclate residual until it does not return Nans
  
  mat <- t(matrix(osa_res, nrow=nrow(osa_res), ncol=ncol(osa_res))) # munge into matrix
  dimnames(mat) <- list(year=years, index=index) # input dimensions 
  reslong <- reshape2::melt(mat, value.name='resid') # reshape
  return(list(reslong))
}

#' Title Get pearson residuals
#'
#' @param obs observed matrix (n_years, n_bins)
#' @param pred predicted matrix (n_years, n_bins)
#' @param iss Input sample size
#' @param iter_wt Iterative weights (if any)
get_pearson_res = function(iter_wt, iss, obs, pred) {
  ess_wt = round(iter_wt * iss, 0)
  # add constants so we can divide finite
  obs = obs + 0.001
  pred = pred + 0.001
  # Calculate Pearson residuals
  pearson_res = (obs-pred)/sqrt(pred* (1 - pred) / ess_wt) 
  pearson_res = reshape2::melt(pearson_res) 
  names(pearson_res) = c("year", "index", "resid")
  return(pearson_res)
}

#' Title Plot residuals
#'
#' @param data Dataframe with columns year, index, resid
#' @param comp_type Composition type 
#' @param res_type Residual type
#' @param ymin lower ylimit on axis
#'
#' @return
#' @export
#'
#' @examples
res_plot = function(data, comp_type, res_type, ymin) {
  require(ggplot2)
  ggplot(data, aes(year, index, size=abs(resid), color=resid>0)) + geom_point() +
    labs(x = "Year", y = comp_type, title = res_type, 
         size = "Absolute Residual", color = "Obs > Pred (Pos Resid)") +
    theme_reg() +
    scale_y_continuous(limits = c(ymin, NA)) +
    guides(color = guide_legend(order = 0),
           size = guide_legend(order = 1))
}