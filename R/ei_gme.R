#' @name ei_gme
#' @title Ecologic Inference applying entropy
#' @description
#' The function ei_gme defines the Shannon entropy function which takes a vector of probabilities as input and returns the negative
#' sum of p times the natural logarithm of p.The function will set the optimization parameters and using the "optim" function an optimal
#' solution is obtained.
#' The function defines the independent variables in the two databases needed, which we call datahp with "n_hp" observations and datahs
#' with "n_hs" observations; and the function of the binary variable of interest y. Then the weights of each observation for the two
#' databases used are defined, if there are no weights available it will be 1.
#' The errors are calculated pondering the support vector of dimension \code{var, 0, -var}. This support vector can be specified by the user.
#' The default support vector is based on variance.We recommend a wider interval with v(-1,0,1) as the maximum.
#' The restrictions are defined to guarantee consistency.
#' The optimization of the Shannon entropy function is solved with the "optim" function local solver "BFGS" and the tolerance by default is settled in
#' 1e-24 but can be specified by the user.The model used in the optimization can be specified too between: "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
#' "Brent". The method by default and the recommended is BFGS
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarize mutate
#' @importFrom stats sd model.matrix model.frame model.response optim terms
#' @importFrom graphics segments
#' @details
#'To solve the optimization upper and lower bounds for p and w are settled, specifically, p and w must be above 0 and lower than 1.
#'In addition, the initial values of p are settled as a uniform distribution and the errors (w) as 1/L.
#'
#' @param fn is the formula that represents the dependent variable in the optimization.
#' In the context of this function, 'fn' is used to define the dependent variable
#' to be optimized by the entropy function.
#' @param datahp The data where the variable of interest y is available and also the independent variables.
#'  Note: The variables and weights used as independent variables must have the same name in 'datahp' and in 'datahs'
#'        The variables in both databases need to match up in content.
#' @param datahs The data with the information of the independent variables as a disaggregated level.
#'        Note: The variables and weights used as independent variables must be the same and must have the same name in 'datahp' and in 'datahs'
#' @param w The weights to be used in this function.
#' @param tol The tolerance to be applied in the optimization function. If the tolerance is not specified, the default tolerance has been set in 1e-24
#' @param v The support vector
#' @param method The method used in the function optim.This can be selected by the user between: "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
#' "Brent". The method by default and the recommended is BFGS
#' @return The function will provide you a dataframe called table with the next information:
#'    \itemize{
#'   \item \strong{probabilities}  Probabilities for each individual to each possibility j of the variable of interest y.
#'   \item \strong{error primal} Errors calculated to the j possibilities of y.
#'   \item \strong{predictions} The prediction for each individual is calculated as the sum of the probability plus the error primal.
#'   The function provides information about the optimization process as :
#'     \item \strong{value_of_entropy} The value of entropy resulting from the optimization.
#'   \item \strong{iterations} Indicates the times the objective function and the gradient has been evaluated during the optimization process
#'   \item \strong{message} Indicates the message if it has been generated in the process of optimization
#'   \item \strong{tol} Indicates the tolerance used in the optimization
#'   \item \strong{method} Indicates the method used in the optimization
#'   \item \strong{v} Indicates the vector of support used in the function
#'   The function provides a dataframe containing the information about lambda:
#'   \item  \strong{lambda} The estimated lambda values.
#'   It is provided an object with the restrictions checked which should be approximately zero.
#'   \item  \strong{check restrictions} Being  g1 the restriction related to the unit probability constraint, g2 to the error unit sum constraint, and g3 to the consistency restriction that implies that the difference between the cross moment in both datasets must be zero.
#'   }
#'   The restriction g3 can be checked thoroughly with the objects by separate.
#'    \itemize{
#'   \item  \strong{cross moments hp} Cross moments in \code{datahp}.
#'   \item  \strong{cross moments hs} Cross moments in \code{datahs}.}
#' @references
#' Fernandez-Vazquez, E., Díaz-Dapena, A., Rubiera-Morollon, F., Viñuela, A., (2020) Spatial Disaggregation of Social Indicators: An Info-Metrics Approach. Social Indicators Research, 152(2), 809–821. https://doi.org/10.1007/s11205-020-02455-z.
#' @examples
#' #In this example we use the data of this package
#' datahp <- financial()
#' datahs <- social()
#' # Setting up our function for the dependent variable.
#' fn               <- datahp$poor_liq ~ Dcollege+Totalincome+Dunemp
#' #Applying the function ei_gme to our databases. In this case datahp
#' #is the data where we have our variable of interest datahs is the data
#' # where we have the information for the disaggregation.
#' #w can be included if we have weights in both surveys
#' #Tolerance in this example is fixed in 1e-20 and v will be (-1,0,1)
#' v=matrix(c(-1, 0, 1), nrow = 1)
#' result  <- ei_gme(fn=fn,datahp=datahp,datahs=datahs,w,tol=1e-20,method="BFGS",v=v)
#' @export
ei_gme <- function(fn,datahp,datahs,w,tol,method,v=NULL){

  #Independent variables
  x_hp             <- model.matrix(fn,datahp)                         #Matrix of independent variables in hp - dimensions n_hs (observations in hp) and k (parameters)

  fn_right         <- fn[-2]                                          #hs does not have dependent variable. It is needed a formula without the left hand side
  x_hs             <- model.matrix(fn_right,datahs)                   #Matrix of independent variables in hs - dimensions n_hp (observations in hs) and k (parameters)

  #loading and rescale of weights (by the total) in each survey
  if ("w" %in% colnames(datahp) && "w" %in% colnames(datahs)) {
    w_hp <- datahp$w / sum(datahp$w)   # vector with n_hp observations
    w_hs <- datahs$w / sum(datahs$w)   # vector with n_hs observations
  } else {
    w_hp <- rep(1, nrow(datahp)) / nrow(datahp)
    w_hs <- rep(1, nrow(datahs)) / nrow(datahs)}

  if (any(is.na(x_hs))) {
    stop("NA (missing) values have been found in the dataset")
  }
  if (any(is.na(x_hp))) {
    stop("NA (missing) values have been found in the dataset")
  }

  #Definition of "y" as a matrix with n_hp observations by j categories.
  form <- model.frame(fn,datahp)
  y          <- model.response(form)                            #First column of matrix y (2 columns in shannon disagg)
  y_prev <- model.matrix(~ factor(y) - 1)

  n_hp             <- dim(x_hp)[1]  ; n_hp                            #observations in hp
  n_hs             <- dim(x_hs)[1]  ; n_hs                            #Observations in hs
  k                <- dim(x_hs)[2]  ; k                               #columns in the matrix of independent variables (including the intercept)
  J                <- ncol(y_prev)    ; J                               #categories or columns in y
  #Support vector for each J

  var<-if (is.numeric(y)) {
    var_result <- var(y, na.rm = TRUE)
  } else if (is.factor(y) || is.character(y)) {
    y_numeric <- as.numeric(factor(y))
    var <- var(y_numeric, na.rm = TRUE)
  } else {
    stop("The variable must be numeric or categorical.")
  }

  if (is.null(v)) {

    v <- matrix(c(-var, 0, var), nrow = 1)
  }
  #Three extreme values (L) to estimate the errors
  l                <- dim(v)[1]                                       #Parameter L is equal to the number of values in vector v
  #Auxiliary matrixes
  L_ones           <- c(1,1,1)                                        #A matrix with L ones.    It is used in the restrictions for w
  N_ones_hs        <- matrix(1,n_hs,1)                                #A matrix with n_hs ones. It is used in the restrictions for p and w
  J_ones           <- matrix(1,J,1)                                   #A matrix with J ones.    It is used in the restrictions for p

  #Define the Shannon entropy function (I multiply it by -1 to maximize it)
  shannon_entropy_dual <- function(lambda_v) {
    lambda         <- matrix(lambda_v, nrow=J ,ncol=, byrow=F)        #Lambda is the parameter to be obtained. Here we transform them into matrix form.
    omega          <- rowSums(t(exp(lambda %*% t(x_hs * w_hs))))     #With J_ones they are multiplied to make the sum in J; t(x_hs * w_hs) this part is Xt(k,n_hs) times the weights of the HS
    psi            <- t(exp(v[1]*lambda %*% t(x_hs * w_hs)) + exp(v[2]*lambda %*% t(x_hs * w_hs)) + exp(v[3]*lambda %*% t(x_hs * w_hs)))#we do the multiplication in three steps because we have to multiply the support vector by the product of lambda and Xt(k,n_hs) by the weights of the HS, CHECK THAT IT IS EQUIVALENT
    #Objective function
    sum(log(omega))  + sum(log(psi)) - sum(t(x_hp) %*%  (y_prev  *  w_hp) * t(lambda))   }


  lambda_ini       <- matrix(0,J,k)
  lambda_v         <- as.vector(lambda_ini)


  if (missing(tol) || is.null(tol)) {
    tol <- 1e-24
  }

  control = list(maxit=1000,reltol=tol)

  if(missing(method)) {
    method <- "BFGS"
  }


  llamar_optim <- function(par, fn, gr = NULL, method, tol, ...) {
    resultado <- stats::optim(par, fn, gr, ...)
    resultado$method <- method  # Guarda el metodo en el resultado
    resultado$tol <- control$reltol  # Guarda la tolerancia en el resultado
    return(resultado)
  }
  res <- llamar_optim(par =lambda_v, fn = shannon_entropy_dual,method=method, control = control) ; res
  tol<-res$tol
  methodused=toString(res$method)
  lambda          <- matrix(res$par, nrow=J ,ncol=, byrow=F)
  #Final estimation of  and psi
  omega       <- rowSums(t(exp(lambda %*% t(x_hs * w_hs))))                                                                          #With J_ones they are multiplied to make the sum in J; t(x_hs * w_hs) this part is Xt(k,n_hs) times the weights of the HS
  psi         <- t(exp(v[1]*lambda %*% t(x_hs * w_hs)) + exp(v[2]*lambda %*% t(x_hs * w_hs)) + exp(v[3]*lambda %*% t(x_hs * w_hs)))  #we do the multiplication in three steps because we have to multiply the support vector by the product of lambda and Xt(k,n_hs) by the weights of the HS, CHECK THAT IT IS EQUIVALENT

  #Estimation of p_dual and errors from lamdbas
  p_dual      <- t(exp(lambda %*% (t(x_hs * w_hs))))/omega
  u_dual_1    <- t(exp(v[1]* lambda %*% t(x_hs * w_hs)))/psi
  u_dual_2    <- t(exp(v[2]* lambda %*% t(x_hs * w_hs)))/psi
  u_dual_3    <- t(exp(v[3]* lambda %*% t(x_hs * w_hs)))/psi
  error_dual      <- u_dual_1 * v[1] + u_dual_2 * v[2] + u_dual_3 * v[3]
  predictions     <- p_dual+ error_dual

  w               <- cbind(u_dual_1,u_dual_2,u_dual_3)

  #Restrictions
  g1              <-  p_dual %*% J_ones - N_ones_hs                        #First  set of restrictions defines the sum of probabilities for each household in hs must be one - n_hs restrictions
  g2              <-  u_dual_1 + u_dual_2 + u_dual_3 - 1                   #Second set of restrictions defines the sum of probabilities for each household in hs must be one in the errors must be one in each j dimension - n_hs * j restricstions

  #Cross moments in the hp. "y" is weighted with "w_hp" by a multiplication element by element (NOT a matrix multiplication)
  cross_moments_hp<- (t(x_hp) %*%  (y_prev  *  w_hp))

  #Cross moments in the hs. "y" is weighted with "w_hs" by a multiplication element by element (NOT a matrix multiplication)
  cross_moments_hs<-  (t(x_hs) %*% ((p_dual * w_hs) + ( error_dual  * w_hs)))
  g3              <- cross_moments_hp - cross_moments_hs   #Last set of restriction defines the cross moments in hp must be equal to hs

  table <- data.frame(datahs$n,predictions, p_dual, error_dual)
  assign_col_names <- function(df, j) {
    n_names<-paste0("n")
    prediction_names <- paste0("predictions_j", 1:j)
    p_names <- paste0("p_dual_j", 1:j)
    e_names <-paste0("e_dual_j", 1:j)
    all<- c("n", prediction_names, p_names,e_names)
    colnames(df) <- all
    return(df)}

  table <- assign_col_names(table, J)


  values<- list(
    entropy=res$value,
    iterations = res$counts,
    message= res$message,
    method=toString(method)
  )


  checkrestrictions<-list(g1=g1,g2=g2,g3=g3)


  table2<-data.frame(lambda)

  generate_output <- function(estimations, values, tol, v, lambda, checkrestrictions, cross_moments_hp, cross_moments_hs, J, fn) {
    output <- structure(list(
      estimations = estimations,
      values = values,
      tol = tol,
      v = v,
      lambda = lambda,
      checkrestrictions = checkrestrictions,
      cross_moments_hp = cross_moments_hp,
      cross_moments_hs = cross_moments_hs,
      J = J,
      fn = fn
    ), class = "shannon")
    return(output)
  }

  output <- generate_output(estimations = table,values =values,tol=tol,v=v, lambda=table2,checkrestrictions=checkrestrictions,cross_moments_hp=cross_moments_hp,cross_moments_hs=cross_moments_hs,J=J,fn=fn)

  class(output) <- "shannon"
  return(output)
}


#' Generate a Plot
#'
#' @description This  function generates a descriptive plot using the result obtained in ei_gme.
#' It illustrates the mean and the confidence interval by disaggregated territorial unit.
#'
#' @param x The output produced by ei_gme
#' @param reg The data column containing the disaggregated territorial units
#' @param ... Additional arguments passed to the plotting function.
#' @return This function provides a graph representing the mean and confidence interval of each disaggregated territorial unit

#' @export
plot.shannon <-  function(x,reg,...){
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("This package requires the 'dplyr' package. Please install and load it.")
  }
  output<-x
  table<-data.frame(reg,output$estimations)
  J   <- output$J

  regmeans<- table %>%
    group_by(reg) %>%
    summarise(across
              (everything(),\(x) mean (x, na.rm = TRUE)))

  ic_lower <- regmeans$reg
  ic_upper <- regmeans$reg

  for (j in 1:J) {
    col_name <- paste0("predictions_j", j)
    regmeans <- regmeans %>%
      dplyr::mutate(
        ic_lower = !!dplyr::sym(col_name) - 1.96 * sd(!!dplyr::sym(col_name)),
        ic_upper = !!dplyr::sym(col_name) + 1.96 * sd(!!dplyr::sym(col_name))
      )
    regmeans <- regmeans %>%
      dplyr::rename_with(~ paste0("ic_lower_", col_name), ic_lower) %>%
      dplyr::rename_with(~ paste0("ic_upper_", col_name), ic_upper)
  }


  ic_lower_col<- paste0("ic_lower_predictions_j", j)
  ic_upper_col <- paste0("ic_upper_predictions_j", j)
   for (j in 1:j) {
    prediction_col <- paste0("predictions_j", j)

    # plot
    plot(regmeans$reg, regmeans[[ prediction_col]],
         type = "p",
         xlab = "Region",
         ylab = prediction_col,
         main = "",
         ylim = c(min(regmeans[, ic_lower_col]), max(regmeans[, ic_upper_col])))



    segments(regmeans$reg, regmeans[[ic_lower_col]],
             regmeans$reg, regmeans[[ic_upper_col]],
             lwd = 2)

    segments(regmeans$reg - 0.1, regmeans[[ic_lower_col]],
             regmeans$reg + 0.1, regmeans[[ic_lower_col]],
             lwd = 2)
    segments(regmeans$reg - 0.1, regmeans[[ic_upper_col]],
             regmeans$reg + 0.1, regmeans[[ic_upper_col]],
             lwd = 2)
  }}

#'Summary
#'
#' @description This function provides a summary of the output obtained with the function ei_gme.
#'
#' @param object The output obtained from ei_gme
#' @param ... Additional arguments passed to the summary function.
#' @return This summary function returns the entropy value and the last iteration in the optimization process.
#'         A dataframe with the means of the estimations for each characteristic j with the predictions the probabilities and the error estimated.
#'         A dataframe with the lambda estimated for each k.
#' \itemize{
#'   \item \code{Iterations}:Indicates the times the objective function and the gradient has been evaluated during the optimization process
#'   \item \code{Entropy value}:The value of entropy resulting from the optimization.
#'   \item \code{mean_estimations}: The predictions, p_dual, and the error for each category j of the variable y
#'   \item \code{lambda}:The estimated lambda values.
#' }
#' @import dplyr
#' @export
summary.shannon <- function(object,...){

  output<-object
  fn=output$fn
  j=output$J

  cat ("Iterations")
  print(output$values$iterations[1])
  cat ("Entropy value")
  print(output$values$entropy[1])
  prediction_col <- paste0("predictions_j", j)
  error_col <- paste0("e_dual_j", 1:j)

  mean_estimations<- output$estimations %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
  print("mean_estimations")
  print(mean_estimations)

  table2<-data.frame( output$lambda)
  colnames(table2) <- c("intercept",attr(terms(fn), "term.labels"))
  rownames(table2)<-( paste0("j", 1:nrow(table2)))
  print ("lambda")
  print(table2)}
