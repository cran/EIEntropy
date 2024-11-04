#'@name ei_gce
#'@title Ecologic Inference applying entropy
#'@description The function ei_gce defines the Kullback-Leibler function which
#'minimises the distance between the distribution of probabilities P and the
#'distribution Q. The distribution Q is based on prior information that we have
#'of our variable of interest previous to the analysis. The function will set
#'the optimization parameters and, using the "nlminb"
#'function, an optimal solution is obtained.
#'The function defines the independent variables in the two databases
#'needed, which we call datahp with "n_hp" observations and datahs with "n_hs"
#'observations; and the function of the variable of interest y. Then the
#'weights of each observation for the two databases used are defined, if there
#'are not weights available it will be 1 by default. The errors are calculated
#'pondering the support vector of dimension \code{var, 0, -var}. This support vector
#'can be specified by the user. The default support vector is based on variance.
#'We recommend a wider interval with v(1,0,-1) as the maximum.
#'The restrictions are defined in order to guarantee consistency. The
#'minimization of Kullback_Leibler distance is solved with "nlminb" function
#'with maximum number of iterations 1000 and with tolerance
#'defined by the user. If the user did not define tolerance it will be 1e-10 by
#'default. For additional details about the methodology see Fernández-Vazquez, et al. (2020)
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by summarize mutate
#' @importFrom stats sd model.matrix model.frame model.response nlminb
#' @details
#'To solve the optimization upper and lower bounds for p and w are settled, specifically, p and w must be above 0 and lower than 1.
#'In addition, the initial values of p are settled as the defined prior and the errors (w) as 1/L.
#'@param fn is the formula that represents the dependent variable in the optimization.
#' In the context of this function, 'fn' is used to define the dependent variable
#' to be optimized by the Kullback-Leibler divergence function.
#' Note: If the dependent variable is categorical the sorting criterion for the columns, and therefore for J, is alphabetical order.
#'@param datahp The data where the variable of interest y is available and also the independent variables.
#'  Note: The variables and weights used as independent variables must have the same name in 'datahp' and in 'datahs'
#'@param datahs The data with the information of the independent variables as a
#'  disaggregated level.
#'  Note: The variables and weights used as independent variables must have the same name in 'datahp' and in 'datahs'. The variables in both databases need to match up in content.
#'@param weights A character string specifying the column name to be used as weights in both `datahp` and `datahs` datasets.
#' If `weights` is provided and exists in both datasets, each dataset's weights will be normalized by the sum of the weights in that dataset.
#' If `weights` is NULL or the specified column does not exist in both datasets, equal weights are applied across all observations.
#'@param tol The tolerance to be applied in the optimization function. If the tolerance is not specified, the default tolerance has been set in 1e-10
#'@param q The prior distribution Q
#'@param v The support vector
#'@param iter The maximum number of iterations allowed for the optimization algorithm to run 
#' Increasing the number of iterations may improve the likelihood of finding an optimal solution, 
#' but can also increases computation time.If the maximum number of iterations is not specified, it will default to 1000
#' @return The function will provide you a dataframe called table with the next information:
#'    \itemize{
#'   \item \strong{probabilities}  Probabilities for each individual to each possibility \code{j} of the variable of interest \code{y}.
#'   \item \strong{error dual}  Errors calculated to the \code{j} possibilities of \code{y}.
#'   \item \strong{predictions}  The prediction for each individual is calculated as the sum of the probability plus the error.
#'   The function provides information about the optimization process as:
#'   \item \strong{divergencekl}The Kullback-Leibler divergence value resulting from the optimization.
#'   \item \strong{iterations} Indicates the times the objective function and the gradient has been evaluated during the optimization process,if any.
#'   \item \strong{message} Indicates the message if it has been generated in the process of optimization.
#'   \item \strong{tol} Indicates the tolerance of the optimization process.
#'   \item \strong{v} Indicates the support vector used in the function.
#'   The function provides a dataframe containing the information about lambda:
#'   \item  \strong{lambda} The estimated lambda values.
#'   It is provided an object with the restrictions checked which should be zero.
#'   \item  \strong{check restrictions} Being  \code{g1}  the restriction related to the unit probability constraint, \code{g2}  to the error unit sum constraint, and \code{g3}  to the consistency restriction that implies that the difference between the cross moment in both datasets must be zero.
#'   The restriction g3 can be checked thoroughly with the objects by separate.
#'   \item  \strong{cross moments hp} Cross moments in \code{datahp}.
#'   \item  \strong{cross moments hs} Cross moments in \code{datahs}.
#'   }
#' @references
#' Fernandez-Vazquez, E., Diaz-Dapena, A., Rubiera-Morollon, F., Viñuela, A., (2020) Spatial Disaggregation of Social Indicators: An Info-Metrics Approach. Social Indicators Research, 152(2), 809–821. https://doi.org/10.1007/s11205-020-02455-z.
#' @examples
#' #In this example we use the data of this package
#' datahp <- financial()
#' datahs <- social()
#' # Setting up our function for the dependent variable.
#' fn               <- datahp$poor_liq ~ Dcollege+Totalincome+Dunemp
#' #In this case we know that the mean probability of being poor is 0.35.With this function
#' #we can add the information as information a priori. This information a priori correspond to the
#' #Q distribution and in this function is called q for the sake of simplicity:
#' q<- c(0.5,0.5)
#' v<- matrix(c(0.2,0,-0.2))
#' #Applying the function ei_gce to our databases. In this case datahp is the
#' # data where we have our variable of interest
#' #datahs is the data where we have the information for the disaggregation.
#' #w can be included if we have weights in both surveys
#' result  <- ei_gce(fn,datahp,datahs,q=q,weights="w",v=v)
#' @export
ei_gce<- function(fn,datahp,datahs,q,weights=NULL,v,tol,iter){

  x_hp             <- model.matrix(fn,datahp)                         #Matrix of independent variables in hp - dimensions n_hs (observations in hp) and k (parameters)
  
  fn_right         <- fn[-2]                                          #hs does not have dependent variable. It is needed a formula without the left hand side
  x_hs             <- model.matrix(fn_right,datahs)                   #Matrix of independent variables in hs - dimensions n_hp (observations in hs) and k (parameters)
  
  #loading and rescale of weights (by the total) in each survey
  if (!is.null(weights) && weights %in% colnames(datahp) && weights %in% colnames(datahs)) {
    # Use the specified column for weights, normalizing by the sum of weights
    w_hp <- datahp[[weights]] / sum(datahp[[weights]], na.rm = TRUE)
    w_hs <- datahs[[weights]] / sum(datahs[[weights]], na.rm = TRUE)
  } else {
    # If no weights or column is missing, apply equal weights
    w_hp <- rep(1 / nrow(datahp), nrow(datahp))
    w_hs <- rep(1 / nrow(datahs), nrow(datahs))
  }
  if (any(is.na(x_hs))) {
    stop("NA (missing) values have been found in the dataset")
  }
  if (any(is.na(x_hp))) {
    stop("NA (missing) values have been found in the dataset")
  }
  
  #Definition of "y" as a matrix with n_hp observations by j categories.
  form <- model.frame(fn,datahp)
  y          <- model.response(form) 
  if (length(unique(y)) == 2 && all(sort(unique(y)) == c(0, 1))) {
    y_factor <- factor(y, levels = c(1, 0))  
  } else {
    y_factor <- factor(y)
  }
  y_prev <- model.matrix(~ y_factor - 1)
  

  # Set the parameters
  n_hp <-dim(x_hp)[1]  ; n_hp                 #observations in hp
  n_hs <- dim(x_hs)[1]  ; n_hs                 #Observations in hs
  k    <- dim(x_hs)[2]  ; k                    #columns in the matrix of independent variables (including the intercept)
  J    <- ncol(y_prev)  ; J                    #categories or columns in y


var<-if (is.numeric(y)) {
      var_result <- var(y, na.rm = TRUE)
    } else if (is.factor(y) || is.character(y)) {
      y_numeric <- as.numeric(factor(y))
      var <- var(y_numeric, na.rm = TRUE)
    } else {
      stop("The variable must be numeric or categorical.")
    }
if (missing(v)|| is.null(v))  {
  v <- matrix(c(var, 0, -var), nrow = 1)
  }else {
    if (length(v) != 3 || nrow(v) != 1) {
      v <- matrix(as.vector(v)[1:3], nrow = 1)
      warning("The matrix `v` was automatically reshaped to dimensions `1 x 3`.")
    }
  }

  #Three values (L) to ponderer in the errors
  l      <- dim(v)[2]
  #Auxiliar matrices
  L_ones <- c(1,1,1)                           #A matrix with L ones.    It is used in the restrictions for w
  N_ones_hs <- matrix(1,n_hs,1)                #A matrix with n_hs ones. It is used in the restrictions for p and w
  J_ones <- matrix(1,J,1)                      #A matrix with J ones.    It is used in the restrictions for p

  #Priors TRHEE PRIORS ONE BY EACH OF L
  if(missing(q)|| is.null(q))  {q <- matrix(1/J, n_hs, J)}
  p_prior    <- q
  w_prior_1  <- matrix(1/l,n_hs,J,byrow=T)
  w_prior_2  <- w_prior_1
  w_prior_3  <- w_prior_1

  # Define the Shannon entropy function (I multiply it by -1 to maximize it)
  kl_dual <- function(lambda_v) {
    lambda        <- matrix(lambda_v, nrow=J ,ncol=, byrow=F)
    omega         <-rowSums(p_prior * t(exp(lambda %*% t(x_hs * w_hs))))
    #ROWSUM FOR SUMMATORY, T FOR DOMAIN TO BE WHAT IT HAS TO BE
    psi           <- w_prior_1 * t(exp(v[1]* (lambda %*% t(x_hs * w_hs)))) + w_prior_2 * t(exp(v[2]* (lambda %*% t(x_hs * w_hs)))) + w_prior_3 * t(exp(v[3]* (lambda %*% t(x_hs * w_hs))))
    #we do the multiplication in three steps because we have to multiply the support vector by the product of lambda and Xt(k,n_hs) by the weights of the HS,
    
    #Objective function
    sum(log(omega))  + sum(log(psi)) - sum(t(x_hp) %*%  (y_prev *  w_hp) * t(lambda)) } #Lambda multiplies element by element in the objetive function
  
  #Initial values
  lambda_ini      <- matrix(0,J,k)
  lambda_v <- as.vector(lambda_ini)

  
  if (missing(tol) || is.null(tol)) {
    tol <- 1e-10
  }
  if (missing(iter) || is.null(iter)) {
    iter <- 1000
  }
  
  llamar_nlp <- function(par, fn, lower = -Inf, upper = Inf,  ...) {
    
    control = list(iter.max = iter, rel.tol = tol)
    resultado <- stats::nlminb(start = par, objective = fn, lower = lower, upper = upper, control = control, ...)
    
    return(resultado)
  }
  
  res <- llamar_nlp(par = lambda_v, fn = kl_dual)


  lambda          <- matrix(res$par, nrow=J ,ncol=, byrow=F)

  #Final estimation of omega and psi
  omega         <- rowSums(p_prior * t(exp(lambda %*% t(x_hs * w_hs))))
  psi           <- w_prior_1 * t(exp(v[1]* (lambda %*% t(x_hs * w_hs)))) + w_prior_2 * t(exp(v[2]* (lambda %*% t(x_hs * w_hs)))) + w_prior_3 * t(exp(v[3]* (lambda %*% t(x_hs * w_hs))))#we do the multiplication in three steps because we have to multiply the support vector by the product of lambda and Xt(k,n_hs) by the weights of the HS, CHECK THAT IT IS EQUIVALENT
  
  #Estimation of p_dual and errors from lambdas
  p_dual   <-  (p_prior*t(exp(lambda %*% t(x_hs * w_hs)))) / omega
  

  u_dual_1 <- w_prior_1 * t(exp(v[1]* (lambda %*% t(x_hs * w_hs)))) /(psi)
  u_dual_2 <- w_prior_2 * t(exp(v[2]* (lambda %*% t(x_hs * w_hs)))) /(psi)
  u_dual_3 <- w_prior_3 * t(exp(v[3]* (lambda %*% t(x_hs * w_hs))))/(psi)
  error_dual <- u_dual_1 * v[1] + u_dual_2 * v[2] + u_dual_3 * v[3]

  predictions <-p_dual+ error_dual

  w <- cbind(u_dual_1,u_dual_2,u_dual_3)


  #Restrictions
  g1 <-  p_dual %*% J_ones - N_ones_hs             #First  set of restrictions defines the sum of probabilities for each household in hs must be one - n_hs restrictions
  g2 <-  u_dual_1 + u_dual_2 + u_dual_3 - 1

  #We multiply each w_j by the vector v. Then, each column of results is join together by colums in a matrix n_hs by J
  error_primal <- error_dual

  #Cross moments in the hp: the mean in each region for each k: X(k,j)
  cross_moments_hp<- (t(x_hp) %*%  (y_prev  *  w_hp))


  #Cross moments in the hs. "y" is weighted with "w_hs" by a multiplication element by element (NOT a matrix multiplication)
  cross_moments_hs <- (t(x_hs) %*% ((p_dual * w_hs) + ( error_primal  * w_hs)))

  g3 <- cross_moments_hp - cross_moments_hs  #Last set of restriction defines the cross moments in hp must be equal to hs
 
  cross_moments_hs <- noquote(apply(cross_moments_hs, c(1, 2), function(x) sprintf("%10.2f", x)))
  cross_moments_hp <- noquote(apply(cross_moments_hp, c(1, 2), function(x) sprintf("%10.2f", x)))
  
  #Set the colnames for the table
  table <- data.frame(w_hs,predictions, p_dual, error_dual)                                                                        #Joining all the relevant information in one table
  category_names <- levels(y_factor)
  assign_col_names <- function(df, y_factor) {
    w_names <-("weights")
    prediction_names <- paste0("predictions_", category_names)  # Usamos los nombres de las categorías
    p_names <- paste0("p_dual_", category_names)
    e_names <- paste0("e_dual_", category_names)
    all<- c("weights", prediction_names, p_names,e_names)
    colnames(df) <- all
    return(df)}
  
  
  colnames(cross_moments_hs) <- category_names
  colnames(cross_moments_hp) <- category_names
  
  table <- assign_col_names(table, y_factor)
  #the values of the optimization to show in the output
  values<- list(
    divergencekl =res$objective,
    iterations = res$iterations,
    message= res$message
  )
  k_names <-  c("(Intercept)", attr(terms(fn), "term.labels"))
  row.names(g3)<-k_names

  checkrestrictions<-list(g1=g1,g2=g2,g3=g3)

  
  checkrestrictions <- lapply(checkrestrictions, function(x) {
  
    if (is.matrix(x) || is.data.frame(x)) {
      apply(x, 2, function(y) round(as.numeric(y), 3))
    } else {
      round(as.numeric(x), 3)
    }
  })

  table2<-data.frame(lambda)


  
  
  # Asignar los nombres de las filas en la tabla
  colnames(table2) <- k_names
  rownames(table2)<-category_names
  
  generate_output <- function(estimations, values, tol, v, lambda, checkrestrictions, cross_moments_hp, cross_moments_hs, J, fn,divergencekl) {
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
  row.names(checkrestrictions$g3)<-k_names
  output <- generate_output(estimations = table,values =values,tol=tol,v=v, lambda=table2,checkrestrictions=checkrestrictions,cross_moments_hp=cross_moments_hp,cross_moments_hs=cross_moments_hs,J=J,fn=fn)

  class(output) <- "kl"
  return(output)
}

#' Generate a Plot
#'
#' @description This function generates a descriptive plot using the results obtained in ei_gce. It illustrates the mean and the confidence interval by disaggregated territorial unit.
#'
#' @param x The output produced by ei_gce
#' @param reg The data column containing the disaggregated territorial units
#' @param ... Additional arguments passed to the plotting function.
#' @return This function provides a graph representing the weighted mean and confidence interval of each disaggregated territorial unit
#' @import dplyr
#' @export
plot.kl <- function(x,reg,...){
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("This package requires the 'dplyr' package. Please install and load it.")
  }
  output<-x
  category_names <- rownames(output$lambda)
  table<-data.frame(reg,output$estimations)
  J=output$J
  regmeans <- table %>%
    group_by(reg) %>%
    summarise(across(
      everything(), 
      ~ weighted.mean(.x, w = weights, na.rm = TRUE)
    ))


  ic_lower <- regmeans$reg
  ic_upper <- regmeans$reg

  for (j in seq_along(category_names)) {
    col_name <- paste0("predictions_", category_names[j])
    ic_lower_col <- paste0("ic_lower_", category_names[j])
    ic_upper_col <- paste0("ic_upper_", category_names[j])
    
    regmeans <- regmeans %>%
      dplyr::mutate(
        !!ic_lower_col := !!dplyr::sym(col_name) - 1.96 * sd(!!dplyr::sym(col_name)),
        !!ic_upper_col := !!dplyr::sym(col_name) + 1.96 * sd(!!dplyr::sym(col_name))
      )
  }


  for (j in seq_along(category_names)){
    prediction_col <- paste0("predictions_", category_names[j])
    ic_lower_col <- paste0("ic_lower_", category_names[j])
    ic_upper_col <- paste0("ic_upper_", category_names[j])
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
  }
  }

#' Summary
#'
#' @description This function provides a summary of the output obtained with the function ei_gce.
#'
#'
#' @param object The output obtained from ei_gce
#' @param ... Additional arguments passed to the summary function.
#' @return This summary function returns the Kullback-Leibler divergence value and the last iteration in the optimization process.
#'         A dataframe with the means of the estimations for each characteristic j with the predictions the probabilities and the error estimated.
#'         A dataframe with the lambda estimated for each k.
#'  \itemize{
#'   \item \code{Iterations}:Indicates the times the objective function and the gradient has been evaluated during the optimization process
#'   \item \code{divergencekl value}:The Kullback-Leibler divergence value resulting from the optimization.
#'   \item \code{mean_estimations}:The weighted mean of predictions, p_dual, and the error for each category j of the variable y
#'   \item \code{lambda}:The estimated lambda values.
#'  }
#' @import dplyr
#' @export
summary.kl <- function(object,...){
    output<-object
    fn=output$fn
    j=output$J
    category_names <- rownames(output$lambda)

    cat ("Iterations")
    print(output$values$iterations[1])
    cat ("Kullback-Leibler divergence value")
    print(output$values$divergencekl[1])

    
    
    mean_estimations <- output$estimations %>%
      summarise(across(everything(), ~ weighted.mean(.x, w = output$estimations$weights, na.rm = TRUE)))
    print("mean_estimations")
    print(mean_estimations)
  print ("lambda")
  print(output$lambda)}








