#'@name ei_gce
#'@title Ecologic Inference applying entropy
#'@description The function ei_gce defines the Kullback-Leibler function which
#'minimises the distance between the distribution of probabilities P and the
#'distribution Q. The distribution Q is based on prior information that we have
#'of our variable of interest previous to the analysis. The function will set
#'the optimization parameters and, using the "nlminb"
#'function, an optimal solution is obtained.
#'The function defines the independent variables in the two databases
#'needed, which we call dataA with "n_A" observations and dataB with "n_B"
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
#'@param fn Is the formula that represents the dependent variable in the optimization.
#' In the context of this function, 'fn' is used to define the dependent variable
#' to be optimized by the Kullback-Leibler divergence function.
#' Note: If the dependent variable is categorical the sorting criterion for the columns, and therefore for J, is alphabetical order.
#'@param dataA The data where the variable of interest y is available and also the independent variables.
#'  Note: The variables and weights used as independent variables must have the same name in 'dataA' and in 'dataB'
#'@param dataB The data which contains information on the independent variables at a disaggregated level.
#'  Note: The variables and weights used as independent variables must have the same name in 'dataA' and in 'dataB'. The variables in both databases need to match up in content.
#'@param weights A character string specifying the column name to be used as weights in both 'dataA'  and 'dataB' datasets.
#' If the argument `weights` is provided and present in both datasets,  the weights in each dataset will be normalized by the sum of the weights within that dataset.
#' If `weights` is NULL or the specified column does not exist in both datasets, equal weights are applied across all observations.
#'@param tol The tolerance to be applied in the optimization function. If the tolerance is not specified, the default tolerance has been set in 1e-10
#'@param q The prior distribution Q
#'@param v The support vector
#'@param iter The maximum number of iterations allowed for the optimization algorithm to run 
#' Increasing the number of iterations may improve the likelihood of finding an optimal solution, 
#' but can also increases computation time.If the maximum number of iterations is not specified, it will default to 1000
#' @return The function will provide you a dataframe called estimations with the next information:
#'    \itemize{
#'   \item \strong{weights} The weights used in the optimization process.
#'   \item \strong{predictions}  The prediction for each individual is calculated as the sum of the probability plus the error.
#'   \item \strong{probabilities}  Probabilities for each individual to each possibility \code{j} of the variable of interest \code{y}.
#'   \item \strong{errors}  Errors calculated to the \code{j} possibilities of \code{y}.
#'   The function provides information about the optimization process as:
#'   \item \strong{divergencekl} The Kullback-Leibler divergence value resulting from the optimization.
#'   \item \strong{iterations} Indicates the times the objective function and the gradient has been evaluated during the optimization process,if any.
#'   \item \strong{message} Indicates the message if it has been generated in the process of optimization.
#'   \item \strong{q} Indicates prior implemented in the optimization.
#'   \item \strong{tol} Indicates the tolerance of the optimization process.
#'   \item \strong{v} Indicates the support vector used in the function.
#'   The function provides a dataframe containing the information about lambda:
#'   \item  \strong{lambda} The estimated lambda values.
#'   It is provided an object with the restrictions checked which should be zero.
#'   \item  \strong{check restrictions} Being  \code{g1}  the restriction related to the unit probability constraint, \code{g2}  to the error unit sum constraint, and \code{g3}  to the consistency restriction that implies that the difference between the cross moment in both datasets must be zero.
#'   The restriction g3 can be checked thoroughly with the objects by separate.
#'   \item  \strong{cross moments A} Cross moments in \code{dataA}.
#'   \item  \strong{cross moments B} Cross moments in \code{dataB}.
#'   }
#' @references
#' Fernandez-Vazquez, E., Diaz-Dapena, A., Rubiera-Morollon, F., Viñuela, A., (2020) Spatial Disaggregation of Social Indicators: An Info-Metrics Approach. Social Indicators Research, 152(2), 809–821. https://doi.org/10.1007/s11205-020-02455-z.
#' @examples
#' #In this example we use the data of this package
#' dataA <- financial()
#' dataB <- social()
#' # Setting up our function for the dependent variable.
#' fn               <- dataA$poor_liq ~ Dcollege+Totalincome+Dunemp
#' #In this case we know that the mean probability of being poor is 0.35.With this function
#' #we can add the information as information a priori. This information a priori correspond to the
#' #Q distribution and in this function is called q for the sake of simplicity:
#' q<- c(0.35,0.65)
#' v<- matrix(c(0.2,0,-0.2))
#' #Applying the function ei_gce to our databases. In this case dataA is the
#' #dataA where we have our variable of interest
#' #dataB is the data where we have the information for the disaggregation.
#' #w can be included if we have weights in both surveys
#' result  <- ei_gce(fn,dataA,dataB,q=q,weights="w",v=v)
#' @export
ei_gce<- function(fn,dataA,dataB,q,weights=NULL,v,tol,iter){

  x_A             <- model.matrix(fn,dataA)                         #Matrix of independent variables in A - dimensions n_B (observations in A) and k (parameters)
  
  fn_right         <- fn[-2]                                          #B does not have dependent variable. It is needed a formula without the left hand side
  x_B             <- model.matrix(fn_right,dataB)                   #Matrix of independent variables in B - dimensions n_A (observations in B) and k (parameters)
  
  #loading and rescale of weights (by the total) in each survey
  if (!is.null(weights) && weights %in% colnames(dataA) && weights %in% colnames(dataB)) {
    # Use the specified column for weights, normalizing by the sum of weights
    w_A <- dataA[[weights]] / sum(dataA[[weights]], na.rm = TRUE)
    w_B <- dataB[[weights]] / sum(dataB[[weights]], na.rm = TRUE)
  } else {
    # If no weights or column is missing, apply equal weights
    w_A <- rep(1 / nrow(dataA), nrow(dataA))
    w_B <- rep(1 / nrow(dataB), nrow(dataB))
  }
  if (any(is.na(x_B))) {
    stop("NA (missing) values have been found in the dataset")
  }
  if (any(is.na(x_A))) {
    stop("NA (missing) values have been found in the dataset")
  }
  
  #Definition of "y" as a matrix with n_A observations by j categories.
  form <- model.frame(fn,dataA)
  y          <- model.response(form) 
  if (length(unique(y)) == 2 && all(sort(unique(y)) == c(0, 1))) {
    y_factor <- factor(y, levels = c(1, 0))  
  } else {
    y_factor <- factor(y)
  }
  y_prev <- model.matrix(~ y_factor - 1)
  

  # Set the parameters
  n_A <-dim(x_A)[1]  ; n_A                 #observations in A
  n_B <- dim(x_B)[1]  ; n_B                 #Observations in B
  k    <- dim(x_B)[2]  ; k                    #columns in the matrix of independent variables (including the intercept)
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
  N_ones_B <- matrix(1,n_B,1)                #A matrix with n_B ones. It is used in the restrictions for p and w
  J_ones <- matrix(1,J,1)                      #A matrix with J ones.    It is used in the restrictions for p

  #Priors TRHEE PRIORS ONE BY EACH OF L
  if (missing(q) || is.null(q)) {
    q <- matrix(1/J, n_B, J)  # Si q falta o es NULL, se inicializa como una matriz con valores 1/J
  } else if (is.atomic(q) && length(q) == J) {  
    # Si q es un vector numérico con longitud J, se expande en una matriz de n_B filas
    q <- matrix(rep(q, each = n_B), nrow = n_B, ncol = J, byrow = FALSE)
  }
  p_prior    <- q
  w_prior_1  <- matrix(1/l,n_B,J,byrow=T)
  w_prior_2  <- w_prior_1
  w_prior_3  <- w_prior_1

  # Define the Shannon entropy function (I multiply it by -1 to maximize it)
  kl_dual <- function(lambda_v) {
    lambda        <- matrix(lambda_v, nrow=J ,ncol=, byrow=F)
    omega         <-rowSums(p_prior * t(exp(lambda %*% t(x_B * w_B))))
    #ROWSUM FOR SUMMATORY, T FOR DOMAIN TO BE WHAT IT HAS TO BE
    psi           <- w_prior_1 * t(exp(v[1]* (lambda %*% t(x_B * w_B)))) + w_prior_2 * t(exp(v[2]* (lambda %*% t(x_B * w_B)))) + w_prior_3 * t(exp(v[3]* (lambda %*% t(x_B * w_B))))
    #we do the multiplication in three steps because we have to multiply the support vector by the product of lambda and Xt(k,n_B) by the weights of the B,
    
    #Objective function
    sum(log(omega))  + sum(log(psi)) - sum(t(x_A) %*%  (y_prev *  w_A) * t(lambda)) } #Lambda multiplies element by element in the objetive function
  
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
     
    # Asegurar que fn no devuelve NA/NaN
    fn_safe <- function(p) {
        val <- fn(p, ...)
        
        if (is.na(val) || is.nan(val) || is.infinite(val)) {
          return(Inf)  
        }
        
        return(val)
    }
    control = list(iter.max = iter, rel.tol = tol)
    resultado <- stats::nlminb(start = par, objective = fn_safe, lower = lower, upper = upper, control = control, ...)
    
    return(resultado)
  }
  
  res <- llamar_nlp(par = lambda_v, fn = kl_dual)


  lambda          <- matrix(res$par, nrow=J ,ncol=, byrow=F)

  #Final estimation of omega and psi
  omega         <- rowSums(p_prior * t(exp(lambda %*% t(x_B * w_B))))
  psi           <- w_prior_1 * t(exp(v[1]* (lambda %*% t(x_B * w_B)))) + w_prior_2 * t(exp(v[2]* (lambda %*% t(x_B * w_B)))) + w_prior_3 * t(exp(v[3]* (lambda %*% t(x_B * w_B))))#we do the multiplication in three steps because we have to multiply the support vector by the product of lambda and Xt(k,n_B) by the weights of the B, CHECK THAT IT IS EQUIVALENT
  
  #Estimation of probability and errors from lambdas
  probability   <-  (p_prior*t(exp(lambda %*% t(x_B * w_B)))) / omega
  

  u_dual_1 <- w_prior_1 * t(exp(v[1]* (lambda %*% t(x_B * w_B)))) /(psi)
  u_dual_2 <- w_prior_2 * t(exp(v[2]* (lambda %*% t(x_B * w_B)))) /(psi)
  u_dual_3 <- w_prior_3 * t(exp(v[3]* (lambda %*% t(x_B * w_B))))/(psi)
  error_dual <- u_dual_1 * v[1] + u_dual_2 * v[2] + u_dual_3 * v[3]

  predictions <-probability+ error_dual

  w <- cbind(u_dual_1,u_dual_2,u_dual_3)


  #Restrictions
  g1 <-  probability %*% J_ones - N_ones_B             #First  set of restrictions defines the sum of probabilities for each household in B must be one - n_B restrictions
  g2 <-  u_dual_1 + u_dual_2 + u_dual_3 - 1

  #We multiply each w_j by the vector v. Then, each column of results is join together by colums in a matrix n_B by J
  error_primal <- error_dual

  #Cross moments in the A: the mean in each region for each k: X(k,j)
  cross_moments_A<- (t(x_A) %*%  (y_prev  *  w_A))


  #Cross moments in the B. "y" is weighted with "w_B" by a multiplication element by element (NOT a matrix multiplication)
  cross_moments_B <- (t(x_B) %*% ((probability * w_B) + ( error_primal  * w_B)))

  g3 <- cross_moments_A - cross_moments_B  #Last set of restriction defines the cross moments in A must be equal to B
 
  cross_moments_B <- noquote(apply(cross_moments_B, c(1, 2), function(x) sprintf("%10.2f", x)))
  cross_moments_A <- noquote(apply(cross_moments_A, c(1, 2), function(x) sprintf("%10.2f", x)))
  
  #Set the colnames for the table
  table <- data.frame(w_B,predictions, probability, error_dual)                                                                        #Joining all the relevant information in one table
  category_names <- levels(y_factor)
  assign_col_names <- function(df, y_factor) {
    w_names <-("weights")
    prediction_names <- paste0("predictions_", category_names)  # Usamos los nombres de las categorías
    p_names <- paste0("probabilities_", category_names)
    e_names <- paste0("errors_", category_names)
    all<- c("weights", prediction_names, p_names,e_names)
    colnames(df) <- all
    return(df)}
  
  
  colnames(cross_moments_B) <- category_names
  colnames(cross_moments_A) <- category_names
  
  table <- assign_col_names(table, y_factor)
  #the values of the optimization to show in the output
  values<- list(
    divergencekl =res$objective,
    iterations = res$iterations,
    message= res$message,
    q=p_prior
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
  
  generate_output <- function(estimations, values, tol, v, lambda, checkrestrictions, cross_moments_A, cross_moments_B, J, fn,divergencekl) {
    output <- structure(list(
      estimations = estimations,
      values = values,
      tol = tol,
      v = v,
      lambda = lambda,
      checkrestrictions = checkrestrictions,
      cross_moments_A = cross_moments_A,
      cross_moments_B = cross_moments_B,
      J = J,
      fn = fn
    ), class = "shannon")
    return(output)
  }
  row.names(checkrestrictions$g3)<-k_names
  output <- generate_output(estimations = table,values =values,tol=tol,v=v, lambda=table2,checkrestrictions=checkrestrictions,cross_moments_A=cross_moments_A,cross_moments_B=cross_moments_B,J=J,fn=fn)

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
#'   \item \code{mean_estimations}:The weighted mean of predictions, probabilities, and the errors for each category j of the variable y
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








