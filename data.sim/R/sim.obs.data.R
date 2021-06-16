#' @title Random linear model observed data generator
#' @description Simulate random observed data based on input data and parameters.
#' @author Alex Bajcz
#'  
#' @param Xs A data frame containing 1 or more columns of X (independent variable) data. The data can be any combination of factors and quantitative data. However, the function will assume that the R-default factor levels, so custom factor levels will be ignored/lost. 
#' @param param.mat Either a vector or a matrix with a single column containing the designed parameter values that should make up the parameter matrix within the linear predictor. If supplied as a vector, it will be internally coerced into a matrix. This should not be a data frame. The column need not be named. Logically, the number of values specified should match the number of parameters in your eventual design matrix (taking into account random effects and factor variable recoding) exactly; if this is not the case, the function will return an error to this effect. So, care needs to be taken to carefully consider the number of values required. In practice, if random slopes/intercepts are assumed, some of the initial values provided to this argument will be overridden with randomly selected values from a distribution with hyperoparameter values you specify. 
#' @param family A single character element corresponding to the family of linear model to be simulated. Can take any of the following values: "gaussian", "gamma", "inverse.gaussian", "poisson", "bernoulli", "binomial", "negative.binomial". Defaults to "gaussian".
#' @param link A single character element corresponding to the link function of the linear model to be simulated. Can take any of the following values: "identity", "log", "inverse.squared", "logit", "sqrt", "inverse", "cloglog". Keep in mind that some family-link combinations would be non-standard and even potentially non-sensical, although the model will attempt to run any combination specified. Defaults to "identity". The only common linking function not currently implemented is the probit. 
#' @param sigma A single numerical value corresponding to the standard deviation of the residual error distribution in a gaussian family model. Must be strictly >0. Defaults to NULL.
#' @param lambda A single numerical value corresponding to the shape parameter of the residual error distribution in an inverse-gaussian family model. Must be strictly >0. Defaults to NULL.
#' @param N.trials A vector of length = nrows(Xs) containing the number of trials involved in generating Binomial-family data that are not Bernoullian. Since these values essentially represent denominators in a positive fraction, they must all be strictly >0. Defaults to NULL.
#' @param scale A single numerical value corresponding to a scale parameter of the residual error distribution in an gamma family model. Must be strictly >0. Defaults to NULL.
#' @param size A single numerical value corresponding to a size parameter of the residual error distribution in a negative binomial family model. Must be strictly >0. Defaults to NULL.
#' @param random.slopes A single logical value (TRUE or FALSE) indicating whether random slopes should be assumed for some covariate based on a specified grouping variable. Defaults to FALSE. If set to TRUE, values for group.slopes, mu.slopes, and sigma.slopes must also be specified. 
#' @param random.intercepts A single logical value (TRUE or FALSE) indicating whether random intercepts should be assumed for a specified grouping variable. Defaults to FALSE. If set to TRUE, values for group.intercepts, mu.intercepts, and sigma.intercepts must also be specified. 
#' @param slope.int.covary A single logical value (TRUE or FALSE) indicating whether random intercepts and random slopes should be assumed to be covariable. Defaults to FALSE. If set to TRUE, random.slopes and random.intercepts must also be TRUE, and additionally a value for var.covar.mat must be supplied. 
#' @param var.covar.mat A strictly square 2x2 matrix representing the desired variance-covariance matrix for drawing random correlated slope and intercept parameter values from a bivariate normal distribution (see Kery 2010, page 161 for details). The diagonal of this matrix should be equal to 1, and every value should be strictly between 0 and 1. In it's current form, the function can only handle random effects based on a single grouping variable. 
#' @param rand.int.term A single character value. Defaults to NULL. When random.slopes = TRUE, a value for this parameter must be specified. The form required is very specific: It must be "X~:X^" where ~ is the column number in Xs of the grouping variable of the random effects and where ^ is the column number in Xs of a continuous covariate whose slopes will vary by levels of the grouping variable (it is also ok to reverse these two numbers). For example, if the grouping variable for which there will be random slopes is in column 7 of Xs and if the covariate for which these slopes will vary by group is in column 2, "X2:X7" or "X7:X2" are both acceptable forms. 
#' @param rand.nongroup.var A single character value. Defaults to NULL. When random.slopes = TRUE, a value for this parameter must be specified. The form required is very specific: It must be "X^" where ^ is the column number in Xs of a continuous covariate whose slopes will vary by levels of the grouping variable (see under rand.int.term for more details). 
#' @param group.intercepts A vector of length up to but not exceeding the length of param.mat. Defaults to NULL. When random.intercepts = TRUE, a value for group.intercepts must be provided. It should correspond to the position values in param.mat that hold the parameter estimates for the grouping variable that will have random intercepts. Keep in mind that the function will remove the global intercept term if random intercepts are assumed, which may affect which values in param.mat need to correspond to the grouping variable intercepts. See the examples for details. In the event that slope.int.covary is TRUE, values for group.intercepts and group.slopes need to be provide. These will then be concatenated internally (in that order), and the resulting vector cannot exceed param.mat in length. 
#' @param mu.intercepts A single numerical value corresponding to the mean of the hyperdistribution from which new random intercept parameter values will be chosen, one for each group in the grouping variable specified. Defaults to NULL but must be provided if random.intercepts is TRUE.
#' @param sigma.intercepts A single numerical value corresponding to the standard deviation of the hyperdistribution from which new random intercept parameter values will be chosen, one for each group in the grouping variable specified. Defaults to NULL but must be provided if random.intercepts is TRUE.
#' @param group.slopes A vector of length up to but not exceeding the length of param.mat. Defaults to NULL. When random.slopes = TRUE, a value for group.slopes must be provided. It should correspond to the position values in param.mat that hold the parameter estimates for the continuous covariate variable that will have random slopes by group. See the examples for details. In the event that slope.int.covary is TRUE, values for group.intercepts and group.slopes need to be provided. These will then be concatenated internally (in that order), and the resulting vector cannot exceed param.mat in length. 
#' @param mu.slopes A single numerical value corresponding to the mean of the hyperdistribution from which new random slope parameter values between Y and 1 continuous covariated specified will be chosen, one for each group in the grouping variable specified. Defaults to NULL but must be provided if random.slopes is TRUE.
#' @param sigma.slopes A single numerical value corresponding to the standard deviation of the hyperdistribution from which new random slope parameter values between Y and 1 continuous covariated specified will be chosen, one for each group in the grouping variable specified. Defaults to NULL but must be provided if random.slopes is TRUE.
#' @return A list containing two vectors, /code{exp.ys} and /code[obs.ys], which correspond to the expected and observed Y values simulated by the function, given the inputs. 
#' @export
#' @details The sim.obs.data function generates observed and expected
#'   "generalized linear regression" response data utilizing independent
#'   variable data, parameter values, model structure, and assumptions you
#'   specify. It can handle most standard model families with their typical
#'   linking functions, any number of covariates (including factors with
#'   multiple levels that will need to be dummy coded), and random slopes and/or
#'   intercepts based on a single grouping variable. However, the mechanics of
#'   complex models will require inputs to be specified in a very specific
#'   manner.
#'
#'   The function will assume there is a global intercept term in the model
#'   unless random effects have been specified (in which case the model will
#'   switch instead of mean encoding). The global intercept term could be
#'   countered by setting the first value in param.mat to 0, but keep in mind
#'   that that would effectively coerce the data to have a y intercept of 0,
#'   which rarely makes sense.
#'
#'   If there are polynomial, interaction, or other complex terms in your model,
#'   these should be precalculated and comprised within the input to Xs. In
#'   other words, if you are trying to fit ~X1 + X2 + X1:X2, there should be a
#'   column in your input to Xs that is already the product of X1 and X2, and
#'   there should also be a corresponding beta term in the param.mat input
#'   vector for this term.
#' @examples ##This example shows a workflow for generating an entire simulated 
#' #data set (multiple X variables plus 1 Y variable) using sim.pred.data() plus 
#' #sim.obs.data that shows the typical syntax for a call to sim.obs.data plus 
#' #many of its features. It's a somewhat trivial case because it calls all 
#' #seven available families (with their standard linking functions) once using 
#' #the same input linear coefficient values even though the likely coefficient 
#' #values will tend to be very different for, as an example, a gamma model 
#' #versus a binomial one.
#'
#' #Generates a predictor data set with two Xs, one numeric and one factorial.
#' Xstest = sim.data::sim.pred.data(col.num = 2, col.length = 25, 
#'          FUNs = c("rnorm", "sample"), factor.YN = c(FALSE, TRUE),
#'          x.sample = 0:2, replace = TRUE, sample.prob = c(0.4, 0.2, 0.4), 
#'          mean = 2, sd = 0.5)
#' print(Xstest) #shows what's we've made.
#'
#' # we'll simulate one Y variable for each model family implemented using its 
#' #standard linking function.
#' family.test = c("gaussian", "gamma", "inverse.gaussian", "poisson", 
#'               "bernoulli", "binomial", "negative.binomial")
#' link.test = link = c("identity", "inverse", "inverse.squared", "log", 
#'             "logit", "logit", "log")
#'
#' #We'll randomly draw parameter values for the specialized families from 
#' #plausible distributions.
#' sigma.test = runif(1, min=0.5, max=5)
#' lambda.test = runif(1, min=0.5, max=5)
#' N.trials.test = round(rnorm(length(Xstest[,1]), 40, 10))
#' scale.test = runif(1, min=0.5, max=5)
#' size.test = runif(1, min=0.5, max=5)
#' param.mat.test1 = c(0.1, 0.2, 2.5, -0.4) #The function will coerce this to a
#' #matrix with 1 column for you.
#'
#' # we'll create empty matrices to store our output Y data.
#' record.obs = matrix(data = NA, nrow=length(Xstest[,1]), ncol=length(family.test))
#' record.exp = matrix(data = NA, nrow=length(Xstest[,1]), ncol=length(family.test))
#'
#' #For each family-link combo, we'll call sim.obs.data to generate a set of 
#' #appropriate Y data. sim.obs.data will ignore inputs that are not essential 
#' #for the family being called.
#' for(trial in 1:length(family.test)) {
#'   test1 = sim.obs.data(Xs = Xstest, param.mat = param.mat.test1, 
#'           family = family.test[trial], link = link.test[trial], 
#'           sigma = sigma.test, lambda = lambda.test, N.trials = N.trials.test,
#'                       scale = scale.test, size = size.test, 
#'                       random.intercepts = TRUE,group.intercepts = 2:4, 
#'                       mu.intercepts = 1, sigma.intercepts = 0.5)
#'  record.obs[,trial] = test1$obs.ys #record the observed data generated.
#'  record.exp[,trial] = test1$exp.ys #record the "expected" data generated.
#' }
#'
#' print(record.obs) ; print(record.exp) #Show what we've made. 



sim.obs.data = function (Xs, param.mat, family = "gaussian", link = "identity", 
                         sigma = NULL, lambda = NULL, N.trials = NULL, scale = NULL, size = NULL,  
                         random.slopes = FALSE, random.intercepts = FALSE, slope.int.covary = FALSE, 
                         var.covar.mat = NULL, rand.int.term = NULL, rand.nongroup.var = NULL,
                         group.intercepts = NULL, mu.intercepts = NULL, sigma.intercepts = NULL, 
                         group.slopes = NULL, mu.slopes = NULL, sigma.slopes = NULL) {

  #One common error is likely to be providing a vector to param.mat even though it needs to be a matrix for matrix multiplication, so we will coerce it here to a matrix.
  param.mat = matrix(param.mat, ncol=1)
  
  #This section checks for a great many predictable errors related to function inputs and returns what I would hope is a helpful error message. For the most part, this should prevent users from ending up with errant output based on a misspecified function call. Some of these may be redundant and could be simplified. 
  if (dim(Xs)[2] < 1) { stop("Make sure to specify at least one X column to Xs!")}
  if (family == "gaussian" & is.null(sigma)) { stop("Family specified as gaussian but no sigma value provided!")}
  if (length(sigma) > 1) { stop("More than 1 sigma value was provided!")}
  if (sigma <= 0 ) { stop("A negative sigma was provided! Sigma must be >0.")}
  if (family == "inverse.gaussian" & is.null(lambda)) { stop("Family specified as inverse gaussian but no lambda value provided!")}
  if (length(lambda) > 1) { stop("More than 1 lambda value was provided!")}
  if (lambda <= 0 ) { stop("A negative lambda was provided! Lambda must be >0.")}
  if (family == "binomial" & is.null(N.trials)) { stop("Family specified as binomial but no N.trials vector provided!")}
  if (length(N.trials) != nrow(Xs)) { stop("The length of N.trials does not equal the number of rows in Xs.")}
  if (any(N.trials <= 0 )) { stop("At least one value in N.trials <= 0, which isn't logical!")}
  if (family == "gamma" & is.null(scale)) { stop("Family specified as gamma but no scale value provided!")}
  if (length(scale) > 1) { stop("More than 1 scale value was provided!")}
  if (scale <= 0 ) { stop("A negative scale was provided! Scale must be >0.")}
  if (family == "negative.binomial" & is.null(size)) { stop("Family specified as negative binomial but no size value provided!")}
  if (length(size) > 1) { stop("More than 1 size value was provided!")}
  if (size <= 0 ) { stop("A negative size was provided! Size must be >0.")}
  if (random.slopes !=TRUE & random.slopes != FALSE) { stop("Invalid value provided for random.slopes!")}
  if (random.slopes == TRUE & (is.null(group.slopes) | is.null(mu.slopes)| is.null(sigma.slopes) | 
      is.null(rand.nongroup.var) | is.null(rand.int.term))) { stop("Ensure that group.slopes, mu.slopes, rand.int.term, rand.nongroup.var, and sigma.slopes are specified if random.slopes is TRUE!")}
  if (random.intercepts !=TRUE & random.intercepts != FALSE) { stop("Invalid value provided for random.intercepts!")}
  if (random.intercepts == TRUE & is.null(group.intercepts) & (is.null(mu.intercepts)| is.null(sigma.intercepts))) { stop("Ensure that group.intercepts, mu.intercepts, and sigma.intercepts are specified if random.intercepts is TRUE!")}
  if (slope.int.covary !=TRUE & slope.int.covary!= FALSE) { stop("Invalid value provided for slope.int.covary!")}
  if (slope.int.covary == TRUE & isFALSE(random.slopes) & (isFALSE(random.intercepts) | is.null(var.covar.mat))) { stop("Ensure that random.intercepts and random.slopes are TRUE and that a value for var.covar.mat is provided if slope.int.covary is TRUE!")}
  if (!is.character(rand.int.term) & !is.null(rand.int.term)) { stop("rand.int.term must be a character!")}
  if (!is.null(var.covar.mat) & length(var.covar.mat) != 4) { stop("var.covar.mat should be a 2x2 matrix!")}
  if (slope.int.covary == TRUE & is.null(var.covar.mat)) { stop("A value for var.covar.mat should be specified if slope.int.covary is true!")}
  if (!is.null(var.covar.mat) && (diag(var.covar.mat) != c(1,1) | any(var.covar.mat < 0) | any(var.covar.mat > 1))) { stop ("All values in var.covar.mat should be strictly between 0 and 1 and the values on the diagonal should be equal to 1!")}
  if (!is.null(rand.nongroup.var) && !is.character(rand.nongroup.var)) { stop("Invalid value provided for rand.nongroup.var! It should be a character and of length 1!") }
  if (!is.null(rand.nongroup.var) & (length(rand.nongroup.var) != 1)) { stop("Invalid value provided for rand.nongroup.var! It should be a character and of length 1!") }
  if (length(group.intercepts) > length(param.mat)) { stop("The length of group.intercepts should not exceed that of param.mat!")}
  if (random.intercepts == TRUE & is.null(group.intercepts)) { stop("A value for group.intercepts must be provided if random.intercepts is TRUE!")}
  if (slope.int.covary == TRUE & (is.null(group.intercepts) | is.null(group.slopes))) { stop("If slope.int.covary is true, values for group.intercepts and group.slopes must both be provided!")}
  if (slope.int.covary == TRUE & length(c(group.intercepts, group.slopes)) > length(param.mat)) { stop("Too many combined parameters provided to group.slopes and group.intercepts! The length should not exceed that of param.mat.")}
  if (length(group.slopes) > length(param.mat)) { stop("The length of group.slopes should not exceed that of param.mat!")}

  #The first step is to "decompose" the columns provided as inputs to Xs so that we can assemble them into a predictable formula. However, we need to make sure we preserve the factor nature of any factor variables.
    for (col in 1:dim(Xs)[2]) { #For each column in Xs...
    if (!is.factor(Xs[,col])) { #If this column isn't a factor...
    assign(paste0("X", col), Xs[,col]) #Make a new object called X*, where * is equal to the current column number.
    } else {
    assign(paste0("X", col), as.factor(Xs[,col])) #If the current column was a factor, make sure the new object is one also. 
    }
  }
  
  #Step 2: Assemble the formula that model.matrix will need as input to generate the design matrix. This step will need to vary a bit depending on whether we're simulating fixed effects, random effects, or both. 
  
  formula1 = paste0("X1") #We assume there is at least one X variable provided and is now called X1. 
  if (dim(Xs)[2] > 1) { #Provided there are additional X variables...
    for (term in 2:dim(Xs)[2]) { #For each additional term...
    formula1 = paste0(formula1, "+", "X", term) #Add this term to the formula character string.
    }
  }
    if(random.intercepts == TRUE) { #If the user specifies random intercepts, we also need to subtract the intercept term from the model so each group has its own unique but dependent intercept.
      formula1 = paste0(formula1, "-1")
    }
    if(random.slopes == TRUE) { #If the user specifies random slopes, we also will need to add in a group:covariate interaction term, as specified to the rand.int.term argument, and subtract out the main effect for the continuous variable, as specified to the rand.nongroup.var argument. 
      formula1 = paste0(formula1, "+", rand.int.term, "-", rand.nongroup.var)
    }
  formula1 = paste0("~", formula1) #Lastly, add our formula operator the front. 

  #With the proper input assembled, we can construct our design matrix.
  design.mat = model.matrix(as.formula(formula1))

  #Since it is very easy to miscount the number of parameters that will result from a model with many factor variables and/or random effects, we now check to make sure the model is specified correctly and abort if not. 
  if(dim(design.mat)[2] != dim(param.mat)[1]) { stop("Model misspecified! The number of columns in the design matrix does not correspond to the number of values in param.mat.")}
  
  #Step 3 is to incorporate any random effects specified by overwriting values in param.mat that correspond to group-level parameters with new draws from the specified hyperdistribution. 
  
  if (random.intercepts == TRUE & slope.int.covary == FALSE) {
    param.mat[group.intercepts,1] = rnorm(length(group.intercepts), mu.intercepts, sigma.intercepts) #A potential problem with this code is that it assumes group.intercepts has been specified "correctly" in both content and length...
  }
  if (random.slopes == TRUE & slope.int.covary == FALSE) {
    param.mat[group.slopes, 1] = rnorm(length(group.slopes), mu.slopes, sigma.slopes) #Same problem as above.
  }
  if (slope.int.covary == TRUE) {
    param.mat[c(group.slopes, group.intercepts), 1] = MASS::mvrnorm(length(c(group.slopes, group.intercepts)), c(mu.intercepts, mu.slopes), var.covar.mat)
  }
  
  #Step 4 is to create our linear predictor by matrix multiplying our design matrix by our parameter matrix.
  linear.predict = design.mat %*% param.mat
  
  #Step 5 is to generate a set of expected Y values based on the linear predictor we just calculated and linking function provided as input. I am *pretty sure* I have programmed these correctly but this should be checked. 
  if(link == "identity") {
    exp.ys = linear.predict
  }
  if (link == "log") {
    exp.ys = exp(linear.predict)
  }
  if (link == "inverse.squared") {
    exp.ys = (1/sqrt(linear.predict))
  }
  if (link == "logit") { #Not really expected Ys but rather estimates of the probability of successes in each trial.
    exp.ys = (exp(linear.predict)/(1 + (exp(linear.predict))))
  }
  if (link == "sqrt") {
    exp.ys = (linear.predict)^2
  }
  if (link == "inverse") {
    exp.ys = (linear.predict)^-1
  }
  if (link == "cloglog") { #Not really expected Ys but rather estimates of the probability of successes in each trial.
    exp.ys = -exp(-exp((linear.predict))-1)
  }
  
  if(is.null(exp.ys)) { stop("Incorrect link specified!") } #This trips most likely if someone specifies a link function that isn't coded.

  
  #Step 6 is to use our expected Y values and our assumed residual random error distribution to generate random noise, based on the distribution family specified as inputs. 
  
  N = length(X1) #The number of expected Y values to draw. More concise than referring to this value directly.
  
  if(family == "gaussian") {
    obs.ys = rnorm(N, exp.ys, sigma) #sigma is required if family is gaussian. 
  }
  if(family == "gamma") {          
    obs.ys = rgamma(N, shape = exp.ys, scale = exp.ys*scale) #scale is required if family is gamma. I think this is a fair implementation here, but the gamma distribution is tough so this should be confirmed.  
  }
  if(family == "inverse.gaussian") {
    obs.ys = SuppDists::rinvGauss(N, exp.ys, lambda = lambda) #lambda is required if family is inverse.gaussian. 
  }
  if(family == "poisson") {
    obs.ys = rpois(N, exp.ys)
  }
  if(family == "bernoulli" | family == "binomial") { #If the family is binomial or bernoulli, to draw new values, we have to provide the number of trials involved in each event. This is 1 for every event in the case of bernoulli (every even is a single "coin-flip").
      if (family == "bernoulli" ) { 
        usethisN = 1 } else {
          usethisN = N.trials }
          obs.ys = rbinom(N, size=usethisN, prob=exp.ys) #N.trials must be provided if family is binomial.
  }
  if(family == "negative.binomial") {
    obs.ys = rnbinom(N, mu = exp.ys, size = size) #I think this is a fair implementation here, but the negative binomial distribution is tough so this should be confirmed.  
  }
  
  if(is.null(obs.ys)) { stop("Incorrect family specified!") } #The most likely trip here is if someone specifies a family that isn't coded.
  
#Step 7 is to output the expected and observed y values generated. 
  return(list(obs.ys = obs.ys, exp.ys = exp.ys))
  
}


