#' @title Random linear model predictor data generator
#' @description Generates random predictor (X) data based on input parameters and specifications.
#' @author Alex Bajcz
#' 
#' @param col.num A single numerical value indicating the number of output variables (columns) to generate. Should be strictly >0 and an integer. 
#' @param col.length A single numerical value indicating the number of rows to generate in each column. Should be strictly >0 and an integer. 
#' @param FUNs A vector of character strings of length col.num corresponding to the sampling function that should be used to create each column. Options are: "sample", "rnorm", "rbinom", "rpois", "runif", "seq", and "rep". 
#' @param factor.YN A vector of logical values (TRUE or FALSE) of length col.num corresponding to whether each created column should be treated as a factor (TRUE). Defaults to FALSE, which will make all columns created numeric. 
#' @param ... Additional named arguments to be passed to the sampling functions as appropriate. Includes x.sample, replace, and sample.prob for sample(); max and min for runif(); mean and sd for rnorm(); size and binom.prob for rbinom(); to and from for seq(); x.rep and each for rep; and lambda for rpois(). If you are using a specific sampling function multiple times (i.e. for making multiple columns), inputs to these extra arguments will need to be supplied in a matrix with rows equal to the number of times you use the associated function (see examples). This enables you to use different input parameter values for two different columns drawn using the same function if you so choose. Alternative, each row can be the same to use the same inputs. 
#' @export
#' @details The idea of this function is to rapidly generate random columns of data based on sampling specifications you specify so that you can create an entire randomized data set (potentially with thousands of columns) in a single function call rather than column by column. It allows you to draw random data using a substantial range of standard randomization functions and even control which input parameters are used each instance a given randomization function is used if it's used more than once. As such, it may have multiple viable use cases. However, it's primary purpose is to generate a series of X data to be fed into sim.obs.data to then generate observed data based on some proposed relationship between the X data and observed data generated. See the documentation for the sim.obs.data function for details. 
#' 
#' @examples 
#' #This example code shows all the basic features of the function as well as the proper input syntax.
#' test = sim.pred.data(col.num = 14, col.length = 50, 
#'    FUNs = c("sample", "sample", "rnorm", "rnorm", "runif", "runif", "rbinom", "rbinom", 
#'       "rpois", "rpois", "seq", "seq", "rep", "rep"),
#'    factor.YN = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, 
#'        FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE),
#'    x.sample = matrix(c(0,1, 10,-10), nrow=2, byrow=T), 
#'    replace=matrix(c(TRUE, TRUE), nrow=2), 
#'    sample.prob=matrix(c(0.5, 0.5, 0.1, 0.9), nrow=2, byrow=T), 
#'    mean = matrix(c(0, 20), nrow=2), sd = matrix(c(5,1), nrow=2), 
#'    min = matrix(c(-10,0), nrow=2, byrow=T), max = matrix(c(0, 10), nrow=2, byrow=T),
#'    size = matrix(c(3, 20), nrow=2), binom.prob=matrix(c(0.5,0.1), nrow=2, byrow=T),
#'    lambda = matrix(c(0.5, 2), nrow=2),
#'    from = matrix(c(0, 100), nrow=2), to = matrix(c(-100, 500), nrow=2),
#'    x.rep = matrix(c(0,1,2,0,50, 100), nrow=2, byrow=T), each = matrix(c(1, 3), nrow=2))
#' print(test)
#' str(test)

##This function will help you generate a set of X data that can then be funneled into a call to sim.obs.data() or perhaps used for some other purpose like as X data for a ggplot.

sim.pred.data = function(col.num, col.length, FUNs, factor.YN = FALSE, ... ) {
  
#We first make an empty data frame to store our generated data. We also need to create an empty vector for column names.
pred.data = data.frame(matrix(NA, ncol=col.num, nrow=col.length))
col.names = numeric(0)

#To allow draws using the same function (e.g. mean) but with different parameters, we'll employ a counter system.
sample.counter = 0
rnorm.counter = 0
runif.counter = 0
rbinom.counter = 0
rpois.counter = 0
seq.counter = 0
rep.counter = 0

#We then have to pack any extra arguments provided to the ellipsis argument to more easily work with them. I got the following code here: 
#https://www.r-bloggers.com/2015/02/r-three-dots-ellipsis/
arguments = list(...)

#A very likely issue is that the extra arguments provided are not matrices but rather vectors. I will assume that this is because the user intends to use each sampling function just once or would rather just repeat the same instructions for each use of each sampling function, so I will coerce any vectors given to matrices with a single row and then force the function later to not allow inputs beyond the last one provided.
vectors.to.change = which(lapply(arguments, is.vector) == TRUE)
if (length(vectors.to.change) > 0) {
  for (i in vectors.to.change) {
    arguments[[i]] = matrix(arguments[[i]], nrow=1, byrow = T)
  }
}

#We should catch some predictable errors next...
if(col.num != length(FUNs)) { stop("Specify a sampling function for each desired output column (col.num)!")}
if(col.num <= 0 | col.num%%1 != 0) { stop("Make sure you specify a positive integer number of columns to generate!")}
if(col.length <= 0 | col.length%%1 != 0) { stop("Make sure you specify a positive integer number of row to generate!")}
if(any(FUNs == "sample") & (is.null(arguments$x.sample) | is.null(arguments$replace) | is.null(arguments$sample.prob))) { stop("To use sample, supply named values for x, replace, and sample.prob")}
if(any(FUNs == "rnorm") & (is.null(arguments$mean) | is.null(arguments$sd))) { stop("To use rnorm, supply named values for mean and sd")}
if(any(FUNs == "runif") & (is.null(arguments$min) | is.null(arguments$max))) { stop("To use runif, supply named values for min and max")}
if(any(FUNs == "rbinom") & (is.null(arguments$size) | is.null(arguments$binom.prob))) { stop("To use rbinom, supply named values for binom.prob and size")}
if(any(FUNs == "rpois") & is.null(arguments$lambda)) { stop("To use rpois, supply named values for lambda")}
if(any(FUNs == "seq") & (is.null(arguments$from) | is.null(arguments$to))) { stop("To use seq, supply named values for from and to")}
if(any(FUNs == "rep") & (is.null(arguments$each) | is.null(arguments$x.rep))) { stop("To use rep, supply named values for each and x.rep")}

#Next, we iterate over every column to be generated. For each one, we draw the appropriate amount and type of data randomly based on the desired function and specifications provided in the call. This section is only limited by our imagination in terms of what random-drawing functions might be implemented here.
  for (col in 1:col.num) {
    
    if(FUNs[col] == "seq") {
      seq.counter = seq.counter + 1
      temp1 = seq(from = arguments$from[max(seq.counter,dim(arguments$from)[1]),], length.out = col.length, to = arguments$to[max(seq.counter,dim(arguments$to)[1]),]) #x, replace, and sample.prob will need to be specified.
    }
    if(FUNs[col] == "rep") {
      rep.counter = rep.counter + 1
      temp1 = rep(x = arguments$x.rep[max(rep.counter, dim(arguments$x.rep)[1]),], length.out = col.length, each = arguments$each[max(rep.counter, dim(arguments$each)[1]),]) #x, replace, and sample.prob will need to be specified.
    }
    if(FUNs[col] == "sample") {
     sample.counter = sample.counter + 1
     temp1 = sample(x = arguments$x.sample[max(sample.counter, dim(arguments$x.sample)[1]),], size = col.length, replace = arguments$replace[max(sample.counter, dim(arguments$replace)[1]),], prob = arguments$sample.prob[max(sample.counter, dim(arguments$sample.prob)[1]),]) #x, replace, and sample.prob will need to be specified.
    }
    if(FUNs[col] == "rnorm") {
      rnorm.counter = rnorm.counter + 1
      temp1 = rnorm(col.length, mean = arguments$mean[max(rnorm.counter, dim(arguments$mean)[1]),], sd = arguments$sd[max(rnorm.counter, dim(arguments$sd)[1]),]) #mean and sd will need to be specified.
    }
    if(FUNs[col] == "runif") {
      runif.counter = runif.counter + 1
      temp1 = runif(col.length, min = arguments$min[max(runif.counter, dim(arguments$min)[1]),], max = arguments$max[max(runif.counter, dim(arguments$max)[1]),]) #min and max will need to be specified.
    }
    if(FUNs[col] == "rbinom") {
      rbinom.counter = rbinom.counter + 1
      temp1 = rbinom(col.length, size = arguments$size[max(rbinom.counter, dim(arguments$size)[1]),], prob = arguments$binom.prob[max(rbinom.counter, dim(arguments$binom.prob)[1]),]) #size and binom.prob will need to be specified.
    }
    if(FUNs[col] == "rpois") {
      rpois.counter = rpois.counter + 1
      temp1 = rpois(col.length, lambda = arguments$lambda[max(rpois.counter, dim(arguments$lambda)[1]),]) #lambda will need to be specified.
    }
   
  #If the end user wants the data in this column to be factor data, this is implemented next.  
    if (isTRUE(any(factor.YN))){
       if(factor.YN[col] == TRUE) {
      temp1 = factor(temp1)
      } 
    }
    
    if(is.null(temp1)) {stop ("Make sure all functions provided to FUNs are implemented!")}
    
  #The next three lines insert the generated data into the appropriate column in our storage data frame, pred.data, and then also names that column something predictable, like "X2" for the second column.  
    col.names[col] = paste0("X", col)
    pred.data[,col] = temp1
    temp1 = NULL #This is to ensure that a temp1 staying in memory doesn't hang around for the next loop.
  }
colnames(pred.data) = col.names
  
#We output our data frame for our predicted data.
  return(data.frame(pred.data))
}

