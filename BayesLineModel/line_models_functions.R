#Functions to cluster 2-dimensional effects into groups defined
# by lines and correlations. 
# Matti Pirinen 5.4.2022

#Consider two effects:  x = effect1 and y = effect2.
#The line models are defined by three parameters
# scale, standard deviation of the (larger) effect
# slope, the slope of line around which the effects are scattered (where slope = Inf means y-axis)
# cor, the non-negative correlation between effects where cor = 1 means the same effect size

#In more detail the model has the following properties:
#
# 1) effects are scattered around line y = slope*x.
#    'slope' can be any real number or Inf (in which case effect x is zero)
#
# 2) the larger of the prior variances of effects is scale^2
#    In detail, if |slope| <= 1, then Var(x) = scale^2 
#               if |slope| > 1, then Var(y) = scale^2 
#
# 3) Correlation 'r' between effects is determined as follows.
#    Consider a distribution scattered around line y = x with correlation 'cor'. 
#    Rotate that distribution by an orthogonal rotation defined by angle
#    theta = atan(slope) - pi/4
#     and use the corresponding distribution, scaled so that maximum variance
#     is scale^2.
#    (NOTE: We do not use directly "correlation cor around line y = slope*x",
#           because the shape of that distribution depends on slope which we 
#           do not want here.)
#
# Examples:
# Set scale = 0 to get the null model. (Values of slope and cor do not matter in this case.) 
# Use scale > 0, slope = 1 and cor = 0 to get independent effects model.
# Use scale > 0, slope = 1 and cor = 1 to get fixed effects model.
# Use scale > 0, slope = 1 and 0 < cor < 1  to get correlated effects model.
#               note that values need to be quite near cor = 1 to get similar values
#               between the two effects with high probability (e.g. cor = 0.99 or 0.999).
#
# To choose scale, note that 95% of the effects will be < scale.

# There are three functions for user to call:
#
# visualize.line.models(scales, slopes, cors, 
#                       model.names = NULL, model.cols = NULL,
#                       legend.position = "bottomright")
# To visualize the lines and 95% highest probability regions of the models.

#line.models(X, SE, scales, slopes, cors, model.names = NULL, 
#            model.priors = rep(1/length(slopes), length(slopes)),r.lkhood = 0)
# To evaluate the model probabilities for each pair of effects separately.  

#line.models.with.proportions(X, SE, scales, slopes, cors, model.names = NULL,
#                             r.lkhood = 0, n.iter = 200, n.burnin = 20)
# To evaluate the model probabilities together with the proportions of effects pairs
# coming from each model. This is a joint estimation of groups across all pairs of effects.


prior.V <- function(scale = 0, slope = 1, cor = 0){
# Returns prior covariance matrix for line model defined by scale, slope and cor.
  
  if(cor > 1) stop("cor > 1")
  if(cor < 0) stop("cor < 0")
  if(scale < 0) stop("scale < 0")
  
  theta = atan(slope) - pi/4
  R = matrix(c(cos(theta), -sin(theta), # R is rotation matrix of angle theta
               sin(theta), cos(theta)), 
             ncol = 2, byrow = T)
  S = R %*% matrix(c(1, cor, cor, 1), ncol = 2) %*% t(R)
  S = scale^2 * S / max(as.vector(S))
  return(S)
 } 

sample.line.model <- function(n, scale, slope, cor, scale.weights){

  if(cor > 1) stop("cor > 1")
  if(cor < 0) stop("cor < 0")
  if(scale < 0) stop("scale < 0")
  if(!is.vector(scale.weights)) stop("scale.weights must be vector.")
  if(any(scale.weights < 0)) stop("scale.weights can have only positive values.")
  J = length(scale.weights)

  V = prior.V(scale = scale, slope = slope, cor = cor)
  
  #Choose components for each sample
  ind = sample(1:J, size = n, prob = scale.weights, replace = T) 
  
  if(abs(V[1]*V[4]-V[2]*V[3]) < 1e-16){A = matrix(c(sqrt(V[1]),sqrt(V[4]),0,0),ncol = 2)}
  else {A = t(chol(V))}
  x = t(A %*% matrix(rnorm(n*2), nrow = 2) / rep(J + 1 - ind, each = 2))
  cat(paste0("Sampling effects with scale=",scale," slope=",slope," cor=",cor),",\n")
  cat(paste0("   scale.weights=",paste(scale.weights,collapse = ","),"."),"\n")
  cat(paste0("Theoretical SD for larger effect:",
            signif(sqrt(sum(scale^2*scale.weights/sum(scale.weights)/rev(1:J)^2)),5)))
  cat(paste0("; observed value:",signif(max(apply(x,2,sd)),5)),".\n")
  return(x)
}

rdirichlet <- function(alpha){
  #Random sampling from Dirichlet distribution with parameter vector alpha
  g = rgamma(length(alpha), shape = alpha, scale = 1)
  return(g/sum(g))
}

log.dmvnorm <- function(x, mu = rep(0, length(x)), S = diag(1, length(x)) ){
  #returns log of density of MV-Normal(mean = mu, var = S) at x 
  K = length(mu)
  stopifnot(all(dim(S) == K))
  stopifnot(length(x) == K)
  chol.S = chol(S) #Cholesky decomposition
  log.det = 2*sum(log(diag(chol.S))) #log of det(S)
  inv.chol.S = solve(t(chol.S)) #inverse of cholesky^T
  return(-K/2*log(2*pi) - 0.5*(log.det + crossprod(inv.chol.S %*% as.numeric(x-mu))))
}


visualize.line.models <- function(scales, slopes, cors, 
                                  model.names = NULL, model.cols = NULL,
                                  legend.position = "bottomright",
                                  xlim = NULL, ylim = NULL, 
                                  xlab = "EFFECT1", ylab = "EFFECT2"){
  #Draws 95% highest probability regions of the models.
  #Note: Does NOT allow specification of 'scale.weights' 
  # but uses the given scales as the single component for each distribution.
  # (This is because HDR of a mixture of Gaussians is not simple to compute.)

  K = length(scales)
  lim = 3*max(scales) #Default: show models within 3 SDs
  if(is.null(xlim)) xlim = c(-lim,lim)
  if(is.null(ylim)) ylim = c(-lim,lim)
  plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  if(is.null(model.names)) model.names = paste0("Mod",1:K)
  if(is.null(model.cols)) model.cols = 2:K
  prob.level = 0.95
  b = qchisq(prob.level, df = 2)
  grid()
  abline(h = 0, lwd = 1)
  abline(v = 0, lwd = 1)

  for(ii in 1:K){
    if(scales[ii] < 1e-16) {
      points(0,0, col = model.cols[ii], pch = 19)
      next
    }
    if(!is.finite(slopes[ii])) {
      abline(v = 0, col = model.cols[ii], lty = 1)}
    else{
      abline(0, slopes[ii], col = model.cols[ii], lty = 1)}
    if(cors[ii] < 1){
      Sigma = prior.V(scale = scales[ii], slope = slopes[ii], cor = cors[ii])
      a = as.numeric(t(solve(Sigma)))
      x.lim = sqrt(abs(4*a[4]*b/( (a[2]+a[3])^2 - 4*a[4]*a[1])))
      x = seq(-x.lim*0.999, x.lim*0.999, length = 500) 
      y.upper = (-(a[2]+a[3])*x + sqrt((a[2]+a[3])^2*x^2 - 4*a[4]*(a[1]*x^2-b)))/2/a[4]
      y.lower = (-(a[2]+a[3])*x - sqrt((a[2]+a[3])^2*x^2 - 4*a[4]*(a[1]*x^2-b)))/2/a[4]
      lines(x, y.upper, lty = 2, lwd = 1.5, col = model.cols[ii])
      lines(x, y.lower, lty = 2, lwd = 1.5, col = model.cols[ii])
    }
  }
  if(!is.null(legend.position)){
    legend(legend.position, pch = 15, leg = model.names, col = model.cols)}
}

visualize.scales <- function(scales, scale.weights = c(1), cols = NULL){
  #Plots the univariate distributions of the effect sizes
  #Input parameters 'scales' and 'scale.weights' as in line.models()
  #Can take in user given colors but there must be equal number of colors as there
  # are models specified by scales vector.
  
  x.ran = 2.5 * max(scales)
  K = length(scales)
  stopifnot(all(scale.weights >= 0))
  
  if(is.vector(scale.weights)){
    scale.weights = matrix(scale.weights, byrow = T,
                              ncol = length(scale.weights), nrow = K)
  }
  stopifnot(is.matrix(scale.weights))
  stopifnot(nrow(scale.weights) == K)
  scale.weights = scale.weights / rowSums(scale.weights)
  n.comps = ncol(scale.weights)
  x = seq(-x.ran, x.ran, length = 1000)
  if(is.null(cols)) cols = c(2:(K+1))
  if(length(cols) != K) stop(paste("Length of vector of cols must be number of models, here",K))
  
  y.max = 0
  print(paste("Plotting the following effect distributions:"))
  for(draw in c(FALSE, TRUE)){
    if(draw) {
      plot(NULL, xlim = c(-x.ran,x.ran), ylim = c(0, y.max),
           xlab = "effect", ylab = "density")
      legend("topleft", col = cols, leg = paste(scales), pch = 15)}
    for(kk in 1:K){
      y = rep(0, length(x))
      for(jj in 1:n.comps){
        y = y + scale.weights[kk,jj]*dnorm(x, 0, scales[kk]/(n.comps - jj + 1))}
      if(draw) {
        lines(x, y, col = cols[kk], lwd = 1.5)
        if(scales[kk] < 1e-300) points(0,0, col = cols[kk], pch =  19)
        print(paste("Scale:",scales[kk]," weights:",paste(scale.weights[kk,],collapse = ",")))}
      y.max= max(c(y ,y.max))
    }
  }
}  


line.models <- function(X, SE, 
                        scales, slopes, cors, 
                        model.names = NULL, 
                        model.priors = rep(1/length(slopes), length(slopes)),
                        r.lkhood = 0, scale.weights = c(1)){
  #Evaluates model probabilities for each data point separately.
  #You can name the models in output using 'model.names'
  
  #Default prior probabilities for models are equal, 1/K, where K is the number of models.
  
  #Default assumption are that two sets of effects are independent. If that is not the case,
  #Use 'r.lkhood' to determine the correlation in the likelihood function.
  
  #Effect distribution for larger effect can be specified as a mixture of Gaussians.
  # If scale.weights is a vector of length k > 1, then model 'm' is a mixture
  # of k Gaussians: 
  # sum_{i=1}^k scale.weights[i]*N(0, prior.V(scales[m]/(k-i+1), slopes[m], cols[m])) 
  # where vector scale.weights has been normalized to sum to 1.
  # If scale.weights is a matrix with k > 1 columns and K rows, 
  # then model 'm' is a mixture of k Gaussians: 
  #  sum_{i=1}^k scale.weights[m,i]*N(0,scales[m]/(k-i+1)) 
  # where rows of matrix scale.weights has been normalized to sum to 1.
  # Default value scale.weights = c(1) means that each model has only one component
  #  whose standard deviation is 'scales[m]'.
  
  #Returns the membership probabilities for each data point.
  
  K = length(slopes) #number of models
  stopifnot(length(scales) == K)
  stopifnot(length(cors) == K)
  stopifnot(r.lkhood >= (-1) && r.lkhood <= 1)
  stopifnot(all(cors <= 1) && all(cors >= 0))
  stopifnot(ncol(X) == 2)
  stopifnot(ncol(X) == 2)
  stopifnot(nrow(X) == nrow(SE))
  stopifnot(all(model.priors >= 0))
  stopifnot(all(scale.weights >= 0))
  
  pis = model.priors/sum(model.priors)
  if(is.vector(scale.weights)){
    scale.weights = matrix(scale.weights, byrow = T,
                              ncol = length(scale.weights), nrow = K)
  }
  stopifnot(is.matrix(scale.weights))
  stopifnot(nrow(scale.weights) == K)
  n.comps = ncol(scale.weights)
  
  pr.V = list() #list of lists of prior matrices of each model and each scale
  for(kk in 1:K) { #over models
    w = scale.weights[kk,]
    w = w/sum(w)
    model.Vs = list()
    for(jj in 1:n.comps){ #over scales
      model.Vs[[jj]] = prior.V(scale = scales[kk]/(n.comps + 1 - jj), 
                               slope = slopes[kk], cor = cors[kk])}
    pr.V[[kk]] = model.Vs 
  }
  n = nrow(X)
  
  logdnorm = matrix(NA, ncol = K, nrow = n)
  for(ii in 1:n){
    V = diag(SE[ii,]) %*% 
      matrix(c(1, r.lkhood, r.lkhood, 1), 2, 2, byrow = T) %*%
      diag(SE[ii,]) #var of likelihood
    scale.weights
    for(kk in 1:K){
      tmp = sapply(1:n.comps, function(jj){
        log.dmvnorm(X[ii,], mu = c(0,0), S = V + pr.V[[kk]][[jj]])})
      tmp = tmp + log(scale.weights[kk,])
      logdnorm[ii,kk] = log(sum(exp(tmp)))                
    }
  }
  p = exp(t(t(logdnorm) + log(pis)))
  res = p/rowSums(p)
  colnames(res) = model.names
  return(res)
}

  
line.models.with.proportions <- 
  function(X, SE,  
           scales, slopes, cors,  
           model.names = NULL,
           r.lkhood = 0,
           n.iter = 200, n.burnin = 20){
  #Gibbs sampler to estimate the model probabilities for each data point 
  # and the proportions of data points in each group.
  #You can name the models in output using 'model.names'
  #Default assumption is that two sets of effects are independent. If that is not the case,
  #use 'r.lkhood' to determine the correlation in the likelihood function.
  #Prior for models is Dirichlet(1/K,...,1/K), where K is the number of models.  
        
  #Returns the membership probabilities for each data point ('groups')
  # and also posterior distribution of the proportions ('params').
    

  K = length(slopes) #number of models
  stopifnot(length(scales) == K)
  stopifnot(length(cors) == K)
  stopifnot(r.lkhood >= (-1) && r.lkhood <= 1)
  stopifnot(all(cors <= 1) && all(cors >= 0))
  stopifnot(ncol(X) == 2)
  stopifnot(ncol(X) == 2)
  stopifnot(nrow(X) == nrow(SE))
  stopifnot(n.burnin >= 0)
  stopifnot(n.iter > 0)
  
  pr.V = list()
  for(kk in 1:K) pr.V[[kk]] = prior.V(scale= scales[kk], slope = slopes[kk], cor = cors[kk])
  n = nrow(X)
  
  logdnorm = matrix(NA, ncol = K, nrow = n)
  for(ii in 1:n){
    V = diag(SE[ii,]) %*% 
      matrix(c(1, r.lkhood, r.lkhood, 1), 2, 2, byrow = T) %*%
      diag(SE[ii,]) #var of likelihood
     for(kk in 1:K){
        logdnorm[ii,kk] = log.dmvnorm(X[ii,], mu = c(0,0), S = V + pr.V[[kk]])
     }
  }
  
  diri.a = rep(1/K, K) #Dirichlet prior parameters

  R.par = matrix(NA, ncol = K, nrow = n.iter) #results for parameters
  R.ind = matrix(0, ncol = K, nrow = n) #results for individual effects
  colnames(R.par) = colnames(R.ind) = model.names

  pis = rep(1,K)/K #initial probabilities, vector pi
  grs = sample(1:K, size = n, prob = pis, replace = T) #initially random grouping among variants

  for(ii in 1:(n.burnin + n.iter)){ #Gibbs sampler 
    
      #count how many are assigned to each group      
      tmp = rle(sort(grs))
      gr.counts = rep(0,K)
      gr.counts[tmp$values] = tmp$length
      
      #sample pis
      pis = rdirichlet(diri.a + gr.counts)
      log.pis = log(pis)
      
      #sample group indicators for variants
      p = exp(t(t(logdnorm) + log.pis)) #no need to normalize for 'sample()'
      grs = apply(p, 1, function(pr){sample(1:K, size = 1, prob = pr)})
      if(ii > n.burnin){ #save results if burnin is over
        R.ind[(1:n) + n*(grs-1)] = 1 + R.ind[(1:n) + n*(grs-1)]
        R.par[ii-n.burnin,] = pis
      }
      if(ii %% 100 == 0) print(paste(ii,"iterations done."))
  }

  res.par = cbind(apply(R.par, 2, mean),
                t(apply(R.par ,2 , function(x){quantile(x, c(0.025, 0.975))})),
                apply(R.par, 2, sd))
  colnames(res.par) = c("mean","95%low","95%up","sd")

  res.ind = R.ind / n.iter
  return(list(params = res.par, groups = res.ind))
}

