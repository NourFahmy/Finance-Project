library(quadprog)
library(Matrix)
senate_picks = read.csv("senate_picks.csv", header = T)
market_returns = read.csv("market_returns.csv", header = T)

sims = dim(market_returns)[1] # length of the data
stocks_count = dim(market_returns)[2]

#Simulate monthly efficient frontiers using past 2 years' data

plot(NULL, xlim=c(.00, .045),ylim=c(-0.022,.03),
     lwd=3,col="red", xlab = "SD of portfolio return",
     ylab = "mean of portfolio return")

for (start in seq(2, sims - 504, by=21)) {
  
  currentStocks = colnames(senate_picks[start,])[apply(senate_picks[start,], 2, function(pick) any(pick>0))]
  
  if (length(currentStocks)>2){
 
    senate = subset(x=market_returns, select=currentStocks)[start:(start+503),]
    sen_stocks = dim(senate)[2]
    
    mu = colMeans(senate)
    sigma = as.matrix(nearPD(cov(senate))$mat)
    
    m = 500 # no. of points to evaluate
    muP = seq(-0.02,0.02,length=m) # target portfolio return
    sdP = rep(0, length(muP)) # sd of portfolio return
    weight = matrix(0,nrow=m,ncol=sen_stocks) # storage for portfolio weights
    for (i in 1:length(muP)) { # find the optimal portfolios
      result = solve.QP(Dmat=2*sigma,dvec=rep(0,sen_stocks),
                        Amat = cbind(rep(1,sen_stocks),mu),bvec=c(1,muP[i]),meq=2)
      sdP[i] = sqrt(result$value)
      weight[i,] = result$solution
    }
    
    GMP = which.min(sdP) # global minimum point
    # efficient frontier
    lines(sdP[GMP:m],muP[GMP:m], type="l", lty = 1,lwd=3, col = "red")
    points(sdP[1:(GMP-1)],muP[1:(GMP-1)], type="l", lty = 2,lwd=3, col = "red")
  }
}

for (start in seq(2, sims - 504, by=21)) {
  market = market_returns[start:(start+503), 2:stocks_count]
  
  mu = colMeans(market)
  sigma = as.matrix(nearPD(cov(market))$mat)
  
  m = 500 # no. of points to evaluate
  muP = seq(-0.02,0.02,length=m) # target portfolio return
  sdP = rep(0, length(muP)) # sd of portfolio return
  weight = matrix(0,nrow=m,ncol=stocks_count-1) # storage for portfolio weights
  for (i in 1:length(muP)) { # find the optimal portfolios
    result = solve.QP(Dmat=2*sigma,dvec=rep(0,stocks_count-1),
                      Amat = cbind(rep(1,stocks_count-1),mu),bvec=c(1,muP[i]),meq=2)
    sdP[i] = sqrt(result$value)
    weight[i,] = result$solution
  }
  
  GMP = which.min(sdP) # global minimum point
  # efficient frontier
  lines(sdP[GMP:m],muP[GMP:m], type="l", lty = 1,lwd=3, col = "blue")
  points(sdP[1:(GMP-1)],muP[1:(GMP-1)], type="l", lty = 2,lwd=3, col = "blue")
}

legend("topleft", legend=c("Senators' Frontier", "Market Frontier"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
