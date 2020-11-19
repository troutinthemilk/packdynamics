library(RColorBrewer)

##A) larger packs gain access to structured prey
beta <- 2
lambda <- .2
k	<- c(0.01, 0.1, 0.5, 1, 2, 5)
resp.vec	<- list()
approx.vec	<- resp.vec
G	<- 5
P	<- 100
b	<- sqrt(k/(k-1))*beta/lambda
W	<- seq(0,10000,length.out=1e4)

beta*W/(lambda*G)

f.response <- function(P, W, G, beta, lambda, k) { P/G*(1  - exp(-(beta*W/(G*lambda))^k))}
f.approx <- function(P, W, G, beta, lambda, k) { P/G*(beta*W/(G*lambda))^k}

for(i in 1:length(k)) {
  resp.vec[[i]] <- f.response(P, W, G, beta, lambda, k[i])
  approx.vec[[i]] <- f.approx(P, W, G, beta, lambda, k[i])
}

col.vec <- brewer.pal(length(k)+1, 'YlOrBr') #'YlGnBu')
col.vec <- col.vec[-1]
par(mfrow=c(1,2))
plot(W, G*resp.vec[[1]], type='l', lwd=3, xlab='Predator density', ylab='Number prey killed', ylim=c(0,100), col=col.vec[1])
for(i in 2:length(k)) {
  lines(W, G*resp.vec[[i]], type='l', lwd=3, col=col.vec[i])
#   lines(W, approx.vec[[i]], type='l', lwd=1, lty=2, col=col.vec[i])
}

plot((W), resp.vec[[1]]/W, type='l', lwd=3, xlab='Predator density', ylab='Number prey killed per predator', ylim=c(1e-7, 1e2), col=col.vec[1], log='y')
for(i in 2:length(k)) {
  lines((W), resp.vec[[i]]/W, type='l', lwd=3, col=col.vec[i])
#   lines(W, approx.vec[[i]], type='l', lwd=1, lty=2, col=col.vec[i])
}
# lines(W, log(P/(1 + .2*P + 10*W)))

legend('topright', legend=k, col=col.vec, lwd=3, title='k')

##C) carrying capacity in G. 
graphics.off()
W	<- seq(0,1000)
a	<- 
f	<- function(G, W, lambda, beta, k) { 
  arg <- (beta*W/(G*lambda))^k
  lambda*G*(1 - exp(-arg)) - beta*W*exp(-arg) 
}

f.powerlaw	<- function(G, W, lambda, beta, k) { 
  arg <- (beta*W/(G*lambda))
  rho <- beta
  P/W*(1 - (1 + arg)^(-k)) - 1/(lambda*G)*(1+arg)^(-k-1) 
}

f.orig <- function(G, W=W, lambda=lambda, beta=beta, k=5) { 
  arg <- (beta*W/(G*lambda))^k
  G*P/W*(1 - exp(-arg))
}

k.val <- 5
plot(W, f(G=W, W=1e2, lambda=lambda, beta=beta, k=k.val), type='l', lwd=2, xlab='G', ylab='Functional response')
Gstar.vec	<- vector('numeric', length(W))
for(i in 2:length(W)) {
  test <- uniroot(f, c(0.0001,1e4), W=W[i], lambda=lambda, beta=beta, k=k.val)
#   test <- uniroot(f.powerlaw, c(0.001,1e4), W=W[i], lambda=lambda, beta=beta, k=0.01)
  Gstar.vec[i] <- test$root
}

plot(W, Gstar.vec, type='l', lwd=2, ylab='Number of groups')
a <- sqrt(k.val/(k.val-1))*beta/lambda
lines(W, a*W, col='red', lty=2)

((b*u*integrate(e^(d*t - 2*t*u + b*e^(-t*u) - b), t)/(d - u) + e^(d*t -
t*u + b*e^(-t*u) - b)/(d - u))*g + c)*e^(-d*t)