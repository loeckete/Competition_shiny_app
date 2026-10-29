# Original implementation of the EcoEvoApps two-species LV competition equations.
default_config <- function() list(N1=20,N2=20,r1=.2,r2=.2,K1=100,K2=100,alpha=.5,beta=.5,time=150)
preset_config <- function(name) {
 p <- default_config()
 if(name=='B') {p$alpha <- .5;p$beta <- 1.5}
 if(name=='C') {p$alpha <- 1.5;p$beta <- .5}
 if(name=='D') {p$alpha <- 1.5;p$beta <- 1.5;p$N1 <- 30;p$N2 <- 20}
 if(name=='E') {p$alpha <- 1.5;p$beta <- 1.5;p$N1 <- 20;p$N2 <- 30}
 p
}
validate_config <- function(p) {
 keys <- names(default_config())
 if(!all(keys %in% names(p)) || any(!vapply(p[keys],function(x) is.numeric(x)&&length(x)==1&&is.finite(x),logical(1)))) stop('Enter a finite number in every field.')
 if(any(unlist(p[c('N1','N2')])<0)||any(unlist(p[c('N1','N2')])>2000)) stop('Initial abundances must be between 0 and 2,000.')
 if(min(p$r1,p$r2)<=0||max(p$r1,p$r2)>2) stop('Growth rates must be greater than 0 and at most 2.')
 if(min(p$K1,p$K2)<10||max(p$K1,p$K2)>2000) stop('Carrying capacities must be between 10 and 2,000.')
 if(min(p$alpha,p$beta)<0||max(p$alpha,p$beta)>3) stop('Competition coefficients must be between 0 and 3.')
 if(p$time<1||p$time>500) stop('Duration must be from 1 to 500 time units.')
 p
}
rhs <- function(y,p) c(N1=p$r1*y[[1]]*(1-(y[[1]]+p$alpha*y[[2]])/p$K1),
 N2=p$r2*y[[2]]*(1-(y[[2]]+p$beta*y[[1]])/p$K2))
classify <- function(p) {
 # Invasion growth factors, with the resident at its single-species K.
 a <- 1-p$alpha*p$K2/p$K1; b <- 1-p$beta*p$K1/p$K2
 eps <- 1e-9
 if(abs(a)<eps&&abs(b)<eps) return('Neutral boundary: coincident isoclines; a continuum of equilibria.')
 if(abs(a)<eps||abs(b)<eps) return('Boundary case: one invasion growth rate is zero. Inspect the trajectory carefully.')
 if(a>0&&b>0) return('Stable coexistence: each species can increase when rare against the other at its carrying capacity.')
 if(a<0&&b<0) return('Priority effects: neither species can invade the other at its carrying capacity. Starting populations can change the winner.')
 if(a>0) return('Species 1 excludes species 2 for positive starting populations in this constant environment.')
 'Species 2 excludes species 1 for positive starting populations in this constant environment.'
}
interior_equilibrium <- function(p) {
 den <- 1-p$alpha*p$beta
 if(abs(den)<1e-10) return(NULL)
 eq <- c(N1=(p$K1-p$alpha*p$K2)/den,N2=(p$K2-p$beta*p$K1)/den)
 if(any(eq<=0)) return(NULL)
 eq
}
logistic <- function(t,N,r,K) if(N==0) rep(0,length(t)) else K/(1+(K/N-1)*exp(-r*t))
simulate_model <- function(p,max_step=.1) {
 p <- validate_config(p)
 if(!is.finite(max_step)||max_step<=0) stop('Invalid integration step.')
 times <- sort(unique(c(seq(0,p$time,by=.25),p$time)))
 out <- matrix(0,length(times),2,dimnames=list(NULL,c('N1','N2')))
 out[1,] <- y <- c(N1=p$N1,N2=p$N2)
 u1 <- max(p$N1,p$K1);u2 <- max(p$N2,p$K2)
 rate <- max(p$r1*(1+(u1+p$alpha*u2)/p$K1),p$r2*(1+(u2+p$beta*u1)/p$K2))
 hmax <- min(max_step,.2/rate)
 if(p$time/hmax>150000) stop('This combination is too costly for a classroom browser. Shorten the duration or reduce the abundance/capacity contrast.')
 for(j in 2:length(times)) {
  n <- ceiling((times[j]-times[j-1])/hmax);h <- (times[j]-times[j-1])/n
  for(k in seq_len(n)) {
   k1 <- rhs(y,p);k2 <- rhs(y+h*k1/2,p);k3 <- rhs(y+h*k2/2,p);k4 <- rhs(y+h*k3,p)
   y <- y+h*(k1+2*k2+2*k3+k4)/6
   if(any(!is.finite(y))||any(y< -1e-8)) stop('Numerical integration failed.')
  }
  out[j,] <- y
 }
 d <- data.frame(time=times,out)
 d$alone1 <- logistic(times,p$N1,p$r1,p$K1);d$alone2 <- logistic(times,p$N2,p$r2,p$K2)
 list(config=p,data=d)
}
