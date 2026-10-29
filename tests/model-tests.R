source('app/model.R')
check <- function(x,msg) if(!isTRUE(x))stop(msg)
ref <- read.csv('tests/reference.csv')
for(case in LETTERS[1:5]) {
 p<-preset_config(case);r<-simulate_model(p);d<-r$data
 rows<-ref[ref$case==case,]
 idx<-match(rows$time,d$time)
 check(max(abs(as.matrix(d[idx,c('N1','N2')])-as.matrix(rows[c('N1','N2')])))<1e-5,paste('Reference',case))
 check(min(d$N1,d$N2)>=0,'Nonnegative populations')
 check(all(d$N1<=d$alone1+1e-7)&&all(d$N2<=d$alone2+1e-7),'Competition compared with absence')
 fine<-simulate_model(p,max_step=.025)$data
 check(max(abs(as.matrix(fine[c('N1','N2')])-as.matrix(d[c('N1','N2')])))<1e-5,'Step refinement')
}
p<-default_config();p$alpha<-p$beta<-0;d<-simulate_model(p)$data
check(max(abs(d$N1-d$alone1))<1e-5,'Independent logistic solution')
p<-default_config();p$N1<-0;d<-simulate_model(p)$data
check(all(d$N1==0),'No spontaneous recolonization')
p<-default_config();check(grepl('Stable coexistence',classify(p)),'Coexistence classification')
check(max(abs(interior_equilibrium(p)-c(100/1.5,100/1.5)))<1e-8,'Interior equilibrium')
check(grepl('Species 1 excludes',classify(preset_config('B'))),'Exclusion 1')
check(grepl('Species 2 excludes',classify(preset_config('C'))),'Exclusion 2')
check(grepl('Priority effects',classify(preset_config('D'))),'Priority effects')
p$alpha<-p$beta<-1;check(grepl('Neutral boundary',classify(p)),'Coincident isoclines')
p$beta<-.5;check(grepl('Boundary case',classify(p)),'Boundary case')
p<-preset_config('D');p$N1<-p$N2<-20;d<-simulate_model(p)$data
check(max(abs(d$N1-d$N2))<1e-9,'Symmetric unstable trajectory')
p$N1<-21;p$time<-500;d<-simulate_model(p)$data
check(tail(d$N1,1)>99&&tail(d$N2,1)<.01,'Perturbation breaks unstable balance')
p<-default_config();p$K1<-0
check(inherits(try(simulate_model(p),silent=TRUE),'try-error'),'Invalid K')
p<-default_config();p$alpha<-NA_real_
check(inherits(try(simulate_model(p),silent=TRUE),'try-error'),'Missing value')
invisible(parse('app/app.R'));cat('Model tests passed.\n')
