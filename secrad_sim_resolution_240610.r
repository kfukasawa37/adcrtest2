###########1/2 resolution scenario
###load packages and data
library(rflsgen)
library(tidyverse)
library(foreach)
library(doParallel)
library(ooplah)



load("secrad_sim_220811.Rdata")	


sourcepath<-"secrad.r"
source(sourcepath)

###resolution down
xcoord0<-floor((xcoord-0.5)/2)*2+1.5
ycoord0<-floor((ycoord-0.5)/2)*2+1.5

for(i in 1:niter){
	temp<-bind_cols(x=xcoord0,y=ycoord0,habitat=habitat[[i]])%>%
			group_by(x,y)%>%summarize(habitat=mean(habitat),.groups="drop")
	habitat[[i]]<-temp$habitat
}

xcoord<-temp$x
ycoord<-temp$y

ncell<-length(xcoord)
nx<-length(unique(xcoord))
ny<-length(unique(ycoord))

dx<-2
dy<-2
area<-dx*dy

#####trap
trapx<-floor((trapx-0.5)/2)*2+1.5
trapy<-floor((trapy-0.5)/2)*2+1.5	

effort_loc<-rep(NA,ntrap)
for(i in 1:ntrap){
	effort_loc[i]<-which(xcoord==trapx[i]&ycoord==trapy[i])
}

############prepare dataset
range.dpar<-c(-3,-1)
range.c1<-c(-2,2)
range.apar<-c(-1,0)
g0par<--3
c0par<-1

sim_list0<-vector("list",niter)
for(i in 1:niter){
	sim_list0[[i]]<-sim_list[[i]]$clone()
}

sim_list<-vector("list",niter)

nt<-5000

set.seed(8128)

for(i in 1:niter){
	cat(i,"\n")
	sim_list[[i]]<-secrad_data$new(coords=cbind(x=xcoord,y=ycoord),
							area=rep(4,ncell),
							grid_cov=data.frame(X=habitat[[i]]),
							resolution=c(x=2,y=2))
	sim_list[[i]]$add_obs(type="poisson",effort=rep(nt,ntrap),effort_loc=effort_loc,effort_occ=rep(1,ntrap),detect=sim_list0[[i]]$obs[[1]]$detect)
	sim_list[[i]]$set_truemodel(envmodel=list(D~1,C~X,A~1),
						indmodel=c(A=FALSE,g0=FALSE),
						occmodel=c(A=FALSE,g0=FALSE),
						dpar=sim_list0[[i]]$truemodel$dpar,
						cpar=sim_list0[[i]]$truemodel$cpar,
						adpar=sim_list0[[i]]$truemodel$adpar,
						g0par=sim_list0[[i]]$truemodel$g0par)
}


#####################estimation

envmodel<-list(D~1,C~X,A~0)
indmodel<-c(A=FALSE,g0=FALSE)
occmodel<-c(A=FALSE,g0=FALSE)
secrad_list<-vector("list",niter)
for(i in 1:niter){
	secrad_list[[i]]<-secrad$new(secrdata=sim_list[[i]])
	secrad_list[[i]]$set_model(envmodel=envmodel,indmodel=indmodel,occmodel=occmodel)
}


res_list<-vector("list",niter)
initpar<-c(-2,1.5,0,-3)
for(i in 1:niter){
	cat(i,"\n")
	cat(sim_list[[i]]$out_truepar(),"\n")
	res_list[[i]]<-try(optim(initpar,secrad_list[[i]]$loglf,method="BFGS",control=list(maxit=1000,trace=2),loglfscale=-1,verbose=F,hessian=T))
	cat(res_list[[i]]$par,"\n")
}

#save.image("secrad_sim_estim_reso_220811.Rdata")	



##########Postprocessing the results
library(rflsgen)
library(tidyverse)
library(foreach)
library(doParallel)
library(ooplah)



truepar<-sapply(sim_list,function(x) x$out_truepar())

truepar_estim<-truepar
truepar_estim[2,]<-truepar[2,]-truepar[4,]
truepar_estim<-truepar_estim[-4,]
par_estim<-sapply(res_list,"[[","par")
par_estim[2,]<-par_estim[2,]

which.max(abs(truepar_estim[2,]-par_estim[2,]))

plot(truepar_estim[1,],par_estim[1,])
points(c(-100,100),c(-100,100),type="l")


hessian_estim<-lapply(res_list,"[[","hessian")
par_se<-hessian_estim%>%lapply(solve)%>%sapply(diag)%>%sqrt

res_tidy<-tibble(iter=rep(1:niter,each=4),parname=rep(c("dens","con0","con1","g0"),100),true=c(truepar_estim),mle=c(par_estim),se=c(par_se))%>%
		mutate(lci=mle+qnorm(0.025,0,1)*se,uci=mle+qnorm(0.975,0,1)*se,bias=mle-true)
		
res_wide<-res_tidy%>%pivot_wider(id_cols="iter",names_from=parname,values_from=c("true","mle","se","lci","uci","bias"))%>%
		bind_cols(t(truepar))
lmres_dens<-lm(bias_dens~dens0+con0+con1+adv0,data=res_wide,weight=1/(se_dens)^2)
lmres_con0<-lm(bias_con0~dens0+con0+con1+adv0,data=res_wide,weight=1/(se_con0)^2)
lmres_con1<-lm(bias_con1~dens0+con0+con1+adv0,data=res_wide,weight=1/(se_con1)^2)
lmres_g0<-lm(bias_g0~dens0+con0+con1+adv0,data=res_wide,weight=1/(se_g0)^2)


ggplot(res_wide,aes(x=con0,y=bias_dens))+geom_point()+theme_cowplot()
ggplot(res_wide,aes(x=con0,y=bias_con0))+geom_point()+theme_cowplot()
ggplot(res_wide,aes(x=adv0,y=bias_con0))+geom_point()+theme_cowplot()
ggplot(res_wide,aes(x=con0,y=bias_con1))+geom_point()+theme_cowplot()
ggplot(res_wide,aes(x=con0,y=bias_g0))+geom_point()+theme_cowplot()

res_tidy<-res_tidy%>%mutate(parname2=factor(res_tidy$parname,levels=c("con0","con1","dens","g0"),
		labels=c('"ln(DC ratio)"','"Effect of landscape"',"'ln(density)'","ln(italic('g'[0]))")))
#labels<-as_labeller(c(`con0`='ln(DC ratio)',`con1`='Effect of landscape on connectivity',`dens`=label_bquote(ln(italic(D))),`g0`='ln(italic(g[0]))'))
plot_estim<-ggplot(res_tidy%>%filter((parname!="g0")),aes(x=true,y=mle,col=parname))+geom_point(size=1)+
		geom_errorbar(aes(x=true,ymax=uci,ymin=lci),size=0.3)+
		geom_abline(intercept = 0, slope = 1)+
		facet_wrap(~parname2,scales = "free",labeller=label_parsed,ncol=2)+
		theme_cowplot()+
		xlab("True")+ylab("Estimates")+
		theme(legend.position="none")
		
plot_estim;ggsave(paste0("plot_estim_reso_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=16,height=16,units="cm")

#save.image("secrad_sim_reso_240610.Rdata")
