###################load packages and data
library(rflsgen)
library(tidyverse)
library(foreach)
library(doParallel)
library(ooplah)
library(cowplot)
library(raster)
library(gdistance)
library(sp)
library(vegan)

load("secrad_sim_220811.Rdata")	


sourcepath<-ifelse(version$os=="mingw32","L:/secradsim/secrad.r","/L/secradsim/secrad.r")
source(sourcepath)

#####resolution down
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


##########Postprocessing the results

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

res_tidy<-res_tidy%>%mutate(parname2=factor(res_tidy$parname,levels=c("con0","con1","dens","g0"),
		labels=c('"ln(DC ratio)"','"Effect of landscape"',"'ln(density)'","ln(italic('g'[0]))")))

plot_estim<-ggplot(res_tidy%>%filter((parname!="g0")),aes(x=true,y=mle,col=parname))+geom_point(size=1)+
		geom_errorbar(aes(x=true,ymax=uci,ymin=lci),size=0.3)+
		geom_abline(intercept = 0, slope = 1)+
		facet_wrap(~parname2,scales = "free",labeller=label_parsed,ncol=2)+
		theme_cowplot()+
		xlab("True")+ylab("Estimates")+
		theme(legend.position="none")
		
plot_estim;ggsave(paste0("plot_estim_reso_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=16,height=16,units="cm")

###SCR-LCP

source("SCRed.pois.r")

truepar<-sapply(sim_list,function(x) x$out_truepar())

truepar_estim<-truepar
truepar_estim[2,]<-truepar[2,]-truepar[4,]
truepar_estim<-truepar_estim[-4,]
par_estim<-sapply(res_list,"[[","par")

SCRedreslist<-vector("list",niter)

for(i in 1:niter){
	SCRedinit<-c(-3-log(2*pi)-truepar_estim[2,i]+log(nt)-log(4),-log(2)-truepar_estim[2,i],log(exp(truepar_estim[1,i])*ncell*area-sim_list[[i]]$nind),0)
	sim_i<-sim_list[[i]]
	tdf1<-data.frame(De=1:ntrap,
				X=sim_i$coords[sim_i$obs[[1]]$effort_loc,"x"],
				Y=sim_i$coords[sim_i$obs[[1]]$effort_loc,"y"],
				Cov=sim_i$obs[[1]]$effort)
	
	
	ss<-data.frame(X=sim_i$coords[,"x"],Y=sim_i$coords[,"y"],Hab=sim_i$grid_cov$X)
	coordinates(ss)<- ~X+Y
 	gridded(ss) <- TRUE
 	r <- raster(ss, "Hab") 
 	XY<-as.matrix(coordinates(r))
 	cat(i,"\n")
	SCRedreslist[[i]]<-try(optim(SCRedinit,SCRed.pois, y=t(sim_i$obs[[1]]$detect),X=tdf1[,c("X","Y")],cov=r,G=XY,ssbuffer=8,directions=4,hessian=T,control=list(maxit=1000)))

}
	
#retry for failed cases
failed<-which((SCRedreslist%>%sapply(class))=="try-error")
maxtry<-100
for(i in failed){
	sim_i<-sim_list[[i]]
	tdf1<-data.frame(De=1:ntrap,
				X=sim_i$coords[sim_i$obs[[1]]$effort_loc,"x"],
				Y=sim_i$coords[sim_i$obs[[1]]$effort_loc,"y"],
				Cov=sim_i$obs[[1]]$effort)
	
	
	ss<-data.frame(X=sim_i$coords[,"x"],Y=sim_i$coords[,"y"],Hab=sim_i$grid_cov$X)
	coordinates(ss)<- ~X+Y
 	gridded(ss) <- TRUE
 	r <- raster(ss, "Hab") 
 	XY<-as.matrix(coordinates(r))
 	loglf<-Inf
 	count<-0
	while(is.infinite(loglf)&(count<maxtry)){
		count<-count+1
		SCRedinit_j<-runif(4,-2,2)
		loglf<-SCRed.pois(SCRedinit_j,y=t(sim_i$obs[[1]]$detect),X=tdf1[,c("X","Y")],cov=r,G=XY,ssbuffer=11)
	}
 	cat(i,"\n")
	SCRedreslist[[i]]<-try(optim(SCRedinit_j,SCRed.pois, y=t(sim_i$obs[[1]]$detect),X=tdf1[,c("X","Y")],cov=r,G=XY,ssbuffer=11,hessian=T,control=list(maxit=1000)))
}

#save.image("secrad_sim_SCRed_reso_220812.Rdata")

SCRed_code<-SCRedreslist%>%sapply("[[","convergence")
SCRed_con<-(-1)*unlist(SCRedreslist%>%lapply("[[","par")%>%lapply("[",4))
SCRed_dens<-log(exp(unlist(SCRedreslist%>%lapply("[[","par")%>%lapply("[",3)))+sapply(sim_list,"[[","nind"))-log(ncell)-log(area)


SCRed_covmat_list<-SCRedreslist%>%
				lapply("[[","hessian")%>%
				lapply(function(...){try(solve(...))})

SCRed_validse<-apply((SCRed_covmat_list%>%sapply(diag))>0,2,all)
SCRed_valid<-SCRed_validse&(SCRed_code==0)			

SCRed_con_se<-SCRed_covmat_list%>%
				lapply(function(...){try("["(...))},4,4)%>%
				lapply(function(x){ifelse(is.numeric(x),x,NA)})%>%
				unlist%>%
				sqrt

SCRed_dens_se<-SCRed_covmat_list%>%
				lapply(function(...){try("["(...))},3,3)%>%
				lapply(function(x){ifelse(is.numeric(x),x,NA)})%>%
				unlist%>%
				sqrt


SCRed_tibble<-tibble(true_con=truepar_estim[3,],SCRed_con=SCRed_con,SCRed_con_se=SCRed_con_se,SCRed_valid=SCRed_valid,true_dens=truepar_estim[1,],SCRed_dens=SCRed_dens,SCRed_dens_se=SCRed_dens_se,secrad_dens=par_estim[1,])%>%
			mutate(true_con_pos=true_con>=0,SCRed_con_pos=SCRed_con>=0,SCRed_dens_bias=SCRed_dens-true_dens,secrad_dens_bias=secrad_dens-true_dens,
					SCRed_n0=exp(SCRed_dens)*ncell-sapply(sim_list,"[[","nind"),true_n0=sim_list0%>%lapply(function(x) ncol(x$sim_result$detect[[1]])-x$nind)%>%unlist,secrad_n0=exp(secrad_dens)*ncell*area-sapply(sim_list,"[[","nind"),
					SCRed_con_lci=qnorm(0.025,SCRed_con,SCRed_con_se),SCRed_con_uci=qnorm(0.975,SCRed_con,SCRed_con_se),
					SCRed_dens_lci=qnorm(0.025,SCRed_dens,SCRed_dens_se),SCRed_dens_uci=qnorm(0.975,SCRed_dens,SCRed_dens_se))
			
plot_SCRLCP_cost<-ggplot()+geom_point(data=SCRed_tibble%>%filter(SCRed_valid),aes(x=true_con,y=SCRed_con),size=1)+
		geom_errorbar(data=SCRed_tibble%>%filter(SCRed_valid),aes(x=true_con,ymax=SCRed_con_uci,ymin=SCRed_con_lci),size=0.2)+
		theme_cowplot()+
		xlab("True effect of landscape\n on connectivity")+
		ylab("(-1)×effect of landscape\n on cost by SCR-LCP")

plot_SCRLCP_cost;ggsave(paste0("plot_SCRLCP_cost_reso_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=9,height=8,units="cm")


plotdata<-SCRed_tibble%>%filter(SCRed_valid)%>%select(true_dens,SCRed_dens,SCRed_dens_lci,SCRed_dens_uci)%>%mutate(model="SCRed")
colnames(plotdata)[2:4]<-c("mle_dens","lci_dens","uci_dens")		
plotdata<-bind_rows(plotdata,res_wide%>%select(true_dens,mle_dens,lci_dens,uci_dens)%>%mutate(model="secrad"))

plot_SCReddens<-ggplot(data=plotdata%>%filter(model=="SCRed"))+geom_point(aes(x=true_dens,y=mle_dens),size=1)+
	geom_abline()+
	geom_errorbar(aes(x=true_dens,ymin=lci_dens,ymax=uci_dens),size=0.2)+
	theme_cowplot()+
	xlab("True ln(density)")+
	ylab("Ln(density) estimated\n by SCR-LCP")
	
plot_SCReddens;ggsave(paste0("plot_dens_SCRLCP_reso_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=9,height=8,units="cm")
