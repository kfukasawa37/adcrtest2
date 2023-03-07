##############################################
##########Simulation study
##############################################
library(rflsgen)
library(tidyverse)
library(foreach)
library(doParallel)
library(ooplah)
library(cowplot)
library(R6)
library(Matrix)
library(cowplot)
library(sf)
library(secr)
library(raster)
library(gdistance)
library(sp)
library(vegan)
library(viridis)
library(smatr)


############generating neutral landscape
niter<-100

nx<-50
ny<-50
ncell<-nx*ny
xcoord<-rep(1:nx,each=ny)
ycoord<-rep(1:ny,nx)


cls1<-flsgen_create_class_targets(
  "habitat",
  NP = c(5,10),
  AREA = c(100,400),
  CA = c(1000,2000),
  MESH = c(0,200),
)

target<-flsgen_create_landscape_targets(ny,nx,list(cls1))

habitat<-target%>%flsgen_structure(nb_solutions = niter,time_limit=100000)%>%
		lapply(flsgen_generate)%>%
		lapply(as.matrix)%>%
		lapply(c)%>%
		lapply("+",0.5)	#continuous habitat:-0.5,patchy habitat:0.5


#save.image("habitat100_220811.Rdata")


sourcepath<-"secrad.r"
source(sourcepath)

##############generating trapping efforts


trapx<-rep(seq(9,41,4),9)
trapy<-rep(seq(9,41,4),each=9)
ntrap<-length(trapx)
effort_loc<-rep(NA,ntrap)
for(i in 1:ntrap){
	effort_loc[i]<-which(xcoord==trapx[i]&ycoord==trapy[i])
}

############simulating dataset
range.dpar<-c(-3,-1)
range.c1<-c(-2,2)
range.apar<-c(-1,0)
g0par<--3
c0par<-1
sim_list<-vector("list",niter)

nt<-5000

set.seed(8128)

for(i in 1:niter){
	cat(i,"\n")
	sim_list[[i]]<-secrad_data$new(coords=cbind(x=xcoord,y=ycoord),
							area=rep(1,ncell),
							grid_cov=data.frame(X=habitat[[i]]),
							resolution=c(x=1,y=1))
	sim_list[[i]]$add_obs(type="poisson",effort=rep(nt,ntrap),effort_loc=effort_loc,effort_occ=rep(1,ntrap))
	sim_list[[i]]$set_truemodel(envmodel=list(D~1,C~X,A~1),
						indmodel=c(A=FALSE,g0=FALSE),
						occmodel=c(A=FALSE,g0=FALSE),
						dpar=runif(1,range.dpar[1],range.dpar[2]),
						cpar=c(c0par,runif(1,range.c1[1],range.c1[2])),
						adpar=runif(1,range.apar[1],range.apar[2]),
						g0par=g0par)
	sim_list[[i]]$set_sim(nt,1,100,100)
	sim_list[[i]]$generate_ind()
}

cluster = makeCluster(detectCores()-1, type = "FORK")
registerDoParallel(cluster)
sim_list2<-foreach(i=1:niter) %dopar% {
	cat(i,"\n")
	sim_list[[i]]$simulate(return=T,verbose=F)
}
stopCluster(cluster)

#save.image("secrad_sim_220811.Rdata")


sim_list<-sim_list2
rm(sim_list2)
for(i in 1:niter){
	sim_list[[i]]$update_detect()
}


#################check for basic identifiability condition (number of individuals detected at multiple locations)	
sanity<-rep(NA,niter)
for(i in 1:niter){
	sanity[i]<-sum(apply(sim_list[[i]]$obs[[1]]$detect>0,2,sum)>=2)
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

hessian_estim<-lapply(res_list,"[[","hessian")
par_se<-hessian_estim%>%lapply(solve)%>%sapply(diag)%>%sqrt

res_tidy<-tibble(iter=rep(1:niter,each=4),parname=rep(c("dens","con0","con1","g0"),100),true=c(truepar_estim),mle=c(par_estim),se=c(par_se))%>%
		mutate(lci=mle+qnorm(0.025,0,1)*se,uci=mle+qnorm(0.975,0,1)*se,bias=mle-true)
		
res_wide<-res_tidy%>%pivot_wider(id_cols="iter",names_from=parname,values_from=c("true","mle","se","lci","uci","bias"))%>%
		bind_cols(t(truepar))

res_tidy<-res_tidy%>%mutate(parname2=factor(res_tidy$parname,levels=c("con0","con1","dens","g0"),
		labels=c(expression(italic(gamma[0])),'"Effect of landscape"',"'ln(Population density)'","ln(italic('g'[0]))")))
plot_estim<-ggplot(res_tidy%>%filter((parname!="g0")),aes(x=true,y=mle,col=parname))+geom_point(size=1)+
		geom_errorbar(aes(x=true,ymax=uci,ymin=lci),size=0.3)+
		geom_abline(intercept = 0, slope = 1)+
		facet_wrap(~parname2,scales = "free",labeller=label_parsed,ncol=2)+
		theme_cowplot()+
		xlab("True")+ylab("Estimated")+
		theme(legend.position="none")
		
plot_estim;ggsave(paste0("plot_estim_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=16,height=16,units="cm")

RMSE_g0<-sqrt(mean((res_wide$mle_g0-res_wide$true_g0)^2))
MBE_g0<-mean(res_wide$mle_g0-res_wide$true_g0)
################plot

i<-96;sim_list[[i]]$ggsecraddata(covname="X",sample=1);sim_list[[i]]$truemodel
i<-8;j<-37;nobs<-length(sim_list[[i]]$sim_result$detect);nindall<-ncol(sim_list[[i]]$sim_result$detect[[1]]);sim_list[[i]]$ggsecradsim(ind=j);sim_list[[i]]$truemodel



##########NULL model
envmodel_null<-list(D~1,C~1,A~0)
indmodel<-c(A=FALSE,g0=FALSE)
occmodel<-c(A=FALSE,g0=FALSE)
secrad_list_null<-vector("list",niter)
for(i in 1:niter){
	secrad_list_null[[i]]<-secrad$new(secrdata=sim_list[[i]])
	secrad_list_null[[i]]$set_model(envmodel=envmodel_null,indmodel=indmodel,occmodel=occmodel)
}


res_list_null<-vector("list",niter)
initpar<-c(-2,1.5,-3)
for(i in 1:niter){
	cat(i,"\n")
	cat(sim_list[[i]]$out_truepar(),"\n")
	res_list_null[[i]]<-try(optim(initpar,secrad_list_null[[i]]$loglf,method="BFGS",control=list(maxit=1000,trace=2),loglfscale=-1,verbose=F,hessian=T))
	cat(res_list_null[[i]]$par,"\n")
}

par_estim_null<-lapply(res_list_null,"[[","par")
dens_estim_null<-sapply(par_estim_null,"[",1)
dens_bias_null<-dens_estim_null-truepar_estim[1,]

con_estim_null<-sapply(par_estim_null,"[",2)
g0_estim_null<-sapply(par_estim_null,"[",3)


dens_se_null<-lapply(res_list_null,"[[","hessian")%>%
		lapply(solve)%>%
		lapply(diag)%>%
		lapply(sqrt)%>%
		sapply("[",1)

plot(truepar_estim[3,],dens_bias_null)
points(truepar_estim[3,],par_estim[1,]-truepar_estim[1,],col=2)


#####SECR
secrres_list<-vector("list",niter)
for(i in 1:niter){
	trapXY<-data.frame(trapID=1:ntrap,x=trapx,y=trapy,usage=sim_list[[i]]$obs[[1]]$effort)
	traptemp<-read.traps(data=trapXY,detector="count",covnames="usage",binary.usage=F)
	detecttemp<-sim_list[[i]]$obs[[1]]$detect
	detectrc<-which(detecttemp!=0,arr.ind=T)
	triplet<-cbind(detectrc,n=detecttemp[detectrc])
	ntri<-nrow(triplet)
	capttemp<-data.frame(session=numeric(0),ID=numeric(0),occasion=numeric(0),trap=numeric(0))
	for(j in 1:ntri){
		capttemp_j<-data.frame(session=rep(1,triplet[j,"n"]),
						ID=rep(triplet[j,"col"],triplet[j,"n"]),
						occasion=rep(1,triplet[j,"n"]),
						trap=rep(triplet[j,"row"],triplet[j,"n"]))
		capttemp<-rbind(capttemp,capttemp_j)
	}
	captdata<-make.capthist(capttemp,traptemp,fmt="trapID",noccasions=1)
	secrres_list[[i]]<-secr.fit(captdata,model=list(D~1,g0~1,sigma~1),buffer=11,ncores=10)
}

secr_dens<-secrres_list%>%lapply(coef)%>%sapply("[",1,1)-log(10000)
secr_dens_lci<-secrres_list%>%lapply(coef)%>%sapply("[",1,3)-log(10000)
secr_dens_uci<-secrres_list%>%lapply(coef)%>%sapply("[",1,4)-log(10000)

secr_g0<-secrres_list%>%lapply(coef)%>%sapply("[",2,1)
secr_sigma<-secrres_list%>%lapply(coef)%>%sapply("[",3,1)



plotdata<-tibble(true_dens=truepar_estim[1,],secr_dens=secr_dens,secr_dens_lci=secr_dens_lci,secr_dens_uci=secr_dens_uci)
plot_secrdens<-ggplot()+geom_point(data=plotdata,aes(x=true_dens,y=secr_dens),size=1)+
	geom_abline()+
	geom_errorbar(data=plotdata,aes(x=true_dens,ymin=secr_dens_lci,ymax=secr_dens_uci),size=0.2)+
	theme_cowplot()+
	xlab("True ln(density)")+
	ylab("Ln(density) estimated by secr")
	
plot_secrdens;ggsave(paste0("plot_dens_secr_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=9,height=8,units="cm")


###SCR-LCP

source("SCRed.pois.r")

truepar<-sapply(sim_list,function(x) x$out_truepar())

truepar_estim<-truepar
truepar_estim[2,]<-truepar[2,]-truepar[4,]
truepar_estim<-truepar_estim[-4,]
par_estim<-sapply(res_list,"[[","par")

SCRedreslist<-vector("list",niter)

for(i in 1:niter){
	SCRedinit<-c(-3-log(2*pi)-truepar_estim[2,i]+log(nt),-log(2)-truepar_estim[2,i],log(exp(truepar_estim[1,i])*ncell-sim_list[[i]]$nind),0)
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
	

SCRed_code<-SCRedreslist%>%sapply("[[","convergence")
SCRed_con<-(-1)*unlist(SCRedreslist%>%lapply("[[","par")%>%lapply("[",4))
SCRed_dens<-log(exp(unlist(SCRedreslist%>%lapply("[[","par")%>%lapply("[",3)))+sapply(sim_list,"[[","nind"))-log(ncell)


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
					SCRed_n0=exp(SCRed_dens)*ncell-sapply(sim_list,"[[","nind"),true_n0=sim_list%>%lapply(function(x) ncol(x$sim_result$detect[[1]])-x$nind)%>%unlist,secrad_n0=exp(secrad_dens)*ncell-sapply(sim_list,"[[","nind"),
					SCRed_con_lci=qnorm(0.025,SCRed_con,SCRed_con_se),SCRed_con_uci=qnorm(0.975,SCRed_con,SCRed_con_se),
					SCRed_dens_lci=qnorm(0.025,SCRed_dens,SCRed_dens_se),SCRed_dens_uci=qnorm(0.975,SCRed_dens,SCRed_dens_se))
			
plot_SCRLCP_cost<-ggplot()+geom_point(data=SCRed_tibble%>%filter(SCRed_valid),aes(x=true_con,y=SCRed_con),size=1)+
		geom_errorbar(data=SCRed_tibble%>%filter(SCRed_valid),aes(x=true_con,ymax=SCRed_con_uci,ymin=SCRed_con_lci),size=0.2)+
		theme_cowplot()+
		xlab("True effect of landscape\n on connectivity")+
		ylab("(-1)×effect of landscape\n on cost by SCR-LCP")

plot_SCRLCP_cost;ggsave(paste0("plot_SCRLCP_cost_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=9,height=8,units="cm")

plotdata<-SCRed_tibble%>%filter(SCRed_valid)%>%select(true_dens,SCRed_dens,SCRed_dens_lci,SCRed_dens_uci)%>%mutate(model="SCRed")
colnames(plotdata)[2:4]<-c("mle_dens","lci_dens","uci_dens")		
plotdata<-bind_rows(plotdata,res_wide%>%select(true_dens,mle_dens,lci_dens,uci_dens)%>%mutate(model="secrad"))

plot_SCReddens<-ggplot(data=plotdata%>%filter(model=="SCRed"))+geom_point(aes(x=true_dens,y=mle_dens),size=1)+
	geom_abline()+
	geom_errorbar(aes(x=true_dens,ymin=lci_dens,ymax=uci_dens),size=0.2)+
	theme_cowplot()+
	xlab("True ln(density)")+
	ylab("Ln(density) estimated\n by SCR-LCP")
	
plot_SCReddens;ggsave(paste0("plot_dens_SCRLCP_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=9,height=8,units="cm")


#compare distance
costdist<-function(alpha2,cov,G1,G2=NULL,directions=16){
    if(is.null(G2)){
     	G2<-G1
     }
    cost <- exp(alpha2 * cov)                                             
    tr <- transition(cost, transitionFunction=function(x) (1/(mean(x))),  
                     direction = directions)                              
    trLayer <- geoCorrection(tr, scl = F)                                 
    D <- costDistance(trLayer,as.matrix(G1),as.matrix(G2))
    return(D)                  
}

commutedist<-function(alpha2,cov,G1,directions=4){
    cost <- exp(alpha2 * cov)                                             
    tr <- transition(cost, transitionFunction=function(x) (1/exp(mean(log(x)))),  
                     direction = directions)                              
    trLayer <- geoCorrection(tr,  type="r",scl = F)                                 
    D <- commuteDistance(trLayer,as.matrix(G1))
    return(D)                  
}


mantelr<-numeric(niter)
for(i in 1:niter){
	sim_i<-sim_list[[i]]
	ss<-data.frame(X=sim_i$coords[,"x"],Y=sim_i$coords[,"y"],Hab=sim_i$grid_cov$X)
	coordinates(ss)<- ~X+Y
 	gridded(ss) <- TRUE
 	r <- raster(ss, "Hab") 
 	XY<-cbind(trapx,trapy)
 	D_ed<-as.matrix(costdist(SCRed_tibble$SCRed_con[i],r,XY,directions=4))
	D_ad<-as.matrix(commutedist((-1)*res_list[[i]]$par[3],r,XY,directions=4))
	mantelr[i]<-mantel(D_ed,D_ad)$statistic
}


i<-1
	sim_i<-sim_list[[i]]
	ss<-data.frame(X=sim_i$coords[,"x"],Y=sim_i$coords[,"y"],Hab=sim_i$grid_cov$X)
	coordinates(ss)<- ~X+Y
 	gridded(ss) <- TRUE
 	r <- raster(ss, "Hab") 
 	G1<-cbind(trapx,trapy)
  	D_ed<-as.matrix(costdist(SCRed_tibble$SCRed_con[i],r,G1,directions=4))
	D_ad<-as.matrix(commutedist((-1)*res_list[[i]]$par[3],r,G1,directions=4))

	par(mfrow=c(1,2))
	j<-5
	image(r)
	arrows(rep(trapx[j],length(trapx)-1),rep(trapy[1],length(trapy)-1),trapx[-j],trapy[-j],lwd=1/D_ed[j,-j]*min(D_ed[j,-j])*10)
	image(r)
	arrows(rep(trapx[j],length(trapx)-1),rep(trapy[1],length(trapy)-1),trapx[-j],trapy[-j],lwd=1/D_ad[j,-j]*min(D_ad[j,-j])*10)

#compare home range
hrtest<-tibble(MSE_SCRed=numeric(niter),MSE_secrad=numeric(niter),OVER_SCRed=numeric(niter),OVER_secrad=numeric(niter))
ingrid<-(xcoord>=9)&(xcoord<=41)&(ycoord>=9)&(ycoord<=41)
for(i in 1:niter){
	cat(i,"\n")
	true_hr<-matrix(exp(actcent(truepar_estim[,i],secrad_list[[i]],Nhat=FALSE,Utility=FALSE)$logp_arr),ncell,ncell)
	secrad_hr<-matrix(exp(actcent(res_list[[i]]$par,secrad_list[[i]],Nhat=FALSE,Utility=FALSE)$logp_arr),ncell,ncell)
	sim_i<-sim_list[[i]]
	ss<-data.frame(X=sim_i$coords[,"x"],Y=sim_i$coords[,"y"],Hab=sim_i$grid_cov$X)
	coordinates(ss)<- ~X+Y
 	gridded(ss) <- TRUE
 	r <- raster(ss, "Hab") 
 	Gall<-cbind(xcoord,ycoord)
  D_ed_all<-as.matrix(costdist(SCRedreslist[[i]]$par[4],r,Gall,directions=4))
	SCRed_hr<-exp(-exp(SCRedreslist[[i]]$par[2])*D_ed_all^2)
	SCRed_hr_norm<-apply(SCRed_hr,1,sum)
	SCRed_hr<-SCRed_hr/outer(SCRed_hr_norm,rep(1,ncell))

	hrtest$MSE_SCRed[i]<-mean(apply((SCRed_hr[ingrid,]-true_hr[ingrid,])^2,1,sum))
	hrtest$MSE_secrad[i]<-mean(apply((secrad_hr[ingrid,]-true_hr[ingrid,])^2,1,sum))
	hrtest$OVER_SCRed[i]<-mean(apply(pmin(SCRed_hr[ingrid,],true_hr[ingrid,]),1,sum))
	hrtest$OVER_secrad[i]<-mean(apply(pmin(secrad_hr[ingrid,],true_hr[ingrid,]),1,sum))
}

plotdata<-bind_cols(Model=rep(c("SCR-LCP","ADCR"),each=niter),MSE=c(hrtest$MSE_SCRed,hrtest$MSE_secrad),OVER=c(hrtest$OVER_SCRed,hrtest$OVER_secrad),true_con=rep(SCRed_tibble$true_con,2),valid=c(SCRed_tibble$SCRed_valid,rep(TRUE,niter)))
	
plot_overlap_SCRed<-ggplot()+geom_point(data=plotdata%>%filter(valid),aes(x=true_con,y=OVER,col=Model))+
	geom_point(data=plotdata%>%filter(valid)%>%filter(Model=="SCR-LCP"),aes(x=true_con,y=OVER,col=Model),size=2)+
	theme_cowplot()+
	geom_hline(yintercept=1)+
	scale_color_manual(values=c("black","gray"))+
	xlab("True effect of landscape")+ylab("Overlap of predicted\n and true home ranges")

plot_overlap_SCRed;ggsave(paste0("plot_overlap_SCRed_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=12,height=8,units="cm")


plotdata<-bind_cols(model=rep(c("SCRed","secrad"),each=niter),MSE=c(hrtest$MSE_SCRed,hrtest$MSE_secrad),OVER=c(hrtest$OVER_SCRed,hrtest$OVER_secrad),true_con=rep(SCRed_tibble$true_con,2),valid=c(SCRed_tibble$SCRed_valid,rep(TRUE,niter)))
plot_overlap<-ggplot()+geom_point(data=plotdata%>%filter(valid)%>%filter(model=="secrad"),aes(x=true_con,y=OVER),size=1)+
	geom_hline(yintercept=1)+
	theme_cowplot()+xlab("True effect of landscape")+ylab("Overlap of predicted\n and true home range")

plot_overlap;ggsave(paste0("plot_overlap_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=9,height=8,units="cm")


mean((plotdata%>%filter(valid)%>%filter(model=="secrad"))$OVER)
quantile((plotdata%>%filter(valid)%>%filter(model=="secrad"))$OVER,c(0.025,0.975))

i<-60;
	true_hr<-matrix(exp(actcent(truepar_estim[,i],secrad_list[[i]],Nhat=FALSE,Utility=FALSE)$logp_arr),ncell,ncell)
	secrad_hr<-matrix(exp(actcent(res_list[[i]]$par,secrad_list[[i]],Nhat=FALSE,Utility=FALSE)$logp_arr),ncell,ncell)
	sim_i<-sim_list[[i]]
	ss<-data.frame(X=sim_i$coords[,"x"],Y=sim_i$coords[,"y"],Hab=sim_i$grid_cov$X)
	coordinates(ss)<- ~X+Y
 	gridded(ss) <- TRUE
 	r <- raster(ss, "Hab") 
 	Gall<-cbind(xcoord,ycoord)
  	D_ed_all<-as.matrix(costdist(SCRedreslist[[i]]$par[4],r,Gall,directions=16))
	SCRed_hr<-exp(-exp(SCRedreslist[[i]]$par[2])*D_ed_all^2)
	SCRed_hr_norm<-apply(SCRed_hr,1,sum)
	SCRed_hr<-SCRed_hr/outer(SCRed_hr_norm,rep(1,ncell))

j<-1150
	plotdata<-tibble(model=rep(c("true","SCRed","secrad"),each=ncell),x=rep(xcoord,3),y=rep(ycoord,3),prob=c(true_hr[j,],SCRed_hr[j,],secrad_hr[j,]))
	ggplot()+geom_tile(data=plotdata,aes(x=x,y=y,fill=prob))+facet_wrap(~model)+
		geom_point(data=data.frame(sim_list[[i]]$coords,Hab=sim_list[[i]]$grid_cov$X)%>%filter(Hab==0.5),aes(x=x,y=y),size=0.005)+
		geom_point(x=xcoord[j],y=ycoord[j])+
		scale_fill_viridis()+
		theme_cowplot()

j<-1200
	plotdata<-tibble(model=rep(c("True","Predicted"),each=ncell),x=rep(xcoord,2),y=rep(ycoord,2),Probability=c(true_hr[j,],secrad_hr[j,]))
	plot_hr<-ggplot()+geom_tile(data=plotdata,aes(x=x,y=y,fill=Probability))+facet_wrap(~model)+
		geom_point(data=data.frame(sim_list[[i]]$coords)[j,,drop=F],aes(x=x,y=y),col="red",pch=4)+
		geom_point(data=data.frame(sim_list[[i]]$coords,Hab=sim_list[[i]]$grid_cov$X)%>%filter(Hab==0.5),aes(x=x,y=y),size=0.005,stroke=0.4)+
		scale_fill_viridis(begin=0.1)+
		theme_cowplot()
		

plot_hr;ggsave(paste0("plot_hr_",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=17,height=8,units="cm")

j<-1522
	plotdata<-tibble(model=rep(c("True","Predicted by SCR-LCP","Predicted by ADCR"),each=ncell),x=rep(xcoord,3),y=rep(ycoord,3),Probability=c(true_hr[j,],SCRed_hr[j,],secrad_hr[j,]))
	plot_hr_SCRed<-ggplot()+geom_tile(data=plotdata,aes(x=x,y=y,fill=Probability))+facet_wrap(~model)+
		geom_point(data=data.frame(sim_list[[i]]$coords)[j,,drop=F],aes(x=x,y=y),col="red",pch=4)+
		geom_point(data=data.frame(sim_list[[i]]$coords,Hab=sim_list[[i]]$grid_cov$X)%>%filter(Hab==0.5),aes(x=x,y=y),size=0.005,stroke=0.4)+
		scale_fill_viridis(begin=0.1)+
		theme_cowplot()
		
plot_hr_SCRed;ggsave(paste0("plot_hr_SCRLCP",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=23,height=8,units="cm")


###Table 1
rmse_dens_all<-SCRed_tibble%>%
  bind_cols(secr_dens_bias=secr_dens-truepar_estim[1,])%>%
  select(secrad_dens_bias,SCRed_dens_bias,secr_dens_bias)%>%
  summarize_all(.funs=~sqrt(mean(.^2)))

rmse_dens_p<-SCRed_tibble%>%
	bind_cols(secr_dens_bias=secr_dens-truepar_estim[1,])%>%
	filter(true_con>(0.5))%>%
	select(secrad_dens_bias,SCRed_dens_bias,secr_dens_bias)%>%
	summarize_all(.funs=~sqrt(mean(.^2)))

rmse_dens_0<-SCRed_tibble%>%
	bind_cols(secr_dens_bias=secr_dens-truepar_estim[1,])%>%
	filter(abs(true_con)<(0.5))%>%
	select(secrad_dens_bias,SCRed_dens_bias,secr_dens_bias)%>%
	summarize_all(.funs=~sqrt(mean(.^2)))

rmse_dens_n<-SCRed_tibble%>%
	bind_cols(secr_dens_bias=secr_dens-truepar_estim[1,])%>%
	filter(true_con<(-0.5))%>%
	select(secrad_dens_bias,SCRed_dens_bias,secr_dens_bias)%>%
	summarize_all(.funs=~sqrt(mean(.^2)))

bind_rows(rmse_dens_all,rmse_dens_n,rmse_dens_0,rmse_dens_p)

bias_dens_all<-SCRed_tibble%>%
  bind_cols(secr_dens_bias=secr_dens-truepar_estim[1,])%>%
  select(secrad_dens_bias,SCRed_dens_bias,secr_dens_bias)%>%
  summarize_all(.funs=~mean(.))

bias_dens_p<-SCRed_tibble%>%
  bind_cols(secr_dens_bias=secr_dens-truepar_estim[1,])%>%
  filter(true_con>(0.5))%>%
  select(secrad_dens_bias,SCRed_dens_bias,secr_dens_bias)%>%
  summarize_all(.funs=~mean(.))

bias_dens_0<-SCRed_tibble%>%
  bind_cols(secr_dens_bias=secr_dens-truepar_estim[1,])%>%
  filter(abs(true_con)<(0.5))%>%
  select(secrad_dens_bias,SCRed_dens_bias,secr_dens_bias)%>%
  summarize_all(.funs=~mean(.))

bias_dens_n<-SCRed_tibble%>%
  bind_cols(secr_dens_bias=secr_dens-truepar_estim[1,])%>%
  filter(true_con<(-0.5))%>%
  select(secrad_dens_bias,SCRed_dens_bias,secr_dens_bias)%>%
  summarize_all(.funs=~mean(.))

bind_rows(bias_dens_all,bias_dens_n,bias_dens_0,bias_dens_p)


#example of dispersal path
i<-76
	sim_i<-sim_list[[i]]
	ss<-data.frame(X=sim_i$coords[,"x"],Y=sim_i$coords[,"y"],Hab=sim_i$grid_cov$X)
	coordinates(ss)<- ~X+Y
 	gridded(ss) <- TRUE
 	r <- exp(-raster(ss, "Hab")*res_list[[i]]$par[3])
 	XY<-SpatialPoints(cbind(trapx,trapy))
    tr <- transition(r, transitionFunction=function(x) (1/exp(mean(log(x)))),  
                 direction = 4)
    sumpath<-r*0
    for(j in 1:ntrap){
 		for(k in 1:ntrap){
            sumpath<-sumpath+passage(tr,XY[j,],XY[k,])
        }
    }
sumpath2<-sumpath
sumpath2[XY]<-NA

plot_current<-ggplot()+geom_raster(data=sumpath2%>%as.data.frame(xy=TRUE),aes(x=x,y=y,fill=layer))+
		geom_point(data=data.frame(sim_list[[i]]$coords,Hab=sim_list[[i]]$grid_cov$X)%>%filter(Hab==0.5),aes(x=x,y=y),size=0.005)+
		scale_fill_viridis()+
		theme_cowplot()+
		labs(fill = "Current")
plot_current;ggsave(paste0("plot_current",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=10,height=8,units="cm")


###Comparison with 1/2 resolution scenario
reso<-new.env()
load("secrad_sim_SCRed_reso_220812.Rdata",envir=reso)


res_wide_reso<-get("res_wide",pos=reso)
res_wide_reso<-res_wide_reso%>%
					select(mle_con1,se_con1,lci_con1,uci_con1,mle_dens,se_dens,lci_dens,uci_dens)
colnames(res_wide_reso)<-paste0("reso_",colnames(res_wide_reso))
res_wide_all<-res_wide%>%select(mle_con1,lci_con1,uci_con1,mle_dens,se_dens,lci_dens,uci_dens)%>%
				bind_cols(res_wide_reso)

SCRed_tibble_reso<-get("SCRed_tibble",pos=reso)
SCRed_tibble_reso<-SCRed_tibble_reso%>%
					select(SCRed_con,SCRed_con_se,SCRed_con_lci,SCRed_con_uci,SCRed_dens,SCRed_dens_lci,SCRed_dens_uci)
colnames(SCRed_tibble_reso)<-paste0("reso_",colnames(SCRed_tibble_reso))

SCRed_tibble_all<-SCRed_tibble%>%
				select(SCRed_con,SCRed_con_lci,SCRed_con_uci,SCRed_dens,SCRed_dens_lci,SCRed_dens_uci)%>%
				bind_cols(SCRed_tibble_reso)

prinaxis_con_SCRed<-sma(SCRed_con~reso_SCRed_con,data=SCRed_tibble_all,V=cbind(SCRed_con_se^2,reso_SCRed_con_se^2))$coef[[1]]
prinaxis_con_ADCR<-sma(mle_con1~reso_mle_con1,data=res_wide_all,V=cbind(se_con1^2,reso_se_con1^2))$coef[[1]]

plot_half_con_ADCR<-ggplot(res_wide_all,aes(x=reso_mle_con1,y=mle_con1))+
  geom_abline(slope=prinaxis_con_ADCR[2,1],intercept=prinaxis_con_ADCR[1,1],size=1,col="gray")+
	geom_errorbar(aes(x=reso_mle_con1,ymin=lci_con1,ymax=uci_con1),size=0.05)+
	geom_errorbar(aes(y=mle_con1,xmin=reso_lci_con1,xmax=reso_uci_con1),size=0.05)+
	geom_abline(slope=1)+
	geom_point(size=0.5)+
	theme_cowplot()+xlab("effects of landscape (1/2 resolution)")+ylab("effects of landscape (baseline)")
	
plot_half_con_ADCR;ggsave(paste0("plot_half_con_ADCR",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=10,height=8,units="cm")


plot_half_con_SCRed<-ggplot(SCRed_tibble_all,aes(x=reso_SCRed_con,y=SCRed_con))+
  geom_abline(slope=prinaxis_con_SCRed[2,1],intercept=prinaxis_con_SCRed[1,1],size=1,col="gray")+
  geom_errorbar(aes(x=reso_SCRed_con,ymin=SCRed_con_lci,ymax=SCRed_con_uci),size=0.05)+
	geom_errorbar(aes(y=SCRed_con,xmin=reso_SCRed_con_lci,xmax=reso_SCRed_con_uci),size=0.05)+
	geom_abline(slope=1)+
  geom_point(size=0.5)+
	theme_cowplot()+xlab("effects of landscape (1/2 resolution)")+ylab("effects of landscape (baseline)")

plot_half_con_SCRed;ggsave(paste0("plot_half_con_SCRed",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=10,height=8,units="cm")



prinaxis_dens_SCRed<-sma(SCRed_dens~reso_SCRed_dens,data=SCRed_tibble_all,V=cbind(SCRed_dens_se^2,reso_SCRed_dens_se^2))$coef[[1]]
prinaxis_dens_ADCR<-sma(mle_dens~reso_mle_dens,data=res_wide_all,V=cbind(se_dens^2,reso_se_dens^2))$coef[[1]]

plot_half_dens_ADCR<-ggplot(res_wide_all,aes(x=reso_mle_dens,y=mle_dens))+
  geom_abline(slope=prinaxis_dens_ADCR[2,1],intercept=prinaxis_dens_ADCR[1,1],size=1,col="gray")+
  geom_errorbar(aes(x=reso_mle_dens,ymin=lci_dens,ymax=uci_dens),size=0.05)+
  geom_errorbar(aes(y=mle_dens,xmin=reso_lci_dens,xmax=reso_uci_dens),size=0.05)+
  geom_abline(slope=1)+
  geom_point(size=0.5)+
  theme_cowplot()+xlab("log popul. density (1/2 resolution)")+ylab("log popul. density (baseline)")

plot_half_dens_ADCR;ggsave(paste0("plot_half_dens_ADCR",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=10,height=8,units="cm")


plot_half_dens_SCRed<-ggplot(SCRed_tibble_all,aes(x=reso_SCRed_dens,y=SCRed_dens))+
  geom_abline(slope=prinaxis_dens_SCRed[2,1],intercept=prinaxis_dens_SCRed[1,1],size=1,col="gray")+
  geom_errorbar(aes(x=reso_SCRed_dens,ymin=SCRed_dens_lci,ymax=SCRed_dens_uci),size=0.05)+
  geom_errorbar(aes(y=SCRed_dens,xmin=reso_SCRed_dens_lci,xmax=reso_SCRed_dens_uci),size=0.05)+
  geom_abline(slope=1)+
  geom_point(size=0.5)+
  theme_cowplot()+xlab("log popul. density (1/2 resolution)")+ylab("log popul. density (baseline)")

plot_half_dens_SCRed;ggsave(paste0("plot_half_dens_SCRed",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=10,height=8,units="cm")

