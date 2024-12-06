library(sf)
library(tidyverse)
library(cowplot)
library(lubridate)
library(raster)
library(gdistance)
library(viridis)
library(secr)

sourcepath<-"secrad.r"
source(sourcepath)

effort<-read_csv("effort_231225.csv")
effort_st<-effort%>%st_as_sf(coords=c("x","y"),crs=3100)

detect<-read_csv("detectmat_231225.csv")

#read griddata
gridpath<-getwd()
ngrid<-1
dataset<-list()
dataset$griddata<-st_read(gridpath,"meshutm_0.5km_buff_land")

#dataset
coords<-dataset$griddata[,c("x","y")]/1000	#km
coords<-coords%>%as_tibble%>%select(-geometry)
dataset$coords<-coords%>%as.matrix

ncell<-nrow(coords)
ncell<-nrow(coords)
xdist<-abs(outer(coords$x,rep(1,ncell))-outer(rep(1,ncell),coords$x))
ydist<-abs(outer(coords$y,rep(1,ncell))-outer(rep(1,ncell),coords$y))
dx<-min(xdist[xdist!=0])
dy<-min(ydist[ydist!=0])

area<-dataset$griddata$area/1000/1000
dataset$area<-area

grid_cov<-dataset$griddata%>%as_tibble()%>%
  mutate(agri=agri_mean,wtr=wtr_mean)%>%
  select(agri,wtr,-geometry)
mu_agri<-mean(grid_cov$agri)
sd_agri<-sd(grid_cov$agri)
mu_wtr<-mean(grid_cov$wtr)
sd_wtr<-sd(grid_cov$wtr)
grid_cov_std<-grid_cov%>%mutate(agri=(agri-mu_agri)/sd_agri,agri=(agri-mu_agri)/sd_agri)
dataset$grid_cov_std<-grid_cov_std
dataset$grid_cov_musd<-data.frame(mu=c(mu_agri,mu_wtr),sd=c(sd_agri,sd_wtr))
resolution<-c(x=dx,y=dy)

dataset$resolution<-resolution

dataset$effort<-effort$effort
effort_loc<-st_intersects(effort_st,dataset$griddata)%>%
  unlist
dataset$effort_loc<-effort_loc

dataset$detect<-as.matrix(detect[,-1])
dataset$effort_occ<-effort$effort_occ

secrdata<-secrad_data$new(coords=dataset$coords,
                          area=dataset$area,
                          grid_cov=dataset$grid_cov_std,
                          resolution=dataset$resolution)
secrdata$add_obs(type="poisson",
                 effort=dataset$effort,
                 effort_loc=dataset$effort_loc,
                 effort_occ=dataset$effort_occ,
                 detect=dataset$detect)

secrdata$ggsecraddata(covname="agri")

secrdata$ggsecraddata(covname="wtr")

secrad_obj<-secrad$new(secrdata=secrdata)
envmodel<-list(D~1,C~agri+wtr,A~0)
indmodel<-c(A=FALSE,g0=FALSE)
occmodel<-c(A=FALSE,g0=FALSE)
secrad_obj$set_model(envmodel=envmodel,indmodel=indmodel,occmodel=occmodel)

initpar<-generate_init(secrad_obj)
initpar["dens_0"]<--1
initpar["conn_0"]<--2
initpar["g0_1"]<--5
secrad_res<-optim(initpar,secrad_obj$loglf,method="BFGS",control=list(maxit=1000,trace=2),loglfscale=-1,verbose=T,hessian=T)

#save.image("secrad_toyama_agri240412.Rdata")

####SCRed
source("SCRed.pois.r")

SCRedinit<-c(-3-log(2*pi)-(-5),-log(2)-(-2),log(exp(-1)*sum(secrdata$area)-secrdata$nind),0.1,0.1)
sim_i<-secrdata
tdf1<-data.frame(De=1:length(secrdata$obs[[1]]$effort),
                 X=sim_i$coords[sim_i$obs[[1]]$effort_loc,"x"],
                 Y=sim_i$coords[sim_i$obs[[1]]$effort_loc,"y"],
                 Cov=sim_i$obs[[1]]$effort)

"[.list2"<-function(obj,x){obj[[x]]}	#for hacking "[" call


SCRed.pois.2v<-function(par,sim_i,tdf1,g0offset=NULL,...){
  cat(par,"\n")
  coef1<-par[length(par)-1]
  coef2<-par[length(par)]
  cov<-unlist(sim_i$grid_cov[,1]+coef2/coef1*sim_i$grid_cov[,2])
  ss<-data.frame(X=sim_i$coords[,"x"],Y=sim_i$coords[,"y"],Hab=cov)
  coordinates(ss)<- ~X+Y
  gridded(ss) <- TRUE
  r <- raster(ss, "Hab") 
  XY<-as.matrix(coordinates(r))
  par2<-par[-length(par)]
  if(!is.null(g0offset)){
    par2<-as.list(par2)
    par2[[1]]<-par2[[1]]+g0offset
    class(par2)<-"list2"
  }
  SCRedLL<-SCRed.pois(par2,y=t(sim_i$obs[[1]]$detect),X=tdf1[,c("X","Y")],cov=r,G=XY,...)
  cat(SCRedLL,"\n")
  return(SCRedLL)
}


SCRedres<-optim(SCRedinit,SCRed.pois.2v,sim_i=sim_i,tdf1=tdf1,ssbuffer=10,directions=4,g0offset=log(tdf1$Cov),hessian=T,control=list(maxit=1000))


#SECR
trapXY<-data.frame(trapID=1:length(secrdata$obs[[1]]$effort),x=secrdata$obs[[1]]$effort_coords[,"x"]*1000,y=secrdata$obs[[1]]$effort_coords[,"y"]*1000,usage=log(secrdata$obs[[1]]$effort))
traptemp<-read.traps(data=trapXY,detector="count",covnames="usage",binary.usage=F)
detecttemp<-secrdata$obs[[1]]$detect
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
secrres<-secr.fit(captdata,model=list(D~1,g0~1,sigma~1),start=c(0,-5,8),buffer=10000,ncores=10)

logdens_secr<-summary(secrres)$coef[1,]
logdens_secr[,c(1,3,4)]<-logdens_secr[,c(1,3,4)]+log(1000)*2-log(10000)	#log density in km2
logdens_secr

derived(secrres)[2,]*1000*1000/10000 #in km2


#save.image("secrad_toyama_agri241206.Rdata")




####

parname<-c("log(density)","gamma_0","beta_agri","beta_water","g0")
mle<-secrad_res$par
se<-sqrt(diag(solve(secrad_res$hessian)))
estimate<-tibble(name=parname,mle=mle,se=se,ci2.5=qnorm(0.025,mle,se),ci97.5=qnorm(0.975,mle,se))

parname<-c("g0","sigma","log(density)","beta_agri","beta_water")
mle<-SCRedres$par
mle[3]<-log((exp(mle[3])+secrdata$nind))-log(prod(apply(apply(tdf1[,c("X","Y")],2,range),2,diff)+20))
se<-sqrt(diag(solve(SCRedres$hessian)))
se[3]<-se[3]*exp(SCRedres$par[3])/(exp(SCRedres$par[3])+secrdata$nind)		#se of log density with delta method
estimateSCRed<-tibble(name=parname,mle=mle,se=se,ci2.5=qnorm(0.025,mle,se),ci97.5=qnorm(0.975,mle,se))


beta_all<-bind_rows(bind_cols(model="ADCR",estimate%>%filter(grepl("^beta",name))),
                    bind_cols(model="SCR_LCP",estimateSCRed%>%filter(grepl("^beta",name)))
)%>%bind_cols(Variable=rep(c("Agriculture","Water"),2))

beta<-beta%>%filter(model=="ADCR")

labs<-c(beta_agri="Agriculture",beta_water="Water")
toyama_coef_ADCR<-ggplot(beta)+geom_point(aes(x=Variable,y=mle),size=4)+
  geom_errorbar(aes(x=Variable,ymin=ci2.5,ymax=ci97.5),size=1)+
  geom_hline(yintercept=0)+
  theme_cowplot()+
  ylab("Effect of landscape on \n permeability by ADCR")

toyama_coef_ADCR;ggsave(paste0("plot_coef_toyama_ADCR",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=8,height=12,units="cm")

beta_SCR_LCP<-beta_all%>%filter(model=="SCR_LCP")

toyama_coef_SCR_LCP<-ggplot(beta_SCR_LCP)+geom_point(aes(x=Variable,y=mle),size=4)+
  geom_errorbar(aes(x=Variable,ymin=ci2.5,ymax=ci97.5),size=1)+
  geom_hline(yintercept=0)+
  theme_cowplot()+
  ylab("Effect of landscape on \n cost by SCR-LCP")

toyama_coef_SCR_LCP;ggsave(paste0("plot_coef_toyama_SCRLCP",format(Sys.time(), "%Y%m%d%H%M"),".pdf"),device="pdf",width=8,height=12,units="cm")

#predicted home range
secrad_hr<-t(matrix(exp(actcent(secrad_res$par,secrad_obj,Nhat=FALSE,Utility=FALSE)$logp_arr),ncell,ncell))
colnames(secrad_hr)<-paste0("X",1:ncell)

plotdata<-tibble(x=coords$x,y=coords$y,as.data.frame(secrad_hr))

ggplot()+geom_tile(data=plotdata,aes(x=x,y=y,fill=X127))+
  geom_point(data=data.frame(secrdata$coords,Hab=secrdata$grid_cov$wtr)%>%filter(Hab>quantile(Hab,0.75)),aes(x=x,y=y),size=0.005)+
  geom_point(data=plotdata[127,],aes(x=x,y=y),color="red",size=2,pch=4)+
  xlim(coords$x[127]+c(-10,10))+ylim(coords$y[128]+c(-10,10))+
  scale_fill_viridis()+
  theme_cowplot()

