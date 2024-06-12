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


parname<-c("log(density)","gamma_0","beta_agri","beta_water","g0")
mle<-secrad_res$par
se<-sqrt(diag(solve(secrad_res$hessian)))
estimate<-tibble(name=parname,mle=mle,se=se,ci2.5=qnorm(0.025,mle,se),ci97.5=qnorm(0.975,mle,se))


beta<-bind_cols(model="ADCR",Variable=c("Agriculture","Water"),estimate%>%filter(grepl("^beta",name)))

ggplot(beta)+geom_point(aes(x=Variable,y=mle),size=4)+
	geom_errorbar(aes(x=Variable,ymin=ci2.5,ymax=ci97.5),size=1)+
	geom_hline(yintercept=0)+
	theme_cowplot()+
	ylab("Coefficient")


#predicted home range
secrad_hr<-t(matrix(exp(actcent(secrad_res$par,secrad_obj,Nhat=FALSE,Utility=FALSE)$logp_arr),ncell,ncell))
colnames(secrad_hr)<-paste0("X",1:ncell)

plotdata<-tibble(x=coords$x,y=coords$y,as.data.frame(secrad_hr))

ggplot()+geom_tile(data=plotdata,aes(x=x,y=y,fill=X127))+
	geom_point(data=data.frame(secrdata$coords,Hab=secrdata$grid_cov$wtr)%>%filter(Hab!=0),aes(x=x,y=y),size=0.005)+
	geom_point(data=plotdata[127,],aes(x=x,y=y),color="red",size=2,pch=4)+
	xlim(coords$x[127]+c(-10,10))+ylim(coords$y[128]+c(-10,10))+
	scale_fill_viridis(name="Probability",begin=0.1)+
	theme_cowplot()


