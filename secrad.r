########################################
###Spatially-explicit capture-recapture model with advection-diffusion home range
########################################

library(sf)
library(tidyverse)
library(units)
library(cowplot)
library(Matrix)
library(inline)
library(RcppEigen)
library(Rcpp)
library(RcppNumerical)
library(R6)

eigenapprox.code<-'
#include <iostream>
#include <RcppEigen.h>
#include <omp.h>
#include <random>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;
using namespace Eigen;
using namespace std;


//Solve linear system using LU decomposition from SparseMatrix and VectorXd (fast)

// [[Rcpp::export]]
VectorXd eigenapprox_solve(const SparseMatrix<double> mat2,const VectorXd b){
    Eigen::VectorXd x;  // out

    Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;  // solver
    solver.compute(mat2);
    if( solver.info() != Eigen::Success ) {
      stop("LU decomposition failed");        
    }else{
        x = solver.solve(b);
        if( solver.info() != Eigen::Success ) {
           stop("solving failed");
        }
    }
    
  return x;
}

//Solve linear system using QR decomposition from SparseMatrix and VectorXd (slow but stable)

// [[Rcpp::export]]
VectorXd eigenapprox_solveQR(const SparseMatrix<double> mat2,const VectorXd b){
    Eigen::VectorXd x;  // out

    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;  // solver
    solver.compute(mat2);

    if( solver.info() != Eigen::Success ) {
      stop("QR decomposition failed");        
    }else{
        x = solver.solve(b);
        if( solver.info() != Eigen::Success ) {
           stop("solving failed");
        }
    }
    
  return x;
}

//Make diagonal SparseMatrix from VectorXd

SparseMatrix<double> vec2diagSp(const VectorXd vec){
    int n=vec.size();
    int i;
    vector<Triplet<double> > t(n);
    for(i=0;i<n;i++){
        t[i]=Triplet<double>(i,i,vec(i));
    }
    SparseMatrix<double> sp1(n,n);
    sp1.setFromTriplets(t.begin(),t.end());
    return(sp1);
}


//Vector of non-zero diagonal element of SparseMatrix

VectorXd sp2diagvec(const SparseMatrix<double> sp1){
    int n=sp1.rows();
    int i;
    VectorXd res=VectorXd::Zero(n);
    for(i=0;i<sp1.outerSize();++i){
        for(SparseMatrix<double>::InnerIterator it(sp1,i); it; ++it){
            if(it.row()==it.col()){
                res[it.row()]=it.value();
            }
        }
    }
    return(res);
}

//logsumexp()
//[[Rcpp::export]]
double sum_lp(VectorXd lp){
	double max_lp = lp.maxCoeff(); 
	return log((lp.array()-max_lp).exp().sum())+max_lp;
}

//[[Rcpp::export]]
VectorXd rowsum_lp(MatrixXd lp){
	int nrow = lp.rows();
	int i;
	double temp;
	Eigen::setNbThreads(1);
	VectorXd res(nrow);
	#pragma omp parallel for private(i,temp)
	for(i=0;i<nrow;i++){
		temp = sum_lp(lp.row(i));
		#pragma omp critical
		res(i) = temp;
    	}
	Eigen::setNbThreads(0);
		
	return res;
}

//[[Rcpp::export]]
VectorXd colsum_lp(MatrixXd lp){
	int ncol = lp.cols();
	int i;
	double temp;
	Eigen::setNbThreads(1);
	VectorXd res(ncol);
	#pragma omp parallel for private(i,temp)
	for(i=0;i<ncol;i++){
		temp = sum_lp(lp.col(i));
		#pragma omp critical
		res(i) = temp;
    	}
	Eigen::setNbThreads(0);
		
	return res;
}


//log combination
double lcomb(const int N,
		const int k){
	return(lgamma(N+1)-lgamma(N-k+1)-lgamma(k+1));
}

//[[Rcpp::export]]
double log1m_exp(double a) {
  using std::log;
  using std::exp;
  if (a > 0)
    return std::numeric_limits<double>::quiet_NaN();
  else if(a==0)
     return R_NegInf;
  else if (a > -0.693147)
    return log(-expm1(a));  // 0.693147 ~= log(2)
  else
    return log1p(-exp(a));
}

//[[Rcpp::export]]
MatrixXd log1m_exp_mat(const MatrixXd a) {
	int nrow = a.rows();
	int ncol = a.cols();
	int i,j;
	double temp;
	MatrixXd res(nrow,ncol);
	Eigen::setNbThreads(1);
	#pragma omp parallel for private(i,j,temp)
	for(i=0;i<nrow;i++){
		for(j=0;j<ncol;j++){
		 	temp = log1m_exp(a(i,j));
			#pragma omp critical
    			res(i,j) = temp;
    		}
    	}
	Eigen::setNbThreads(0);
	return(res);    	
}

//[[Rcpp::export]]
double dbin_c(const int k,
			const int N,
			const double lp,
			const bool logprob=true){
	double res;
	if(((k==0)&(exp(lp)==0))|((k==N)&(exp(lp)==1))){
		res = 0;
	}else if(((k!=0)&(exp(lp)==0))|((k!=N)&(exp(lp)==1))){
		res = R_NegInf;
	}else{
		res = lcomb(N,k)+lp*k+log1m_exp(lp)*(N-k);
	}
	if(logprob==false){
		res = exp(res);
	}
	return(res);
}

//[[Rcpp::export]]
double dpoisson_c(const int k,
			const double loglambda,
			const bool logprob=true){
	double res;
	if((k==0)&(exp(loglambda)==0)){
		res = 0;
	}else if((k!=0)&(exp(loglambda)==0)){
		res = R_NegInf;
	}else{
		res = loglambda*k-exp(loglambda)-lgamma(k+1);
	}
	if(logprob==false){
		res = exp(res);
	}
	return(res);
}


//matrix extended binomial prob

//[[Rcpp::export]]
MatrixXd dbin_mat(const MatrixXi k,
				const MatrixXi N,
				const MatrixXd lp,
				const bool logprob=true){
	int i,j;
	int nmu=lp.rows();
	int neffort=lp.cols();
	MatrixXd res(nmu,neffort);
	Eigen::setNbThreads(1);
	#pragma omp parallel for private(i,j)
	for(i=0;i<nmu;i++){
		for(j=0;j<neffort;j++){
			#pragma omp critical
			res(i,j)=dbin_c(k(i,j),N(i,j),lp(i,j),logprob);
		}
	}
	Eigen::setNbThreads(0);
	return(res);
}

//matrix extended binomial prob

//[[Rcpp::export]]
MatrixXd dpoisson_mat(const MatrixXi k,
				const MatrixXd loglambda,
				const bool logprob=true){
	int i,j;
	int nmu=loglambda.rows();
	int neffort=loglambda.cols();
	MatrixXd res(nmu,neffort);
	Eigen::setNbThreads(1);
	#pragma omp parallel for private(i,j)
	for(i=0;i<nmu;i++){
		for(j=0;j<neffort;j++){
			#pragma omp critical
			res(i,j) = dpoisson_c(k(i,j),loglambda(i,j),logprob);
		}
	}
	Eigen::setNbThreads(0);
	return(res);
}

//[[Rcpp::export]]
MatrixXd pcap_bin(const MatrixXi detect,
				const VectorXi effort,
				const VectorXi effort_occ,
				const MatrixXd lp_mat,
				const MatrixXi srv,
				const VectorXi ind_cov,
				const bool logprob=true){
	int i,j,k,colnum;
	double temp;
	int nmu=lp_mat.rows();
	int neffort=effort_occ.size();
	int nind=ind_cov.size();
	MatrixXd res = MatrixXd::Zero(nmu,nind);
	Eigen::setNbThreads(1);
	#pragma omp parallel for private(i,j,k,colnum,temp)
	for(i=0;i<nmu;i++){
		for(j=0;j<nind;j++){
			for(k=0;k<neffort;k++){
				colnum = k+neffort*(ind_cov(j)-ind_cov.minCoeff());
				temp = dbin_c(detect(k,j),effort(k),lp_mat(i,colnum)+log(srv(j,effort_occ(k)-1)),logprob);
				#pragma omp atomic
				res(i,j) += temp;
			}
		}
	}
	Eigen::setNbThreads(0);
	return(res);
}

//[[Rcpp::export]]
MatrixXd pcap_poisson(const MatrixXi detect,
				const VectorXi effort_occ,
				const MatrixXd loglambda_mat,
				const MatrixXi srv,
				const VectorXi ind_cov,
				const bool logprob=true){
	int i,j,k,colnum;
	double temp;
	int nmu=loglambda_mat.rows();
	int neffort=effort_occ.size();
	int nind=ind_cov.size();
	MatrixXd res= MatrixXd::Zero(nmu,nind);
	Eigen::setNbThreads(1);
	#pragma omp parallel for private(i,j,k,colnum,temp)
	for(i=0;i<nmu;i++){
		for(j=0;j<nind;j++){
			for(k=0;k<neffort;k++){
				colnum = k+neffort*(ind_cov(j)-ind_cov.minCoeff());
				temp = dpoisson_c(detect(k,j),loglambda_mat(i,colnum)+log(srv(j,effort_occ(k)-1)),logprob);
				#pragma omp atomic
				res(i,j) += temp;
			}
		}
	}
	Eigen::setNbThreads(0);
	return(res);
}

//[[Rcpp::export]]
MatrixXd pcap_poisson_debug(const MatrixXi detect,
				const VectorXi effort_occ,
				const MatrixXd loglambda_mat,
				const MatrixXi srv,
				const VectorXi ind_cov,
				const bool logprob=true){
	int i,j,k,colnum,rescol;
	double temp;
	int nmu=loglambda_mat.rows();
	int neffort=effort_occ.size();
	int nind=ind_cov.size();
	MatrixXd res= MatrixXd::Zero(nmu,nind*neffort);
	Eigen::setNbThreads(1);
	#pragma omp parallel for private(i,j,k,colnum,rescol,temp)
	for(i=0;i<nmu;i++){
		for(j=0;j<nind;j++){
			for(k=0;k<neffort;k++){
				colnum = k+neffort*(ind_cov(j)-ind_cov.minCoeff());
				rescol = k+neffort*j;
				temp = dpoisson_c(detect(k,j),loglambda_mat(i,colnum)+log(srv(j,effort_occ(k)-1)),logprob);
				#pragma omp critical
				res(i,rescol) = temp;
			}
		}
	}
	Eigen::setNbThreads(0);
	return(res);
}

//solutions of advection reaction model for all the potential activity center

// [[Rcpp::export]]
MatrixXd advdiff_core(const SparseMatrix<double> neighmatx,
                        const SparseMatrix<double> neighmaty,
                        const SparseMatrix<double> temp0,
                        const SparseMatrix<double> sumtemp0,
                        const VectorXd xcoord,
                        const VectorXd ycoord,
                        const VectorXd mux,
                        const VectorXd muy,
                        const VectorXd loga,
                        const double dx,
                        const double dy,
                        double delta=0.001,
                        bool logprob=true){
    int i,j;
    int nmu=mux.size();
    int ncoord=xcoord.size();
    MatrixXd res(nmu,ncoord);
    //VectorXd res;

    Eigen::setNbThreads(1);
    SparseMatrix<double, ColMajor> advneighx, advneighy, difdiagx, difdiagy,temp1,advdiff;
    VectorXd sumtemp1,eigenvec;
    VectorXd Ones_ncoord=VectorXd::Ones(ncoord);
    bool cond;
    double denom;

    vector<Triplet<double> > t_mass(ncoord);
    SparseMatrix<double,ColMajor> mass(ncoord,ncoord);
    for(j=0;j<ncoord;j++){
        t_mass[j]=Triplet<double>(ncoord-1,j,1.0);
    }
    mass.setFromTriplets(t_mass.begin(),t_mass.end());

    VectorXd b=VectorXd::Zero(ncoord);
    b[ncoord-1]=1;
    std::vector<int> neighmatx_row(neighmatx.nonZeros());
    std::vector<int> neighmatx_col(neighmatx.nonZeros());
    int k=0;
    for(j=0;j<neighmatx.outerSize();j++){
    	for(SparseMatrix<double>::InnerIterator it(neighmatx,j); it; ++it){
          neighmatx_row[k] = it.row();
          neighmatx_col[k] = it.col();
          k++;
         }
      }
      
    std::vector<int> neighmaty_row(neighmaty.nonZeros());
    std::vector<int> neighmaty_col(neighmaty.nonZeros());
    k=0;
    for(j=0;j<neighmaty.outerSize();j++){
    	for(SparseMatrix<double>::InnerIterator it(neighmaty,j); it; ++it){
          neighmaty_row[k] = it.row();
          neighmaty_col[k] = it.col();
          k++;
         }
      }
    
      
    #pragma omp parallel for private(i,j,advneighx,advneighy,cond, difdiagx, difdiagy,temp1,sumtemp1,advdiff,denom,eigenvec)
    for(i=0;i<nmu;i++){
        //creating advection matrix for x
        advneighx=neighmatx;
        for(j=0;j<neighmatx.nonZeros();++j){
            cond=(neighmatx.coeff(neighmatx_row[j],neighmatx_col[j])==1) && (((xcoord[neighmatx_row[j]]<xcoord[neighmatx_col[j]]) && (xcoord[neighmatx_row[j]]<mux[i])) || ((xcoord[neighmatx_row[j]]>xcoord[neighmatx_col[j]]) && (xcoord[neighmatx_row[j]]>mux[i])));
            if(cond){
               advneighx.coeffRef(neighmatx_row[j],neighmatx_col[j])=0;
               }
            }
        
        //creating advection matrix for y
        advneighy=neighmaty;
        for(j=0;j<neighmaty.nonZeros();++j){
            cond=(neighmaty.coeff(neighmaty_row[j],neighmaty_col[j])==1) && (((ycoord[neighmaty_row[j]]<ycoord[neighmaty_col[j]]) && (ycoord[neighmaty_row[j]]<muy[i])) || ((ycoord[neighmaty_row[j]]>ycoord[neighmaty_col[j]]) && (ycoord[neighmaty_row[j]]>muy[i])));
            if(cond){
                advneighy.coeffRef(neighmaty_row[j],neighmaty_col[j])=0;
                }
            }
        
        
        //advection matrix (upwind)
        difdiagx=vec2diagSp((xcoord.array()-mux[i]).abs().matrix());
        difdiagy=vec2diagSp((ycoord.array()-muy[i]).abs().matrix());
        temp1=(advneighx*difdiagx/dx+advneighy*difdiagy/dy)*exp(loga[i]);
        sumtemp1=(temp1.transpose())*Ones_ncoord;
        temp1=temp1-vec2diagSp(sumtemp1);

        //advection-diffusion linear system
        advdiff=temp0-sumtemp0+temp1;
        denom=1/abs(sp2diagvec(advdiff).minCoeff())*delta;
        advdiff=advdiff*denom;
        advdiff.row(ncoord-1)*=0;
        advdiff=advdiff+mass;

        try
        {
            eigenvec=eigenapprox_solve(advdiff,b);
        }
        catch(...)
        {
            try
            {
                eigenvec=eigenapprox_solveQR(advdiff,b);
            }
            catch(...)
            {
                eigenvec=VectorXd::Constant(ncoord,NAN);
            }
        }
        eigenvec=eigenvec/eigenvec.sum();
        for(j=0;j<ncoord;j++){
            if(eigenvec[j]<0){
                eigenvec[j]=0;
            }
        }
        eigenvec=(eigenvec.array()).log().matrix();
        eigenvec=(eigenvec.array()-sum_lp(eigenvec)).matrix();

        if(logprob==false){
            eigenvec=eigenvec.array().exp().matrix();
        }
        #pragma omp critical
        res.block(i,0,1,ncoord)=eigenvec.transpose();

        //if(i==(nmu-1)){
        //   res=eigenvec;
        //}
    }
    
    Eigen::setNbThreads(0);

    return(res);
}

// [[Rcpp::export]]
Rcpp::List advdiff_t_core(const VectorXd T,
			const VectorXd p0,
			const double mux,
			const double muy,
			const SparseMatrix<double> neighmatx,
			const SparseMatrix<double> neighmaty,
			const SparseMatrix<double> temp0,
			const SparseMatrix<double> sumtemp0,
			const VectorXd xcoord,
			const VectorXd ycoord,
			const double loga,
			const double dx,
			const double dy,
			const double nstep=10){
			
	std::random_device seed_gen;
	std::mt19937 engine(seed_gen());
	int i,j,k;
	int niter = T.size();
	int ncoord=p0.size();
	Rcpp::List res(2);
	MatrixXd resprob(ncoord,niter);
	VectorXi resloc(niter);
	VectorXd pt(ncoord);
	double dt;
	SparseMatrix<double, ColMajor> advneighx, advneighy, difdiagx, difdiagy,temp1,advdiff;
	VectorXd sumtemp1;
	VectorXd Ones_ncoord=VectorXd::Ones(ncoord);
	SparseMatrix<double, ColMajor> Imat =vec2diagSp(Ones_ncoord) ;
	bool cond;

	std::vector<int> neighmatx_row(neighmatx.nonZeros());
	std::vector<int> neighmatx_col(neighmatx.nonZeros());
	k=0;
	for(j=0;j<neighmatx.outerSize();j++){
		for(SparseMatrix<double>::InnerIterator it(neighmatx,j); it; ++it){
			neighmatx_row[k] = it.row();
			neighmatx_col[k] = it.col();
			k++;
		}
	}
	
	std::vector<int> neighmaty_row(neighmaty.nonZeros());
	std::vector<int> neighmaty_col(neighmaty.nonZeros());
	k=0;
	for(j=0;j<neighmaty.outerSize();j++){
		for(SparseMatrix<double>::InnerIterator it(neighmaty,j); it; ++it){
			neighmaty_row[k] = it.row();
			neighmaty_col[k] = it.col();
			k++;
		}
	}
	
	
	//creating advection matrix for x
	advneighx=neighmatx;
	for(j=0;j<neighmatx.nonZeros();++j){
		cond=(neighmatx.coeff(neighmatx_row[j],neighmatx_col[j])==1) && (((xcoord[neighmatx_row[j]]<xcoord[neighmatx_col[j]]) && (xcoord[neighmatx_row[j]]<mux)) || ((xcoord[neighmatx_row[j]]>xcoord[neighmatx_col[j]]) && (xcoord[neighmatx_row[j]]>mux)));
		if(cond){
			advneighx.coeffRef(neighmatx_row[j],neighmatx_col[j])=0;
		}
	}

	//creating advection matrix for y
	advneighy=neighmaty;
	for(j=0;j<neighmaty.nonZeros();++j){
		cond=(neighmaty.coeff(neighmaty_row[j],neighmaty_col[j])==1) && (((ycoord[neighmaty_row[j]]<ycoord[neighmaty_col[j]]) && (ycoord[neighmaty_row[j]]<muy)) || ((ycoord[neighmaty_row[j]]>ycoord[neighmaty_col[j]]) && (ycoord[neighmaty_row[j]]>muy)));
		if(cond){
			advneighy.coeffRef(neighmaty_row[j],neighmaty_col[j])=0;
		}
	}
	
	//advection matrix
	difdiagx=vec2diagSp((xcoord.array()-mux).abs().matrix());
	difdiagy=vec2diagSp((ycoord.array()-muy).abs().matrix());
	temp1=(advneighx*difdiagx/dx+advneighy*difdiagy/dy)*exp(loga);
	sumtemp1=(temp1.transpose())*Ones_ncoord;
	temp1=temp1-vec2diagSp(sumtemp1);

	//advection-diffusion linear system
	std::vector<double> mprob(ncoord,1);
	vector<std::discrete_distribution<std::size_t>> d={};
	pt=p0;
	dt=1/nstep;
	advdiff=Imat+(temp0-sumtemp0+temp1)*dt;

	for(j=0;j<niter;j++){
	      for(i=0;i<T[j]*nstep;i++){
			pt = advdiff*pt;
		}
		for(k=0;k<ncoord;k++){
			mprob[k]=pt(k);
		}
		d.push_back(discrete_distribution<std::size_t> (mprob.begin(),mprob.end()));
		resloc(j) = (int)d[j](engine);
		resprob.col(j) = pt;
		pt=VectorXd::Zero(ncoord);
		pt(resloc(j)) = 1;
	}
	res(0) = (resloc.array()+1).matrix();
	res(1) =resprob;
	return(res);
}

'        
                    
sourceCpp(code=eigenapprox.code)        


omitzero.secrdata<-function(secrdata){
	detect_append<-secrdata$obs[[1]]$detect
	nobs<-length(secrdata$obs)
	if(nobs>1){
		for(i in 2:nobs){
			detect_append<-rbind(detect_append,secrdata$obs[[i]]$detect)
		}
	}
	cond<-colSums(detect_append)!=0
	for(i in 1:nobs){
		secrdata$obs[[i]]$detect<-secrdata$obs[[i]]$detect[,cond]
	}
	secrdata$nind<-sum(cond)
	if(!is.null(secrdata$ind_cov)){
		secrdata$ind_cov<-secrdata$ind_cov[cond]
	}
	return(secrdata)
}


advdiff.eigen<-function(mucoords,gridcoords,logc,loga=NULL,neighmat,resolution,logprob=T){
	mux<-mucoords[,1]
	muy<-mucoords[,2]
	
	xcoord<-gridcoords[,1]
	ycoord<-gridcoords[,2]
	
	neighmatx<-neighmat[[1]]
	neighmaty<-neighmat[[2]]
	
	dx<-resolution[1]
	dy<-resolution[2]

	nmesh<-length(xcoord)
	nmu<-length(mux)
	Imat<-as(Diagonal(x=rep(1,nmesh)),"dgCMatrix")

	if(is.null(loga)){
		loga<-rep(0,length(nmu))
	}
	
	#diffusion
	difenv<-Imat
	diag(difenv)<-sqrt(exp(logc))
	temp0<-difenv%*%(neighmatx/dx^2+neighmaty/dy^2)%*%difenv
	sumtemp0<-Imat
	diag(sumtemp0)<-c(as.matrix(t(rep(1,nmesh))%*%t(temp0)))

	res<-advdiff_core(neighmatx=neighmatx,
					neighmaty=neighmaty,
					temp0=temp0,
					sumtemp0=sumtemp0,
					xcoord=xcoord,
					ycoord=ycoord,
					mux=mux,
					muy=muy,
					loga=loga,
					dx=dx,dy=dy,
					logprob=logprob)

	return(res)
}

#log mean number of detection conditional on activity center
loglambda.ad<-function(type="binom",g0par,effort,effort_loc,effort_occtype,logparr.ad,logprob=T){
	nmu<-dim(logparr.ad)[1]
	neffort<-length(effort_loc)
	triplet<-cbind(rep(1:nmu,neffort),rep(effort_loc,each=nmu),rep(1,nmu*neffort),rep(effort_occtype,each=nmu))
	logparr.ad.reduce<-matrix(logparr.ad[triplet],nrow=nmu,ncol=neffort)
	if(type=="binom"){
		res<-outer(rep(1,nmu),g0par[1,,1])+logparr.ad.reduce
	}else if(type=="poisson"){
		res<-g0par[1,,1]+outer(rep(1,nmu),log(effort))+logparr.ad.reduce
	}else{
		stop("Observation model not implemented.\n")
	}
	if(logprob==F){
		res<-exp(res)
	}
	return(res)
}

logpzero.ad<-function(type="binom",loglambda,effort,logprob=T){
	nmu<-nrow(loglambda)
	neffort<-ncol(loglambda)
	zeromat<-matrix(0,nmu,neffort)
	if(type=="binom"){
		effortmat<-outer(rep(1,nmu),effort)
		res<-dbin_mat(zeromat,effortmat,log1m_exp_mat(-exp(loglambda)),logprob=logprob)
	}else if(type=="poisson"){
		res<-dpoisson_mat(zeromat,loglambda,logprob=logprob)
	}else{
		stop("Observation model not implemented.\n")
	}
	return(res)
}

#caluculate log 1-pdot from list of log(1-pdot_obs)
logpzero.x<-function(logpzerolist_occ_grp){
	nobs<-length(logpzerolist_occ_grp)
	ncell<-dim(logpzerolist_occ_grp[[1]])[1]
	res<-rep(0,ncell)
	for(i in 1:nobs){
		if(dim(logpzerolist_occ_grp[[i]])[2]!=0){
			res<-res+rowSums(logpzerolist_occ_grp[[i]])
		}
	}
	return(res)
}

#log1m_exp for multidimentional array
log1m_exp_arr<-function(arr){
	dimension<-dim(arr)
	ndim<-length(dimension)
	mat<-array(arr,dim=c(dimension[1],prod(dimension[-1])))
	out<-log1m_exp_mat(mat)
	res<-array(out,dim=dimension)
	return(res)
}

dbin_arr<-function(k,N,lp,logprob=T){
	cond<-all((dim(k)==dim(N))&(dim(N)==dim(lp)))
	if(!cond){
		stop("dimension mismatch.")
	}
	dimension<-dim(k)
	ndim<-length(dimension)
	dimleft<-prod(dimension[-1])
	kmat<-array(k,dim=c(dimension[1],dimleft))
	Nmat<-array(N,dim=c(dimension[1],dimleft))
	lpmat<-array(lp,dim=c(dimension[1],dimleft))
	out<-dbin_mat(kmat,Nmat,lpmat,logprob=logprob)
	res<-array(out,dim=dimension)
	return(res)
}


dpoisson_arr<-function(k,loglambda,logprob=T){
	cond<-all(dim(k)==dim(loglambda))
	if(!cond){
		stop("dimension mismatch.")
	}
	dimension<-dim(k)
	ndim<-length(dimension)
	dimleft<-prod(dimension[-1])
	kmat<-array(k,dim=c(dimension[1],dimleft))
	loglambdamat<-array(loglambda,dim=c(dimension[1],dimleft))
	out<-dpoisson_mat(kmat,loglambdamat,logprob=logprob)
	res<-array(out,dim=dimension)
	return(res)
}

logpcapthist<-function(type="binom",detect,loglambda,effort,effort_occ,srv,ind_cov){
	loglambda_dim<-dim(loglambda)
	loglambda_mat<-array(loglambda,dim=c(loglambda_dim[1],prod(loglambda_dim[-1])))
	if(type=="binom"){
		gc();gc()
		lp_mat<-log1m_exp_mat(-exp(loglambda_mat))
		res<-pcap_bin(detect,effort,effort_occ,lp_mat,srv,ind_cov)
		gc();gc()
	}else if(type=="poisson"){
		gc();gc()
		res<-pcap_poisson(detect,effort_occ,loglambda_mat,srv,ind_cov)
		gc();gc()
	}else{
		stop("observation model not implemented.\n")
	}
	return(res)
}

###animal movement simulation
advdiff.eigen.t<-function(t,initloc,mucoords,gridcoords,logc,loga,neighmat,resolution,nstep=10000,returnp=FALSE,steady=FALSE){

	nt<-length(t)
	difft<-diff(c(0,t))
	
	xcoord<-gridcoords[,1]
	ycoord<-gridcoords[,2]
	
	if(length(mucoords)==2){
		mux<-mucoords[1]
		muy<-mucoords[2]
	}else{
		mux<-xcoord[mucoords]
		muy<-ycoord[mucoords]
	}

	if(length(initloc)==2){
		initloc<-which((xcoord==initloc[1])&(ycoord==initloc[2]))
	}
	
	ncells<-length(xcoord)
	
	neighmatx<-neighmat[[1]]
	neighmaty<-neighmat[[2]]
	
	dx<-resolution[1]
	dy<-resolution[2]

	nmesh<-length(xcoord)
	nmu<-length(mux)
	Imat<-as(Diagonal(x=rep(1,nmesh)),"dgCMatrix")

	p0<-rep(0,ncells)
	p0[initloc]<-1
	
	#diffusion
	difenv<-Imat
	diag(difenv)<-sqrt(exp(logc))
	temp0<-difenv%*%(neighmatx/dx^2+neighmaty/dy^2)%*%difenv
	sumtemp0<-Imat
	diag(sumtemp0)<-c(as.matrix(t(rep(1,nmesh))%*%t(temp0)))
	if(returnp){
		pmat<-matrix(NA,ncells,nt)
	}
	res<-rep(NA,nt)
	advdiffout<-advdiff_t_core(T=difft,
						p0,
						mux=mux,
						muy=muy,
						neighmatx=neighmatx,
						neighmaty=neighmaty,
						temp0=temp0,
						sumtemp0=sumtemp0,
						xcoord=xcoord,
						ycoord=ycoord,
						loga=loga,
						dx=dx,dy=dy,
						nstep=nstep)

	res<-advdiffout[[1]]
	if(any(is.nan(advdiffout[[2]]))){
		stop("Probability NaN produced.\n")
	}else if(any(advdiffout[[2]]<0)){
		stop("Negative probability produced.\n")
	}
	if(returnp){
		attr(res,"pmat")<-advdiffout[[2]]
	}
	if(steady){
		pbar <- exp(advdiff_core(neighmatx,neighmaty,temp0,sumtemp0,xcoord,ycoord,mux,muy,loga,dx,dy))
		attr(res,"pbar")<-pbar
	}
	return(res)
}




##Class 'secrad'
secrad<-
R6Class("secrad",
	public=list(
		secrdata=NA,
		envmodel=list(D~1,C~1,A~0),
		indmodel=c(A=FALSE,g0=FALSE),
		occmodel=c(A=FALSE,g0=FALSE),
		adtemp=list(adpar=NA,adout=NA),
		chtemp=list(chpar=NA,chout=NA),
		lastpar=NULL,
		initialize=function(secrdata){
			self$secrdata<-secrdata
		},
		set_model=function(envmodel,indmodel,occmodel){
			self$envmodel<-envmodel
			self$indmodel<-indmodel
			self$occmodel<-occmodel
		},
		del_temp=function(){
			self$adtemp=list(adpar=NA,adout=NA)
			self$chtemp=list(chpar=NA,chout=NA)
		},
		loglf=function(par,loglfscale=1,verbose=F){
			secrdata<-self$secrdata
			envmodel<-self$envmodel
			indmodel<-self$indmodel
			occmodel<-self$occmodel
			self$lastpar<-par
			####prepare data for estimation
			ncell<-secrdata$ncell
			nocc<-secrdata$nocc
			nind<-secrdata$nind
			
			#design matrices of grid
			gridcov<-secrdata$grid_cov
			if(!is.null(gridcov)){
				if(any(is.na(gridcov))){
					stop("NA must not be included in covariates.\n")
				}
			}
			modellist<-lapply(envmodel,"[",-2)
			for(i in 1:length(modellist)){
				names(modellist)[i]<-as.character(envmodel[[i]][2])
			}
			denv.mat<-model.matrix(modellist$D,data=gridcov)
			denv.offset<-model.offset(model.frame(modellist$D,data=gridcov))
			if(is.null(denv.offset)){
				denv.offset<-rep(0,ncell)
			}
			
			cenv.mat<-model.matrix(modellist$C,data=gridcov)
			cenv.offset<-model.offset(model.frame(modellist$C,data=gridcov))
			if(is.null(cenv.offset)){
				cenv.offset<-rep(0,ncell)
			}
			
			woaenv<-F
			aenv.mat<-model.matrix(modellist$A,data=gridcov)
			if(length(aenv.mat)==0){
				aenv.mat<-matrix(0,nrow=ncell,ncol=1)
				woaenv<-T
			}else if(is.element("(Intercept)",colnames(aenv.mat))){
				stop("Intercept cannot be included in the advection submodel.")
			}
			aenv.offset<-model.offset(model.frame(modellist$A,data=gridcov))
			if(is.null(aenv.offset)){
				aenv.offset<-rep(0,ncell)
			}
			
			#individual covariate
			if(any(indmodel)){
				indcov<-secrdata$ind_cov
			}else{
				indcov<-rep(0,nind)
			}
			
			#occasion covariate
			if(any(occmodel)){
				occcov<-secrdata$occ_cov
			}else{
				occcov<-rep(0,nocc)
			}
				
			###partitioning parameters
			#density
			ndpar<-ncol(denv.mat)
			dpar<-par[1:ndpar]
			#connectivity
			ncpar<-ncol(cenv.mat)
			cpar<-par[(1:ncpar)+ndpar]
			#advection
			napar<-ifelse(woaenv,0,ncol(aenv.mat))
			if(napar!=0){
				apar<-apar2<-par[(1:napar)+ndpar+ncpar]
			}else{
				apar<-numeric(0)
				apar2<-0
			}
			#g0
			nobs<-length(secrdata$obs)
			obstype<-unlist(lapply(secrdata$obs,"[[","type"))
			npar_type<-c(1,1);names(npar_type)<-c("binom","poisson")
			ng0par_obs<-npar_type[obstype]
			g0par<-vector("list",nobs)
			ng0par<-0
			for(i in 1:nobs){
				g0par[[i]]<-par[(1:ng0par_obs[i])+ng0par+ndpar+ncpar+napar]
				ng0par<-ng0par+ng0par_obs[i]
			}
			#occasion
			noccgrp<-length(unique(occcov))
			noccpar<-sum(occmodel)*(noccgrp-1)
			if(noccpar!=0){
				occpar<-par[(1:noccpar)+ndpar+ncpar+napar+ng0par]
				names(occpar)<-paste0(rep(names(occmodel)[occmodel],each=noccgrp-1),"_",rep(levels(factor(occcov))[-1],sum(occmodel)))
			}else{
				occpar<-numeric(0)
			}
			
			nindgrp<-length(unique(indcov))	
			nindpar<-sum(indmodel)*(nindgrp-1)
			nmixpar<-nindgrp-1
			if(nindpar!=0){
				indpar<-par[(1:nindpar)+ndpar+ncpar+napar+ng0par+noccpar]
				names(indpar)<-paste0(rep(names(indmodel)[indmodel],each=nindgrp-1),"_",rep(levels(factor(indcov))[-1],sum(indmodel)))
				mixpar<-par[(1:nmixpar)+nindpar+ndpar+ncpar+napar+ng0par+noccpar]
			}else{
				indpar<-numeric(0)
				mixpar<-numeric(0)
			}
		
			if(!is.null(secrdata$auxiliary)){
				aux<-secrdata$auxiliary
			}else{
				aux<-NULL
			}
			#survival matrix(not implemented)
			srv<-matrix(1,nrow=secrdata$nind,ncol=nocc)
		#	if(!is.null(aux)){
		#		auxtype<-unlist(lapply(secrdata$auxiliary,"[","type"))
		#		if(is.element("removal",auxtype)){
		#			numremoval<-which(auxtype=="removal")
		#			for(i in 1:length(numremoval)){
		#				whichremove<-which(secrdata$auxiliary[[numremoval[i]]]$detect==1,arr.ind=T)
		#				doublet<-cbind(ind=whichremove[,2],occ=secrdata$auxiliary[[numremoval[i]]]$effort_occ[whichremove[,1]])
		#				for(j in 1:nrow(doublet)){
		#					if(doublet[j,2]<nocc){
		#						srv[doublet[j,1],(doublet[j,2]+1):nocc]<-0
		#					}
		#				}
		#			}
		#		}
		#	}
			
			#
			if(verbose){
				cat("[parameters]\n","density ",dpar,"\n","connectivity ", cpar,
					"\n","advection ",apar,"\n","detectability ",unlist(g0par),
					"\n","occasion ",occpar,
					"\n","individual ",indpar,"\n","mixing parameter ",mixpar,"\n")
			}
			#log mean density
			lnD<-c(denv.mat%*%dpar+denv.offset)
			
			#individual covariate of advection
			grepA<-grep("^A",names(indpar))
		
			if(length(grepA)!=0){
				apar_ind<-c(0,indpar[grepA])
			}else{
				apar_ind<-0
			}
			
			#individual covariate of detection
			grepg0<-grep("^g0",names(indpar))
			if(length(grepg0)!=0){
				g0par_ind<-c(0,indpar[grepg0])
			}else{
				g0par_ind<-0
			}
			ng0par_ind<-length(g0par_ind)
		
			#occasion covariate of advection
			grepAocc<-grep("^A",names(occpar))
		
			if(length(grepAocc)!=0){
				apar_occ<-c(0,occpar[grepAocc])
			}else{
				apar_occ<-0
			}
			
			#occasion covariate of detection
			grepg0occ<-grep("^g0",names(occpar))
			if(length(grepg0occ)!=0){
				g0par_occ<-c(0,occpar[grepg0occ])
			}else{
				g0par_occ<-0
			}
			ng0par_occ<-length(g0par_occ)
		
			
			#list of g0 coef
			g0list<-vector("list",ng0par)
			for(i in 1:nobs){
				g0list[[i]]<-outer(matrix(1,ng0par_ind,ng0par_occ),g0par[[i]])
				g0list[[i]][,,1]<-g0list[[i]][,,1]+outer(g0par_ind,rep(1,ng0par_occ))+outer(rep(1,ng0par_ind),g0par_occ)
			}
		
			#advection-diffusion model
			if(any(is.na(self$adtemp$adpar))){
				self$adtemp$adpar<-array(NA,dim=c(ncpar+napar+2,length(apar_ind),length(apar_occ)))
				self$adtemp$adout<-array(NA,dim=c(ncell,ncell,length(apar_ind),length(apar_occ)))
			}
			logp_arr<-array(NA,dim=c(ncell,ncell,length(apar_ind),length(apar_occ)))
			for(i in 1:length(apar_ind)){
				for(j in 1:length(apar_occ)){
					flag<-TRUE
					adpar<-c(cpar,apar,apar_ind[i],apar_occ[j])
					if(all(!is.na(self$adtemp$adpar[,i,j]))){
						adpar_pre<-self$adtemp$adpar[,i,j]
						if(all(adpar==adpar_pre)){
							flag<-FALSE
						}
					}
					if(flag){
						logc<-c(cenv.mat%*%matrix(cpar))+cenv.offset
						loga<-c(aenv.mat%*%matrix(apar2))+aenv.offset+apar_ind[i]+apar_occ[j]
						logp_arr[,,i,j]<-advdiff.eigen(mucoords=secrdata$coords,
										gridcoords=secrdata$coords,
										logc=logc,
										loga=loga,
										neighmat=secrdata$neighmat,
										resolution=secrdata$resolution,
										logprob=T)
						self$adtemp$adpar[,i,j]<-adpar
						self$adtemp$adout[,,i,j]<-logp_arr[,,i,j]
					}else{
						logp_arr[,,i,j]<-self$adtemp$adout[,,i,j]
					}
				}
			}
				
			#mean number of detection (or cloglog of prob for binimial)
			loglambdalist<-vector("list",nobs)
			nindgrp_lambda<-pmax(length(apar_ind),length(g0par_ind))
			for(i in 1:nobs){
				loglambdalist[[i]]<-array(NA,dim=c(ncell,length(secrdata$obs[[i]]$effort),nindgrp_lambda))
				A_effort_occtype<- ifelse(rep(occmodel["A"],length(secrdata$obs[[i]]$effort_occ)),
								occcov[secrdata$obs[[i]]$effort_occ]+1,
								rep(1,length(secrdata$obs[[i]]$effort_occ)))
				g0_effort_occtype<-ifelse(rep(occmodel["g0"],length(secrdata$obs[[i]]$effort_occ)),
								occcov[secrdata$obs[[i]]$effort_occ]+1,
								rep(1,length(secrdata$obs[[i]]$effort_occ)))
				for(j in 1:nindgrp_lambda){
					if(indmodel["A"]){
						logp<-logp_arr[,,j,,drop=F]
					}else{
						logp<-logp_arr[,,1,,drop=F]
					}
					g0list_indocc<-g0list[[i]][j,g0_effort_occtype,,drop=F]
					loglambdalist[[i]][,,j]<-loglambda.ad(type=obstype[i],g0list_indocc,secrdata$obs[[i]]$effort,secrdata$obs[[i]]$effort_loc,A_effort_occtype,logp,logprob=T)
				}
			}
			
			####calculation of pdot
			#log probability of non-detection of each individual group and occasion group
			logpzerolist<-vector("list",nobs)
			for(i in 1:nobs){
				logpzerolist[[i]]<-array(NA,dim=c(ncell,length(secrdata$obs[[i]]$effort),nindgrp))
				for(j in 1:nindgrp){
						logpzerolist[[i]][,,j]<-logpzero.ad(type=obstype[i],loglambdalist[[i]][,,j],secrdata$obs[[i]]$effort)
				}
			}
			
			#log 1-pdot for each occasion conditional on X
			logpzero_occ_x_arr<-array(NA,dim=c(ncell,nocc,nindgrp))
			for(i in 1:nocc){
				logpzerolist_occ<-vector("list",nobs)
				for(j in 1:nobs){
					cond<-secrdata$obs[[j]]$effort_occ==i
					logpzerolist_occ[[j]]<-logpzerolist[[j]][,cond,,drop=F]
				}
				
				for(k in 1:nindgrp){
					logpzerolist_occ_grp<-lapply(logpzerolist_occ,"[",,,k)
					logpzero_occ_x_arr[,i,k]<-logpzero.x(logpzerolist_occ_grp)
				}
			}
			
		
			#ln(pdot) for each individual
			logpzero_occ_x_ind<-logpzero_occ_x_arr[,,indcov-min(indcov)+1,drop=F]
			logpzero_x_ind<-(logpzero_occ_x_ind*outer(rep(1,ncell),t(srv)))%>%
						aperm(perm=c(1,3,2))%>%
						rowSums(dim=2)
			logpdot_x_ind<-log1m_exp_mat(logpzero_x_ind)			
		
			#lnD+ln(area)+ln(pdot) for each individual
			logd_area_pdot_x_ind<-outer(lnD+log(secrdata$area),rep(1,nind))+logpdot_x_ind
			
			#logsumexp(lnD+ln(area)+ln(pdot)) for each individual
			logd_area_pdot_ind<-colsum_lp(logd_area_pdot_x_ind)
			
			#####likelihood of individual capture history conditional on detection at least once
			#multinomial coefficient
			detectlist<-lapply(secrdata$obs,"[[","detect")
			detect_append<-detectlist[[1]]
			if(nobs>1){
				for(i in 2:nobs){
					detect_append<-rbind(detect_append,detectlist[[i]])
				}
			}
			detect_pattern<-apply(detect_append,2,paste,collapse="_")
			nc<-table(detect_pattern)
			names(nc)<-NULL
			logmulticoef<-lgamma(nind+1)-sum(lgamma(nc+1))
			
			#log prob of capture history conditional on X
			if(any(is.na(self$chtemp$chpar))){
				self$chtemp$chpar<-array(NA,dim=c(ncpar+napar+length(apar_ind)+length(apar_occ)+length(g0list[[1]]),nobs))
				self$chtemp$chout<-array(NA,dim=c(ncell,nind,nobs))
			}
			logpcapthist_arr<-array(0,dim=c(ncell,nind,nobs))
			for(i in 1:nobs){
				flagch<-TRUE
				chpar<-c(cpar,apar,apar_ind,apar_occ,c(g0list[[i]]))
				if(all(!is.na(self$chtemp$chpar[,i]))){
					chpar_pre<-self$chtemp$chpar[,i]
					if(all(chpar==chpar_pre)){
						flagch<-FALSE
					}
				}
				
				if(flagch){
					logpcapthist_arr[,,i]<-logpcapthist(type=obstype[i],
							detect=secrdata$obs[[i]]$detect,
							loglambda=loglambdalist[[i]],
							effort=secrdata$obs[[i]]$effort,
							effort_occ=secrdata$obs[[i]]$effort_occ,
							srv=srv,
							ind_cov=indcov)
					self$chtemp$chpar[,i]<-chpar
					self$chtemp$chout[,,i]<-logpcapthist_arr[,,i]
				}else{
					logpcapthist_arr[,,i]<-self$chtemp$chout[,,i]
				}
			}
			
			logpcapthist_mat<-rowSums(logpcapthist_arr,dim=2)
			
			#multinomial log-likelihood(pdot is cancelled out)
			loglfmulti<-logmulticoef+sum(colsum_lp(logpcapthist_mat+outer(lnD+log(secrdata$area),rep(1,nind))-outer(rep(1,ncell),logd_area_pdot_ind)))
			
			####Poisson process likelihood (removal not implemented yet)
			logpdot_x_grp <- logpzero_occ_x_arr%>%
						aperm(perm=c(1,3,2))%>%
						rowSums(dim=2)%>%
						log1m_exp_mat()
			logmix<-c(0,mixpar)-sum_lp(c(0,mixpar))
			lambda_grp <- colsum_lp(outer(lnD+log(secrdata$area),rep(1,nindgrp))+outer(rep(1,ncell),logmix)+logpdot_x_grp)
			loglfpois<-dpois(table(indcov),exp(lambda_grp),log=T)
			names(loglfpois)<-NULL
			res<-(sum(loglfpois)+loglfmulti)*loglfscale
			if(verbose){
				cat("[loglf]","\n",res,"\n",
					"[population size]","\n",sum(exp(lnD)*secrdata$area),"\n")
			}
			return(res)
		}
	)
)
		
##Class 'secrad_data'
secrad_data<-
R6Class("secrad_data",
	public=list(
		nind=0,
		nocc=0,
		ncell=0,
		obs=as.list(numeric(0)),
		auxiliary=as.list(numeric(0)),
		area=1,
		ind_cov=NULL,
		grid_cov=NULL,
		occ_cov=NULL,
		coords=cbind(x=1,y=1),
		neighmat=list(x=as(Matrix(0),"dgCMatrix"),y=as(Matrix(0),"dgCMatrix")),
		resolution=c(x=1,y=1),
		truemodel=NULL,
		sim_settings=NULL,
		trueind=NULL,
		sim_result=NULL,
		initialize=function(coords,area=NULL,grid_cov=NULL,resolution=NULL){
			ncell<-nrow(coords)
			if(nrow(grid_cov)!=ncell){
				stop("Lengths of coords and grid_cov are not equal")
			}
			self$ncell<-ncell
			self$coords<-coords
			self$area<-area
			self$grid_cov<-grid_cov
			if(!is.null(resolution)){
				self$resolution<-resolution
			}
			
			distx<-abs(outer(rep(1,ncell),coords[,"x"])-outer(coords[,"x"],rep(1,ncell)))
			disty<-abs(outer(rep(1,ncell),coords[,"y"])-outer(coords[,"y"],rep(1,ncell)))
			neighmatx<-as(Matrix(((distx<resolution["x"]*1.5)&(disty<resolution["y"]*0.5))*1),"dgCMatrix")
			diag(neighmatx)<-0
			neighmaty<-as(Matrix(((disty<resolution["y"]*1.5)&(distx<resolution["x"]*0.5))*1),"dgCMatrix")
			diag(neighmaty)<-0
			self$neighmat=list(x=neighmatx,y=neighmaty)			
		},
		add_obs=function(type="poisson",effort,effort_loc,effort_occ,detect=NULL){
			neffort<-length(effort)
			if(neffort!=length(effort_loc)|neffort!=length(effort_occ)){
				stop("Lengths of effort, effort_loc and effort_occ are not equal.")
			}
			maxocc<-max(effort_occ)
			self$nocc<-pmax(self$nocc,maxocc)
			nobs<-length(self$obs)
			effort_coords<-cbind(x=self$coords[effort_loc,"x"],y=self$coords[effort_loc,"y"])
			if(is.null(detect)){
				detect<-matrix(numeric(0),nrow=neffort,ncol=0)
			}
			self$obs[[nobs+1]]<-list(type=type,detect=detect,effort=effort,effort_loc=effort_loc,effort_occ=effort_occ,effort_coords=effort_coords)
			self$nind=ncol(detect)
		},
		set_occ_cov=function(occ_cov){
			nocc_cov<-length(occ_cov)
			if(nocc_cov!=self$nocc){
				stop("Length of occasion covariate and number of occasions are not equal.")
			}
			self$occ_cov=occ_cov
		},
		set_ind_cov=function(ind_cov){
			nind_cov<-length(ind_cov)
			if(nind_cov!=self$nind){
				stop("Length of individual covariate and number of individuals are not equal.")
			}
			self$ind_cov=ind_cov
		},
		set_truemodel=function(envmodel,indmodel,occmodel,dpar,cpar,adpar,g0par,indpar=numeric(0),occpar=numeric(0),mixpar=numeric(0)){
			nenvmodel<-length(envmodel)
			for(i in 1:nenvmodel){
				names(envmodel)[i]<-as.character(envmodel[[i]][2])
			}
			self$truemodel=list(envmodel=envmodel,indmodel=indmodel,occmodel=occmodel,dpar=dpar,cpar=cpar,adpar=adpar,g0par=g0par,indpar=indpar,occpar=occpar,mixpar=mixpar)
		},
		set_sim=function(timeperocc,stepspertime,stepad,timeburnin){
			self$sim_settings=list(timeperocc=timeperocc,stepspertime=stepspertime,stepad=stepad,timeburnin=timeburnin)
		},
		generate_ind=function(){
			if(is.null(self$truemodel)){
				stop("True model is required.")
			}
			densmodel<-self$truemodel$envmodel$D[-2]
			if(!is.null(self$grid_cov)){
				denv.mat<-model.matrix(densmodel,data=self$grid_cov)
				denv.offset<-model.offset(model.frame(densmodel,data=self$grid_cov))
				if(is.null(denv.offset)){
					denv.offset<-rep(0,self$ncell)
				}
				lnD<-c(denv.mat%*%self$truemodel$dpar+denv.offset)
			}else{
				lnD<-rep(self$truemodel$dpar[1],ncell)
			}
			lambda<-exp(lnD)*self$area
			lambdatotal<-sum(lambda)
			N<-rpois(1,lambdatotal)
			ind_loc<-apply(rmultinom(N,1,lambda/lambdatotal)==1,2,which)
			if(length(self$truemodel$mixpar)!=0){
				rpmix<-exp(c(0,self$truemodel$mixpar))
				ind_cov<-apply(rmultinom(N,1,rpmix/sum(rpmix))==1,2,which)-1
			}else{
				ind_cov<-rep(0,N)
			}
			self$trueind<-cbind(ind_loc=ind_loc,ind_cov=ind_cov)
		},
		simulate=function(return=FALSE,verbose=TRUE){
			nobs<-length(self$obs)
			nocc<-self$nocc
			nindall<-nrow(self$trueind)
			cmodel<-self$truemodel$envmodel$C[-2]
			amodel<-self$truemodel$envmodel$A[-2]
			if(!is.null(self$grid_cov)){
				cenv.mat<-model.matrix(cmodel,data=self$grid_cov)
				cenv.offset<-model.offset(model.frame(cmodel,data=self$grid_cov))
				if(is.null(cenv.offset)){
					cenv.offset<-rep(0,self$ncell)
				}
				logc<-c(cenv.mat%*%self$truemodel$cpar+cenv.offset)
				aenv.mat<-model.matrix(amodel,data=self$grid_cov)
				aenv.offset<-model.offset(model.frame(amodel,data=self$grid_cov))
				if(is.null(aenv.offset)){
					aenv.offset<-rep(0,self$ncell)
				}
				loga<-c(aenv.mat%*%self$truemodel$adpar+aenv.offset)
			}else{
				logc<-rep(self$truemodel$cpar[1],self$ncell)
				loga<-rep(self$truemodel$adpar[1],self$ncell)
			}
			if(self$truemodel$occmodel["A"]){
				a.occeff<-c(0,self$truemodel$occpar["A"])[self$occ_cov+1]
			}else{
				a.occeff<-rep(0,self$nocc)
			}
			if(self$truemodel$occmodel["g0"]){
				g0.occeff<-c(0,self$truemodel$occpar["g0"])[self$occ_cov+1]
			}else{
				g0.occeff<-rep(0,self$nocc)
			}
			if(self$truemodel$indmodel["A"]){
				a.indeff<-c(0,self$truemodel$indpar["A"])[self$trueind[,"ind_cov"]+1]
			}else{
				a.indeff<-rep(0,nindall)
			}
			if(self$truemodel$indmodel["g0"]){
				g0.indeff<-c(0,self$truemodel$indpar["g0"])[self$trueind[,"ind_cov"]+1]
			}else{
				g0.indeff<-rep(0,nindall)
			}
			
			stepspertime<-self$sim_settings$stepspertime
			timeperocc<-self$sim_settings$timeperocc
			nocc<-self$nocc
			timeburnin<-self$sim_settings$timeburnin
			stepad<-self$sim_settings$stepad
				
			typevec<-unlist(lapply(self$obs,"[[","type"))
			nobs<-length(self$obs)
			res<-vector("list",nobs)
			for(k in 1:nobs){
				res[[k]]<-matrix(NA,nrow=length(self$obs[[k]]$effort),ncol=nindall)
			}
			traj<-matrix(NA,nrow=stepspertime*timeperocc*nocc+1,ncol=nindall)
			occlab<-c(NA,rep(1:nocc,each=stepspertime*timeperocc))
			for(i in 1:nindall){
				if(verbose){
					cat(i,"\n")
				}
				loga.ind<-loga[self$trueind[i,"ind_loc"]]+a.indeff[i]+a.occeff[1]
				
				#burnin
				if(timeburnin!=0){
					loc_p0<-advdiff.eigen.t(t=timeburnin,initloc=self$trueind[i,"ind_loc"],mucoords=self$trueind[i,"ind_loc"],gridcoords=self$coords,logc=logc,loga=loga.ind,neighmat=self$neighmat,resolution=self$resolution,nstep=stepad,returnp=F,steady=F)
				}else{
					loc_p0<-self$trueind[i,"ind_loc"]
				}
				traj[1,i]<-loc_p0
				#simulate movement
				nvisitmat<-matrix(0,self$ncell,self$nocc)
				for(j in 1:self$nocc){
					loga.ind<-loga[self$trueind[i,"ind_loc"]]+a.indeff[i]+a.occeff[j]
					t<-seq(0,timeperocc,length=stepspertime*timeperocc+1)
					loc_occ<-advdiff.eigen.t(t=t,initloc=loc_p0,mucoords=self$trueind[i,"ind_loc"],gridcoords=self$coords,logc=logc,loga=loga.ind,neighmat=self$neighmat,resolution=self$resolution,nstep=stepad,returnp=F,steady=F)
					traj[1+(1:(stepspertime*timeperocc))+(stepspertime*timeperocc)*(j-1),i]<-loc_occ[-1]
					nvisit0<-table(loc_occ[-1])
					nvisitmat[as.numeric(names(nvisit0)),j]<-nvisit0
				}
				#simulate detection
				for(k in 1:nobs){
					g0.occ<-self$truemodel$g0par[k]+g0.indeff[i]+g0.occeff[self$obs[[k]]$effort_occ]
					nvisit_effort<-nvisitmat[cbind(self$obs[[k]]$effort_loc,self$obs[[k]]$effort_occ)]
					if(typevec[k]=="poisson"){
						lambda<-exp(g0.occ)*self$obs[[k]]$effort*nvisit_effort/timeperocc/stepspertime
						res[[k]][,i]<-rpois(length(self$obs[[k]]$effort),lambda)
					}else if(typevec[k]=="binom"){
						lambda<-exp(g0.occ)*nvisit_effort/timeperocc/stepspertime
						p<-1-exp(-lambda)
						res[[k]][,i]<-rbinom(length(self$obs[[k]]$effort),self$obs[[k]]$effort,p)
					}
				}
			}
			#output
			self$sim_result<-vector("list",2);names(self$sim_result)<-c("detect","trajectory")
			self$sim_result$detect<-res
			self$sim_result$trajectory<-vector("list",2);names(self$sim_result$trajectory)<-c("trajectory","occasion")
			self$sim_result$trajectory$trajectory<-traj
			self$sim_result$trajectory$occasion<-occlab
			if(return){
				return(self)
			}
		},
		update_detect=function(){
			if(is.null(self$sim_result)){
				stop("Simulation results required.")
			}
			nobs<-length(self$obs)
			nindall<-nrow(self$trueind)
			inddetect<-rep(FALSE,nindall)
			for(i in 1:nobs){
				inddetect<-inddetect|(apply(self$sim_result$detect[[i]],2,sum)!=0)
			}
			self$ind_cov<-self$trueind[inddetect,"ind_cov"]
			self$nind<-sum(inddetect)
			for(i in 1:nobs){
				self$obs[[i]]$detect<-self$sim_result$detect[[i]][,inddetect]
			}
		},
		out_truepar=function(){
			if(is.null(self$truemodel)){
				stop("True model required.\n")
			}
			dpar<-self$truemodel$dpar
			names(dpar)<-paste0("dens",0:(length(dpar)-1))
			cpar<-self$truemodel$cpar
			names(cpar)<-paste0("con",0:(length(cpar)-1))
			apar<-self$truemodel$adpar
			names(apar)<-paste0("adv",0:(length(apar)-1))
			g0par<-self$truemodel$g0par
			names(g0par)<-paste0("det",1:length(g0par))
			indpar<-self$truemodel$indpar
			if(length(indpar)!=0){
				names(indpar)<-paste0("ind",0:(length(indpar)-1))
			}
			occpar<-self$truemodel$occpar
			if(length(occpar)!=0){
				names(occpar)<-paste0("occ",1:(length(occpar)))
			}
			mixpar<-self$truemodel$mixpar
			if(length(mixpar)!=0){
				names(mixpar)<-paste0("mix",1:(length(mixpar)))
			}
			return(c(dpar,cpar,apar,g0par,indpar,occpar,mixpar))
		},			
		ggsecraddata=function(covname=NULL,jitter=0,sample=1){
			if(is.null(self$obs)){
				stop("Observations required.")
			}
			if(is.null(covname)){
				z<-rep(0,self$ncell)
				covname=""
			}else{
				z<-unlist(self$grid_cov[,covname])
			}
			tiledata<-tibble(x=self$coords[,"x"],y=self$coords[,"y"],z=z)
			
			nobs<-length(self$obs)
			nind<-self$nind
			alldetect<-tibble(obs=numeric(0),ind=numeric(0),x=numeric(0),y=numeric(0),n=numeric(0),occ=numeric(0))
			alltrap<-tibble(obs=numeric(0),x=numeric(0),y=numeric(0),occ=numeric(0))
			for(i in 1:nobs){
				focaldata<-self$obs[[i]]$detect
				focalloc<-self$obs[[i]]$effort_loc
				focalocc<-self$obs[[i]]$effort_occ
				trapx<-self$coords[focalloc,"x"]
				trapy<-self$coords[focalloc,"y"]
				alltrap<-bind_rows(alltrap,tibble(obs=rep(i,length(trapx)),x=trapx,y=trapy,occ=focalocc))
		
				for(j in 1:nind){
					n<-focaldata[focaldata[,j]>0,j]
					loc<-focalloc[focaldata[,j]>0]
					x<-self$coords[loc,"x"]
					y<-self$coords[loc,"y"]
					obstype<-rep(i,length(x))
					indnum<-rep(j,length(x))
					occnum<-focalocc[focaldata[,j]>0]
					alldetect<-bind_rows(alldetect,tibble(obs=obstype,ind=indnum,x=x,y=y,n=n,occ=occnum))
				}
			}
			sampleind<-sample(1:nind,round(nind*sample))
			alldetect<-alldetect%>%filter(is.element(ind,sampleind))
			inddetect<-alldetect%>%st_as_sf(coords=c("x","y"))%>%
					st_jitter(factor=jitter)%>%
					suppressWarnings(st_centroid)%>%group_by(ind)%>%
					summarize(union=st_union(geometry))%>%
					mutate(chull=st_convex_hull(union))
		
			ggplot()+geom_tile(aes(x=x,y=y,fill=z),data=tiledata, width=self$resolution["x"], height=self$resolution["y"])+
				geom_point(data=alltrap,aes(x=x,y=y),size=3,pch=0,stroke=1.5)+
				geom_sf(data=inddetect,aes(col=as.character(ind)),show.legend = FALSE)+
				geom_sf(data=inddetect%>%st_drop_geometry%>%dplyr::select(ind,chull)%>%st_as_sf,aes(col=as.character(ind)),fill=NA,show.legend = FALSE)+
				scale_fill_gradientn(colours = rev(terrain.colors(20))[-1])+
				labs(fill = covname)+
				theme_cowplot()
		},
		ggsecradsim=function(ind=1,tile=NULL,occ=1,ptrans="identity"){
			if(is.null(self$sim_result)){
				stop("Simulation results required.")
			}
			mux<-self$coords[self$trueind[ind,"ind_loc"],"x"]
			muy<-self$coords[self$trueind[ind,"ind_loc"],"y"]
			nindall<-ncol(self$sim_result$detect[[1]])
			nobs<-length(self$sim_result$detect)
			
			if(is.null(tile)){
				cmodel<-self$truemodel$envmodel$C[-2]
				amodel<-self$truemodel$envmodel$A[-2]
				if(!is.null(self$grid_cov)){
					cenv.mat<-model.matrix(cmodel,data=self$grid_cov)
					cenv.offset<-model.offset(model.frame(cmodel,data=self$grid_cov))
					if(is.null(cenv.offset)){
						cenv.offset<-rep(0,self$ncell)
					}
					logc<-c(cenv.mat%*%self$truemodel$cpar+cenv.offset)
					aenv.mat<-model.matrix(amodel,data=self$grid_cov)
					aenv.offset<-model.offset(model.frame(amodel,data=self$grid_cov))
					if(is.null(aenv.offset)){
						aenv.offset<-rep(0,self$ncell)
					}
					loga<-c(aenv.mat%*%self$truemodel$adpar+aenv.offset)
				}else{
					logc<-rep(self$truemodel$cpar[1],self$ncell)
					loga<-rep(self$truemodel$adpar[1],self$ncell)
				}
				if(self$truemodel$occmodel["A"]){
					a.occeff<-c(0,self$truemodel$occpar["A"])[self$occ_cov+1]
				}else{
					a.occeff<-rep(0,self$nocc)
				}
				if(self$truemodel$indmodel["A"]){
					a.indeff<-c(0,self$truemodel$indpar["A"])[self$trueind[,"ind_cov"]+1]
				}else{
					a.indeff<-rep(0,nindall)
				}
				resolution<-self$resolution
				names(resolution)<-NULL
				loga.ind<-loga[self$trueind[ind,"ind_loc"]]+a.indeff[ind]+a.occeff[occ]
				z<-c(advdiff.eigen(mucoords=cbind(x=mux,y=muy),
								gridcoords=self$coords,
								logc=logc,loga=loga.ind,
								neighmat=self$neighmat,resolution=resolution,logprob=F))
				covname<-"prob"
				tiletrans<-ptrans
			}else{
				z<-unlist(self$grid_cov[,covname])
				tiletrans<-"identity"
			}
			tiledata<-tibble(x=self$coords[,"x"],y=self$coords[,"y"],z=z)
			
			traj<-as_tibble(self$coords[self$sim_result$trajectory$trajectory[,ind],])
			
			inddetect<-tibble(obs=numeric(0),x=numeric(0),y=numeric(0),n=numeric(0),occ=numeric(0))
			alltrap<-tibble(obs=numeric(0),x=numeric(0),y=numeric(0),occ=numeric(0))
			for(i in 1:nobs){
				focaldata<-self$sim_result$detect[[i]][,ind]
				focalloc<-self$obs[[i]]$effort_loc
				focalocc<-self$obs[[i]]$effort_occ
				trapx<-self$coords[focalloc,"x"]
				trapy<-self$coords[focalloc,"y"]
				alltrap<-bind_rows(alltrap,tibble(obs=rep(i,length(trapx)),x=trapx,y=trapy,occ=focalocc))
				
				if(any(focaldata>0)){
					n<-focaldata[focaldata>0]
					loc<-focalloc[focaldata>0]
					x<-self$coords[loc,"x"]
					y<-self$coords[loc,"y"]
					obstype<-rep(i,length(x))
					indnum<-rep(ind,length(x))
					occnum<-focalocc[focaldata>0]
					inddetect<-bind_rows(inddetect,tibble(obs=obstype,ind=indnum,x=x,y=y,n=n,occ=occnum))
				}
			}
			ggplot()+geom_tile(aes(x=x,y=y,fill=z),data=tiledata, width=self$resolution["x"], height=self$resolution["y"])+
				geom_point(data=alltrap,aes(x=x,y=y),size=3,pch=0,stroke=1.5)+
				geom_path(data=traj,aes(x=x,y=y),alpha=0.5)+
				geom_point(data=inddetect,aes(x=x,y=y,size=n),col="blue")+
				scale_fill_gradientn(colours = rev(terrain.colors(20))[-1],trans=tiletrans)+
				labs(fill = covname)+
				theme_cowplot()
		}	

	)
)

##plot 'secrad_data'
ggsecraddata<-function(addata,covname=NULL,jitter=0,sample=1){
	if(is.null(covname)){
		z<-rep(0,addata$ncell)
		covname=""
	}else{
		z<-addata$grid_cov[,covname]
	}
	tiledata<-tibble(x=addata$coords[,"x"],y=addata$coords[,"y"],z=z)
	
	nobs<-length(addata$obs)
	nind<-addata$nind
	alldetect<-tibble(obs=numeric(0),ind=numeric(0),x=numeric(0),y=numeric(0),n=numeric(0),occ=numeric(0))
	alltrap<-tibble(obs=numeric(0),x=numeric(0),y=numeric(0),occ=numeric(0))
	for(i in 1:nobs){
		focaldata<-addata$obs[[i]]$detect
		focalloc<-addata$obs[[i]]$effort_loc
		focalocc<-addata$obs[[i]]$effort_occ
		trapx<-addata$coords[focalloc,"x"]
		trapy<-addata$coords[focalloc,"y"]
		alltrap<-bind_rows(alltrap,tibble(obs=rep(i,length(trapx)),x=trapx,y=trapy,occ=focalocc))
		
		for(j in 1:nind){
			n<-focaldata[focaldata[,j]>0,j]
			loc<-focalloc[focaldata[,j]>0]
			x<-addata$coords[loc,"x"]
			y<-addata$coords[loc,"y"]
			obstype<-rep(i,length(x))
			indnum<-rep(j,length(x))
			occnum<-focalocc[focaldata[,j]>0]
			alldetect<-bind_rows(alldetect,tibble(obs=obstype,ind=indnum,x=x,y=y,n=n,occ=occnum))
		}
	}
	sampleind<-sample(1:nind,round(nind*sample))
	alldetect<-alldetect%>%filter(is.element(ind,sampleind))
	inddetect<-alldetect%>%st_as_sf(coords=c("x","y"))%>%
			st_jitter(factor=jitter)%>%
			suppressWarnings(st_centroid)%>%group_by(ind)%>%
			summarize(union=st_union(geometry))%>%
			mutate(chull=st_convex_hull(union))

	ggplot()+geom_tile(aes(x=x,y=y,fill=z),data=tiledata)+
		geom_point(data=alltrap,aes(x=x,y=y),size=3,pch=0,stroke=1.5)+
		geom_sf(data=inddetect,aes(col=as.character(ind)),show.legend = FALSE)+
		geom_sf(data=inddetect%>%st_drop_geometry%>%dplyr::select(ind,chull)%>%st_as_sf,aes(col=as.character(ind)),fill=NA,show.legend = FALSE)+
		scale_fill_gradientn(colours = rev(terrain.colors(20))[-1])+
		labs(fill = covname)+
		theme_cowplot()
}

ggsecradsim<-function(addata,ind=1,tile=NULL,occ=1,ptrans="identity"){
	mux<-addata$coords[addata$trueind[ind,"ind_loc"],"x"]
	muy<-addata$coords[addata$trueind[ind,"ind_loc"],"y"]
	
	if(is.null(tile)){
		cmodel<-addata$truemodel$envmodel$C[-2]
		amodel<-addata$truemodel$envmodel$A[-2]
		if(!is.null(addata$grid_cov)){
			cenv.mat<-model.matrix(cmodel,data=addata$grid_cov)
			cenv.offset<-model.offset(model.frame(cmodel,data=addata$grid_cov))
			if(is.null(cenv.offset)){
				cenv.offset<-rep(0,addata$ncell)
			}
			logc<-c(cenv.mat%*%addata$truemodel$cpar+cenv.offset)
			aenv.mat<-model.matrix(amodel,data=addata$grid_cov)
			aenv.offset<-model.offset(model.frame(amodel,data=addata$grid_cov))
			if(is.null(aenv.offset)){
				aenv.offset<-rep(0,addata$ncell)
			}
			loga<-c(aenv.mat%*%addata$truemodel$adpar+aenv.offset)
		}else{
			logc<-rep(addata$truemodel$cpar[1],addata$ncell)
			loga<-rep(addata$truemodel$adpar[1],addata$ncell)
		}
		if(addata$truemodel$occmodel["A"]){
			a.occeff<-c(0,addata$truemodel$occpar["A"])[addata$occ_cov+1]
		}else{
			a.occeff<-rep(0,addata$nocc)
		}
		if(addata$truemodel$indmodel["A"]){
			a.indeff<-c(0,addata$truemodel$indpar["A"])[addata$trueind[,"ind_cov"]+1]
		}else{
			a.indeff<-rep(0,nindall)
		}
		resolution<-addata$resolution
		names(resolution)<-NULL
		loga.ind<-loga[addata$trueind[ind,"ind_loc"]]+a.indeff[ind]+a.occeff[occ]
		z<-c(advdiff.eigen(mucoords=cbind(x=mux,y=muy),
						gridcoords=addata$coords,
						logc=logc,loga=loga.ind,
						neighmat=addata$neighmat,resolution=resolution,logprob=F))
		covname<-"prob"
		tiletrans<-ptrans
	}else{
		z<-addata$grid_cov[,covname]
		tiletrans<-"identity"
	}
	tiledata<-tibble(x=addata$coords[,"x"],y=addata$coords[,"y"],z=z)
	
	traj<-as_tibble(addata$coords[addata$sim_result$trajectory$trajectory[,ind],])
	
	inddetect<-tibble(obs=numeric(0),x=numeric(0),y=numeric(0),n=numeric(0),occ=numeric(0))
	alltrap<-tibble(obs=numeric(0),x=numeric(0),y=numeric(0),occ=numeric(0))
	for(i in 1:nobs){
		focaldata<-addata$sim_result$detect[[i]][,ind]
		focalloc<-addata$obs[[i]]$effort_loc
		focalocc<-addata$obs[[i]]$effort_occ
		trapx<-addata$coords[focalloc,"x"]
		trapy<-addata$coords[focalloc,"y"]
		alltrap<-bind_rows(alltrap,tibble(obs=rep(i,length(trapx)),x=trapx,y=trapy,occ=focalocc))
		
		if(any(focaldata>0)){
			n<-focaldata[focaldata>0]
			loc<-focalloc[focaldata>0]
			x<-addata$coords[loc,"x"]
			y<-addata$coords[loc,"y"]
			obstype<-rep(i,length(x))
			indnum<-rep(ind,length(x))
			occnum<-focalocc[focaldata>0]
			inddetect<-bind_rows(inddetect,tibble(obs=obstype,ind=indnum,x=x,y=y,n=n,occ=occnum))
		}
	}
	ggplot()+geom_tile(aes(x=x,y=y,fill=z),data=tiledata)+
		geom_point(data=alltrap,aes(x=x,y=y),size=3,pch=0,stroke=1.5)+
		geom_path(data=traj,aes(x=x,y=y),alpha=0.5)+
		geom_point(data=inddetect,aes(x=x,y=y,size=n),col="blue")+
		scale_fill_gradientn(colours = rev(terrain.colors(20))[-1],trans=tiletrans)+
		labs(fill = covname)+
		theme_cowplot()
}


generate_init<-function(secrad){
	envmodel<-secrad$envmodel
	indmodel<-secrad$indmodel
	occmodel<-secrad$occmodel
	
	modellist<-lapply(envmodel,"[",-2)
	for(i in 1:length(modellist)){
		names(modellist)[i]<-as.character(envmodel[[i]][2])
	}
	dname<-paste0("dens_",c("0",attr(terms(modellist[["D"]]),"term.labels")))
	cname<-paste0("conn_",c("0",attr(terms(modellist[["C"]]),"term.labels")))
	acovname<-attr(terms(modellist[["A"]]),"term.labels")
	if(length(acovname)!=0){
		aname<-paste0("adv_",acovname)
	}else{
		aname<-character(0)
	}
	
	g0name<-paste0("g0_",1:length(secrad$secrdata$obs))
	
	aindname<-character(0)
	g0indname<-character(0)
	mixname<-character(0)

	if(any(indmodel)){
		indcov<-secrad$secrdata$ind_cov
		indlevel<-levels(factor(indcov))[-1]
		mixname<-paste0("mix_",indlevel)
		
		if(indmodel["A"]){
			aindname<-paste0("adv_ind_",indlevel)
		}
		if(indmodel["g0"]){
			g0indname<-paste0("g0_ind_",indlevel)
		}
	}

	if(occmodel["A"]){
		occcov<-secrad$secrdata$occ_cov
		occlevel<-levels(factor(occcov))[-1]
		aoccname<-paste0("adv_occ_",occlevel)
	}else{
		aoccname<-character(0)
	}
	if(occmodel["g0"]){
		occcov<-secrad$secrdata$occ_cov
		occlevel<-levels(factor(occcov))[-1]
		g0occname<-paste0("g0_occ_",occlevel)
	}else{
		g0occname<-character(0)
	}
	nameall<-c(dname,cname,aname,g0name,aindname,g0indname,aoccname,g0occname,mixname)
	npar<-length(nameall)
	res<-rep(0,npar)
	names(res)<-nameall
	return(res)
}
		
###output posterior activity center
actcent<-function(par,secrad_obj,Nhat=FALSE,Utility=FALSE){
			secrdata<-secrad_obj$secrdata
			envmodel<-secrad_obj$envmodel
			indmodel<-secrad_obj$indmodel
			occmodel<-secrad_obj$occmodel
			####prepare data for estimation
			ncell<-secrdata$ncell
			nocc<-secrdata$nocc
			nind<-secrdata$nind
			
			#design matrices of grid
			gridcov<-secrdata$grid_cov
			if(is.null(gridcov)){
				if(any(is.na(gridcov))){
					stop("NA must not be included in covariates.\n")
				}
			}
			modellist<-lapply(envmodel,"[",-2)
			for(i in 1:length(modellist)){
				names(modellist)[i]<-as.character(envmodel[[i]][2])
			}
			denv.mat<-model.matrix(modellist$D,data=gridcov)
			denv.offset<-model.offset(model.frame(modellist$D,data=gridcov))
			if(is.null(denv.offset)){
				denv.offset<-rep(0,ncell)
			}
			
			cenv.mat<-model.matrix(modellist$C,data=gridcov)
			cenv.offset<-model.offset(model.frame(modellist$C,data=gridcov))
			if(is.null(cenv.offset)){
				cenv.offset<-rep(0,ncell)
			}
			
			woaenv<-F
			aenv.mat<-model.matrix(modellist$A,data=gridcov)
			if(length(aenv.mat)==0){
				aenv.mat<-matrix(0,nrow=ncell,ncol=1)
				woaenv<-T
			}else if(is.element("(Intercept)",colnames(aenv.mat))){
				stop("Intercept cannot be included in the advection submodel.")
			}
			aenv.offset<-model.offset(model.frame(modellist$A,data=gridcov))
			if(is.null(aenv.offset)){
				aenv.offset<-rep(0,ncell)
			}
			
			#individual covariate
			if(any(indmodel)){
				indcov<-secrdata$ind_cov
			}else{
				indcov<-rep(0,nind)
			}
			
			#occasion covariate
			if(any(occmodel)){
				occcov<-secrdata$occ_cov
			}else{
				occcov<-rep(0,nocc)
			}
				
			###partitioning parameters
			#density
			ndpar<-ncol(denv.mat)
			dpar<-par[1:ndpar]
			#connectivity
			ncpar<-ncol(cenv.mat)
			cpar<-par[(1:ncpar)+ndpar]
			#advection
			napar<-ifelse(woaenv,0,ncol(aenv.mat))
			if(napar!=0){
				apar<-apar2<-par[(1:napar)+ndpar+ncpar]
			}else{
				apar<-numeric(0)
				apar2<-0
			}
			#g0
			nobs<-length(secrdata$obs)
			obstype<-unlist(lapply(secrdata$obs,"[[","type"))
			npar_type<-c(1,1);names(npar_type)<-c("binom","poisson")
			ng0par_obs<-npar_type[obstype]
			g0par<-vector("list",nobs)
			ng0par<-0
			for(i in 1:nobs){
				g0par[[i]]<-par[(1:ng0par_obs[i])+ng0par+ndpar+ncpar+napar]
				ng0par<-ng0par+ng0par_obs[i]
			}
			#occasion
			noccgrp<-length(unique(occcov))
			noccpar<-sum(occmodel)*(noccgrp-1)
			if(noccpar!=0){
				occpar<-par[(1:noccpar)+ndpar+ncpar+napar+ng0par]
				names(occpar)<-paste0(rep(names(occmodel)[occmodel],each=noccgrp-1),"_",rep(levels(factor(occcov))[-1],sum(occmodel)))
			}else{
				occpar<-numeric(0)
			}
			
			nindgrp<-length(unique(indcov))	
			nindpar<-sum(indmodel)*(nindgrp-1)
			nmixpar<-nindgrp-1
			if(nindpar!=0){
				indpar<-par[(1:nindpar)+ndpar+ncpar+napar+ng0par+noccpar]
				names(indpar)<-paste0(rep(names(indmodel)[indmodel],each=nindgrp-1),"_",rep(levels(factor(indcov))[-1],sum(indmodel)))
				mixpar<-par[(1:nmixpar)+nindpar+ndpar+ncpar+napar+ng0par+noccpar]
			}else{
				indpar<-numeric(0)
				mixpar<-numeric(0)
			}
		
			if(!is.null(secrdata$auxiliary)){
				aux<-secrdata$auxiliary
			}else{
				aux<-NULL
			}
			#survival matrix(not implemented)
			srv<-matrix(1,nrow=secrdata$nind,ncol=nocc)
		#	if(!is.null(aux)){
		#		auxtype<-unlist(lapply(secrdata$auxiliary,"[","type"))
		#		if(is.element("removal",auxtype)){
		#			numremoval<-which(auxtype=="removal")
		#			for(i in 1:length(numremoval)){
		#				whichremove<-which(secrdata$auxiliary[[numremoval[i]]]$detect==1,arr.ind=T)
		#				doublet<-cbind(ind=whichremove[,2],occ=secrdata$auxiliary[[numremoval[i]]]$effort_occ[whichremove[,1]])
		#				for(j in 1:nrow(doublet)){
		#					if(doublet[j,2]<nocc){
		#						srv[doublet[j,1],(doublet[j,2]+1):nocc]<-0
		#					}
		#				}
		#			}
		#		}
		#	}
			
			#
			#log mean density
			lnD<-c(denv.mat%*%dpar+denv.offset)
			
			#individual covariate of advection
			grepA<-grep("^A",names(indpar))
		
			if(length(grepA)!=0){
				apar_ind<-c(0,indpar[grepA])
			}else{
				apar_ind<-0
			}
			
			#individual covariate of detection
			grepg0<-grep("^g0",names(indpar))
			if(length(grepg0)!=0){
				g0par_ind<-c(0,indpar[grepg0])
			}else{
				g0par_ind<-0
			}
			ng0par_ind<-length(g0par_ind)
		
			#occasion covariate of advection
			grepAocc<-grep("^A",names(occpar))
		
			if(length(grepAocc)!=0){
				apar_occ<-c(0,occpar[grepAocc])
			}else{
				apar_occ<-0
			}
			
			#occasion covariate of detection
			grepg0occ<-grep("^g0",names(occpar))
			if(length(grepg0occ)!=0){
				g0par_occ<-c(0,occpar[grepg0occ])
			}else{
				g0par_occ<-0
			}
			ng0par_occ<-length(g0par_occ)
		
			
			#list of g0 coef
			g0list<-vector("list",ng0par)
			for(i in 1:nobs){
				g0list[[i]]<-outer(matrix(1,ng0par_ind,ng0par_occ),g0par[[i]])
				g0list[[i]][,,1]<-g0list[[i]][,,1]+outer(g0par_ind,rep(1,ng0par_occ))+outer(rep(1,ng0par_ind),g0par_occ)
			}
		
			#advection-diffusion model
			logp_arr<-array(NA,dim=c(ncell,ncell,length(apar_ind),length(apar_occ)))
			for(i in 1:length(apar_ind)){
				for(j in 1:length(apar_occ)){
					adpar<-c(cpar,apar,apar_ind[i],apar_occ[j])
					logc<-c(cenv.mat%*%matrix(cpar))+cenv.offset
					loga<-c(aenv.mat%*%matrix(apar2))+aenv.offset+apar_ind[i]+apar_occ[j]
					logp_arr[,,i,j]<-advdiff.eigen(mucoords=secrdata$coords,
									gridcoords=secrdata$coords,
									logc=logc,
									loga=loga,
									neighmat=secrdata$neighmat,
									resolution=secrdata$resolution,
									logprob=T)
					}
			}
				
			#mean number of detection (or cloglog of prob for binimial)
			loglambdalist<-vector("list",nobs)
			nindgrp_lambda<-pmax(length(apar_ind),length(g0par_ind))
			for(i in 1:nobs){
				loglambdalist[[i]]<-array(NA,dim=c(ncell,length(secrdata$obs[[i]]$effort),nindgrp_lambda))
				A_effort_occtype<- ifelse(rep(occmodel["A"],length(secrdata$obs[[i]]$effort_occ)),
								occcov[secrdata$obs[[i]]$effort_occ]+1,
								rep(1,length(secrdata$obs[[i]]$effort_occ)))
				g0_effort_occtype<-ifelse(rep(occmodel["g0"],length(secrdata$obs[[i]]$effort_occ)),
								occcov[secrdata$obs[[i]]$effort_occ]+1,
								rep(1,length(secrdata$obs[[i]]$effort_occ)))
				for(j in 1:nindgrp_lambda){
					if(indmodel["A"]){
						logp<-logp_arr[,,j,,drop=F]
					}else{
						logp<-logp_arr[,,1,,drop=F]
					}
					g0list_indocc<-g0list[[i]][j,g0_effort_occtype,,drop=F]
					loglambdalist[[i]][,,j]<-loglambda.ad(type=obstype[i],g0list_indocc,secrdata$obs[[i]]$effort,secrdata$obs[[i]]$effort_loc,A_effort_occtype,logp,logprob=T)
				}
			}
			
			####calculation of pdot
			#log probability of non-detection of each individual group and occasion group
			logpzerolist<-vector("list",nobs)
			for(i in 1:nobs){
				logpzerolist[[i]]<-array(NA,dim=c(ncell,length(secrdata$obs[[i]]$effort),nindgrp))
				for(j in 1:nindgrp){
						logpzerolist[[i]][,,j]<-logpzero.ad(type=obstype[i],loglambdalist[[i]][,,j],secrdata$obs[[i]]$effort)
				}
			}
			
			#log 1-pdot for each occasion conditional on X
			logpzero_occ_x_arr<-array(NA,dim=c(ncell,nocc,nindgrp))
			for(i in 1:nocc){
				logpzerolist_occ<-vector("list",nobs)
				for(j in 1:nobs){
					cond<-secrdata$obs[[j]]$effort_occ==i
					logpzerolist_occ[[j]]<-logpzerolist[[j]][,cond,,drop=F]
				}
				
				for(k in 1:nindgrp){
					logpzerolist_occ_grp<-lapply(logpzerolist_occ,"[",,,k)
					logpzero_occ_x_arr[,i,k]<-logpzero.x(logpzerolist_occ_grp)
				}
			}
			
		
			#ln(pdot) for each individual
			logpzero_occ_x_ind<-logpzero_occ_x_arr[,,indcov-min(indcov)+1,drop=F]
			logpzero_x_ind<-(logpzero_occ_x_ind*outer(rep(1,ncell),t(srv)))%>%
						aperm(perm=c(1,3,2))%>%
						rowSums(dim=2)
			logpdot_x_ind<-log1m_exp_mat(logpzero_x_ind)			
		
			#lnD+ln(area)+ln(pdot) for each individual
			logd_area_pdot_x_ind<-outer(lnD+log(secrdata$area),rep(1,nind))+logpdot_x_ind
			
			#logsumexp(lnD+ln(area)+ln(pdot)) for each individual
			logd_area_pdot_ind<-colsum_lp(logd_area_pdot_x_ind)
			
			#####likelihood of individual capture history conditional on detection at least once
			#multinomial coefficient
			detectlist<-lapply(secrdata$obs,"[[","detect")
			detect_append<-detectlist[[1]]
			if(nobs>1){
				for(i in 2:nobs){
					detect_append<-rbind(detect_append,detectlist[[i]])
				}
			}
			detect_pattern<-apply(detect_append,2,paste,collapse="_")
			nc<-table(detect_pattern)
			names(nc)<-NULL
			logmulticoef<-lgamma(nind+1)-sum(lgamma(nc+1))
			
			#log prob of capture history conditional on X
			logpcapthist_arr<-array(0,dim=c(ncell,nind,nobs))
			for(i in 1:nobs){
				chpar<-c(cpar,apar,apar_ind,apar_occ,c(g0list[[i]]))
				
				logpcapthist_arr[,,i]<-logpcapthist(type=obstype[i],
						detect=secrdata$obs[[i]]$detect,
						loglambda=loglambdalist[[i]],
						effort=secrdata$obs[[i]]$effort,
						effort_occ=secrdata$obs[[i]]$effort_occ,
						srv=srv,
						ind_cov=indcov)
			}
			
			logpcapthist_mat<-rowSums(logpcapthist_arr,dim=2)
			logpost_mu<-logpcapthist_mat-outer(rep(1,ncell),apply(logpcapthist_mat,2,sum_lp))
			map_mu<-apply(logpost_mu,2,which.max)
			
			map_hr<-array(NA,dim=c(nind,ncell,length(apar_occ)))
			for(i in 1:nind){
				map_hr[i,,]<-logp_arr[map_mu[i],,indcov[i]+1,]
			}
			
			logpcapthist_0<-array(0,dim=c(ncell,length(apar_ind),nobs))
			for(j in 1:length(apar_ind)){
				for(i in 1:nobs){
					chpar<-c(cpar,apar,apar_ind,apar_occ,c(g0list[[i]]))
				
					logpcapthist_0[,j,i]<-logpcapthist(type=obstype[i],
							detect=matrix(0,length(secrdata$obs[[i]]$effort),1),
							loglambda=loglambdalist[[i]],
							effort=secrdata$obs[[i]]$effort,
							effort_occ=secrdata$obs[[i]]$effort_occ,
							srv=srv,
							ind_cov=unique(indcov)[j])
				}
			}
			
			logpcapthist_0_mat<-rowSums(logpcapthist_0,dim=2)
			logpost_0<-logpcapthist_0_mat-outer(rep(1,ncell),apply(logpcapthist_0_mat,2,sum_lp))


			Nhat<-sum(exp(lnD)*secrdata$area)
			Nud<-Nhat-nind
			if(Nhat){
				Nhat_sp<-matrix(0,ncell,length(apar_ind))
				for(i in 1:nind){
					Nhat_sp[,indcov[i]+1]<-Nhat_sp[,indcov[i]+1]+exp(logpost_mu[,i])
				}
				pmix<-exp(c(0,mixpar))/sum(exp(c(0,mixpar)))
				for(i in 1:length(apar_ind)){
					Nhat_sp[,i]<-Nhat_sp[,i]+exp(logpost_0[,i])*Nud*pmix[i]
				}
			}else{
				Nhat_sp<-NULL
			}


			if(Utility){
				Utility_sp<-array(0,dim=c(ncell,length(apar_ind),length(apar_occ)))
				for(i in 1:nind){
					Utility_sp[,indcov[i]+1,]<-Utility_sp[,indcov[i]+1,]+apply(exp(logp_arr[,,indcov[i]+1,]+outer(logpost_mu[,i],matrix(1,ncell,length(apar_occ)))),c(2,3),sum)
				}
				for(i in 1:length(apar_ind)){
					Utility_sp[,i,]<-Utility_sp[,i,]+apply(exp(logp_arr[,,i,]+outer(logpost_0[,i],matrix(1,ncell,length(apar_occ)))),c(2,3),sum)*Nud*pmix[i]
				}
			}else{
				Utility_sp<-NULL
			}
			res<-list(logp_arr=logp_arr,logpost_mu=logpost_mu,map_hr=map_hr,Nhat=Nhat_sp,Utility=Utility_sp)
			return(res)
}



#ggplot()+geom_tile(data=coords%>%bind_cols(N=Nhat_sp),aes(x=x,y=y,fill=log(Utility_sp[,1,1])),width=secrdata$resolution["x"], height=secrdata$resolution["y"])+
#		scale_fill_gradientn(colours = rev(terrain.colors(20))[-1])


