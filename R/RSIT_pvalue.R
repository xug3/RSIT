#' Calculate  P-values 
#'
#' This function tests a SNPs for a given SNP set for a given rsit estimated null model.
#'
#' @param TYall The output of function "RSIT_est()"
#' @param G The genotype matrix, an m*q matrix where m is the number of subjects and q is the total number genetic variants. 
#' @param impute.method choose the iputation method when there is missing genotype. Optional options are: 'random', 'fixed' or 'bestguess'.
#' 
#' @return This function returns a vector of the pvalue of RSIT-B,RSIT-S,RSIT-V,RSIT-A,RSIT-O, and RSIT-E.
#' 
#' @export



RSIT_pvalue=function(TYall,G,impute.method='fixed'){
  G[G==9]<-NA
  N_MISS<-sum(is.na(G))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G<-Impute(G,'fixed')
  }
  
  var.G<-t(t(G)-colMeans(G));
  SNP.list<-colnames(G)
  MAF<-colMeans(as.matrix(G[,colMeans(var.G^2)!=0]), na.rm=T)/2;
  minormaf=which(MAF>0.5)
  if(length(minormaf)>0) G[,minormaf]=2-G[,minormaf]
  MAF<-colMeans(as.matrix(G[,colMeans(var.G^2)!=0]), na.rm=T)/2;
  G<-as.matrix(var.G[,colMeans(var.G^2)!=0]); 
  #G<-as.matrix(G[,colMeans(var.G^2)!=0]); 
  SNP.name<-SNP.list[colMeans(var.G^2)!=0]
  
  weights = rep(0,ncol(G));
  index.common = MAF>=0.05
  index.rare = MAF<0.05
  scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
  weights[index.common] = scaler*dbeta(MAF[index.common],0.5,0.5)
  weights[index.rare] = dbeta(MAF[index.rare],1,25)
  
  
  WG<-t(t(G)*weights)
  Var_WG=cov(WG)
  TSRIT=sum((t(WG)%*%TYall$TY)^2)
  Sigma=Var_WG*TYall$varTY
  lambdaS=eigen(Sigma,symmetric = TRUE,only.values = TRUE)$values
  RSIT_S = generalchisq(lambdaS, TSRIT)
  
  
  wG=G%*%weights
  var_g<-var(wG)
  score=t(TYall$TY)%*%wG
  std.retro=sqrt(var_g*TYall$varTY)
  score.retro<-c(score/std.retro)
  RSIT_B<-pchisq((score.retro)^2,df=1, lower.tail=F) 
  scoreO=t(WG)%*%TYall$TY
  RSIT_O= Get.p.eigen.skato2(scoreO,Sigma,rho=(0:10)/10)$p    
  
  WG2=WG-rowSums(WG)%*%t(colSums(Var_WG))/sum(Var_WG)
  TSRIT2=sum((t(WG2)%*%TYall$TY)^2)
  Sigma2=(Var_WG-rowSums(Var_WG)%*%t(colSums(Var_WG))/sum(Var_WG))*TYall$varTY
  lambdaS2=eigen(Sigma2,symmetric = TRUE,only.values = TRUE)$values
  RSIT_S2 = generalchisq(lambdaS2, TSRIT2)
  Fisher=-2*log(RSIT_B*RSIT_S2)
  RSIT_E<-pchisq(Fisher,df=4, lower.tail=F) 
  pV=RGET_pvalue(TYall,G)
  weightV=weights*(MAF*(1-MAF))^0.5
  RSIT_V=aca.pvalue(pV,weightV)
  RSIT_A=aca.pvalue(c(RSIT_B,RSIT_S,RSIT_V),rep(1/3,3))
  pvalue=c(RSIT_B,RSIT_S,RSIT_V,RSIT_A,RSIT_O,RSIT_E)
  return(pvalue)
}



RGET_pvalue=function(TYall,G,impute.method='fixed'){
  G[G==9]<-NA
  N_MISS<-sum(is.na(G))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G<-Impute(G,'fixed')
  }
  
  var.G<-t(t(G)-colMeans(G));
  SNP.list<-colnames(G)
  MAF<-colMeans(as.matrix(G[,colMeans(var.G^2)!=0]), na.rm=T)/2;
  minormaf=which(MAF>0.5)
  if(length(minormaf)>0) G[,minormaf]=2-G[,minormaf]
  MAF<-colMeans(as.matrix(G[,colMeans(var.G^2)!=0]), na.rm=T)/2;
  #G<-as.matrix(G[,colMeans(var.G^2)!=0]); 
  #SNP.name<-SNP.list[colMeans(var.G^2)!=0]
  
  var_g<-apply(G, 2, var)
  score=t(TYall$TY)%*%G
  std.retro=sqrt(var_g*TYall$varTY)
  score.retro<-c(score/std.retro)
  pval.retro<-pchisq((score.retro)^2,df=1, lower.tail=F)     
  return(pval.retro)
}


Get.p.eigen.skato2 <- function(score, V, rho, method = "saddle") {
  n.r <- length(rho);lambdas <- vector("list", n.r);pval <- qval <- rep(NA, n.r)
  n.p <- length(score); U  <- score;
  Q <- (1-rho)*sum(U^2)+rho*sum(U)^2
  Burden.score <- Burden.var <- Burden.pval <- SKAT.pval <- NA
  for(i in 1:n.r) {
    if(rho[i]==1) {
      Burden.score <- sum(U)
      Burden.var <- sum(V)
      Burden.pval <- pchisq(Burden.score^2/Burden.var, df=1, lower.tail=FALSE)
      lambdas[[i]] <- Burden.var
      pval[i] <- Burden.pval
      next
    }
    if(rho[i]!=0) {
      R.M <- matrix(rho[i], n.p, n.p)
      diag(R.M) <- 1
      R.M.chol <- t(chol(R.M, pivot = TRUE))
      V.temp <- crossprod(R.M.chol, crossprod(V, R.M.chol))
    } else V.temp <- V
    lambda <- eigen(V.temp, only.values = TRUE, symmetric=TRUE)$values
    lambdas[[i]] <- lambda[lambda > 0]
    pval[i] <- generalchisq(lambdas[[i]],Q[i])
    if(rho[i]==0) SKAT.pval <- pval[i]
  }
  minp <- min(pval)
  for(i in 1:n.r) {
    df <- sum(lambdas[[i]]^2)^2/sum(lambdas[[i]]^4)
    qval[i] <- (qchisq(minp, df, lower.tail = FALSE)-df)/sqrt(2*df)*sqrt(2*sum(lambdas[[i]]^2))+sum(lambdas[[i]])
  }
  ZMZ <- tcrossprod(rowSums(V))/sum(V)
  V.temp <- V - ZMZ
  lambda <- eigen(V.temp, only.values = TRUE, symmetric = TRUE)$values
  lambda <- lambda[lambda > 0]
  muq <- sum(lambda)
  varq <- sum(lambda^2) * 2 + sum(ZMZ * V.temp) * 4
  df <- sum(lambda^2)^2/sum(lambda^4)
  tau <- rho * sum(V) + sum(V %*% V)/sum(V) * (1 - rho)
  re <- tryCatch({
    integrate(function(x){
      t1 <- tau %x% t(x)
      re<-pchisq((apply((qval - t1)/(1-rho),2,min) - muq)/sqrt(varq)*sqrt(2*df) + df, df=df) * dchisq(x,df=1)
      return(re)
    }, lower = 0, upper = 40, subdivisions = 2000, abs.tol = 10^-25)
  }, error=function(e) NA)
  return(list(p = min(1-re[[1]], minp*n.r), minp = minp, minp.rho = rho[which.min(pval)],
              Burden.score=Burden.score, Burden.var=Burden.var, Burden.pval=Burden.pval,
              SKAT.pval=SKAT.pval))
}

aca.pvalue  <- function(pvalue, weights){
  rm.idx = which(is.na(pvalue)|pvalue>(1-10^-5));
  if(length(rm.idx)>0){
    pvalue = pvalue[-rm.idx]; weights = weights[-rm.idx]}
  T_aca <- c(weights^2)%*%tan((0.5-pvalue)*pi)
  # cat("T_aca: ", T_aca, '\n')
  p.aca = 0.5-atan(T_aca/sum(weights^2))/pi
}




generalchisq=function(lambda,Q){ #function of evaluate pvalue of mixture of chisquare
  muq=sum(lambda)
  sigmaq=sqrt(2*sum(lambda^2))
  s1=sum(lambda^3)/sum(lambda^2)^1.5
  s2=sum(lambda^4)/sum(lambda^2)^2
  if((s1^2)>s2){
    a=1/(s1-sqrt(s1^2-s2))
    delta=s1*a^3-a^2
    l=a^2-2*delta
  }else{
    delta=0
    l=sum(lambda^2)^3/sum(lambda^3)^2
  }
  pchisq((sum(Q)-muq)/sigmaq*sqrt(2*(l+2*delta))+l+delta,df=l, ncp=delta,lower.tail=FALSE)}


Impute<-function(Z, impute.method){
  if(is.vector(Z)){
    if(impute.method =="random"){
      IDX<-which(is.na(Z))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX])/2
        Z[IDX]<-rbinom(length(IDX),2,maf1)
      }
    } else if(impute.method =="fixed"){
      IDX<-which(is.na(Z))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX])/2
        Z[IDX]<-2 * maf1
      }
    } else if(impute.method =="bestguess") {
      IDX<-which(is.na(Z))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX])/2
        Z[IDX]<-round(2 * maf1)
      }
    } else {
      stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
    }
  } else{
    p<-dim(Z)[2]
    if(impute.method =="random"){
      for(i in 1:p){
        IDX<-which(is.na(Z[,i]))
        if(length(IDX) > 0){
          maf1<-mean(Z[-IDX,i])/2
          Z[IDX,i]<-rbinom(length(IDX),2,maf1)
        }
      }
    } else if(impute.method =="fixed"){
      for(i in 1:p){
        IDX<-which(is.na(Z[,i]))
        if(length(IDX) > 0){
          maf1<-mean(Z[-IDX,i])/2
          Z[IDX,i]<-2 * maf1
        }
      }
    } else if(impute.method =="bestguess") {
      for(i in 1:p){
        IDX<-which(is.na(Z[,i]))
        if(length(IDX) > 0){
          maf1<-mean(Z[-IDX,i])/2
          Z[IDX,i]<-round(2 * maf1)
        }
      }
    } else {
      stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
    }
  }
  return(as.matrix(Z))
}
