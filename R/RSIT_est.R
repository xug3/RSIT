#'  NULL model estimation
#'
#' This function estimate the parameters and residuals for the NULL model in RSIT
#'
#' @param Y phenotype vector
#' @param X Covariate matrix denoting the covariate variables
#' @param E environment vector
#' @param family String, the phenotype model, two optional values: 'Guassian', 'binomial'
#' @param GRM takes m-by-m genetic correlation matrix or kinship matrix.
#' @return This function returns a list object with model parameters and residuals of the NULL model
#' @export

RSIT_est=function(Y,X,E,family="Guassian", GRM = NULL){ 
    if(!is.null(GRM)){
        dGRM=sqrt(diag(GRM))
        for(i in 1:dim(GRM)[1]) {
             GRM[,i]=GRM[,i]/dGRM[i]
             GRM[i,]=GRM[i,]/dGRM[i]
        }
    }
    if(family=="Guassian"){
      mod=lm(Y~X+E)
      }else if(family=="binomial"){ 
          mod=glm(Y~X+E,family=binomial(link="logit"))
          } else {stop("error distribution, please choose Guassian or binomial")}
     res=Y-mod$fit
     if(!is.null(GRM)){
            phires=GRM%*%res
     } else{
         phires=res
     }
    
    TY=E*res-res*sum(E*res*phires)/sum(res*phires)

    if(!is.null(GRM)){
     phiTY=GRM%*%TY
     } else{
     phiTY=TY
     }
     varTY=sum(TY*phiTY)
     TYall=list(TY=TY,varTY=varTY)
     return(TYall)
}
