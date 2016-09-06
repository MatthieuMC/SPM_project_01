####Generate the a posteriori of Pi according to the observations (Sum_Obs) 

A_Posteriori_dose<-function(dose, Sum_Obs, mat_prior_model, c, Vec_norm, epsilon)
{
  if(epsilon==0){
    FInc<-mat_prior_model[,dose]
    norunif<-((alpha-epsilon)^(1-dose))*((1-alpha-epsilon)^(dose-sizeD))
    norbeta<-prod(beta(FInc*c+Sum_Obs[1,]+1,(1-FInc)*c+Sum_Obs[2,]+1))*dbeta(alpha,alpha*c+Sum_Obs[1,dose]+1,(1-alpha)*c+Sum_Obs[2,dose]+1)
    coef1<-1
    coef2<-1
    if(dose>1){coef1<-prod(pbeta(alpha-epsilon,FInc[1:(dose-1)]*c+Sum_Obs[1,1:(dose-1)]+1,(1-FInc[1:(dose-1)])*c+Sum_Obs[2,1:(dose-1)]+1))}
    if(dose<6){coef2<-prod(1-pbeta(alpha+epsilon,FInc[(dose+1):sizeD]*c+Sum_Obs[1,(dose+1):sizeD]+1,(1-FInc[(dose+1):sizeD])*c+Sum_Obs[2,(dose+1):sizeD]+1))}
    resDose<-norunif*norbeta*coef1*coef2*(1/Vec_norm[dose])
  } else {
    FInc<-mat_prior_model[,dose]
    sizeD<-length(FInc)
    norunif<-((alpha-epsilon)^(1-dose))*((1-alpha-epsilon)^(dose-sizeD))*(2*epsilon)
    norbeta<-prod(beta(FInc*c+Sum_Obs[1,]+1,(1-FInc)*c+Sum_Obs[2,]+1))
    coef1<-1
    coef2<-1
    if(dose>1){coef1<-prod(pbeta(alpha-epsilon,FInc[1:(dose-1)]*c+Sum_Obs[1,1:(dose-1)]+1,(1-FInc[1:(dose-1)])*c+Sum_Obs[2,1:(dose-1)]+1))}
    if(dose<6){coef2<-prod(1-pbeta(alpha+epsilon,FInc[(dose+1):sizeD]*c+Sum_Obs[1,(dose+1):sizeD]+1,(1-FInc[(dose+1):sizeD])*c+Sum_Obs[2,(dose+1):sizeD]+1))}
    coef3<-pbeta(alpha+epsilon,FInc[dose]*c+Sum_Obs[1,dose]+1,(1-FInc[dose])*c+Sum_Obs[2,dose]+1)-pbeta(alpha-epsilon,FInc[dose]*c+Sum_Obs[1,dose]+1,(1-FInc[dose])*c+Sum_Obs[2,dose]+1)
    resDose<-norunif*norbeta*coef1*coef2*coef3*(1/Vec_norm[dose])
    
  }
  return(resDose)
}


A_Posteriori<-function( Sum_Obs,mat_prior_model, c, Vec_norm,Pi,epsilon)
{
  Pi_n<-sapply(D,A_Posteriori_dose,Sum_Obs,mat_prior_model=mat_prior_model,c=c,Vec_norm=Vec_norm,epsilon=epsilon)
  Pi_n<-Pi_n*Pi/sum(Pi_n*Pi)
  return(Pi_n)
}

#A_Posteriori( Sum_Obs, mat_prior_model, c, Vec_norm, Pi, epsilon)

