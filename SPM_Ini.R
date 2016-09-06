####Generate the vector for the normalization of the prior model 


Normalize_dose<-function(dose,mat_prior_model,c,epsilon)
{
  if( epsilon==0 ){
    
    FInc<-mat_prior_model[,dose]
    norunif<-(alpha^(1-dose))*((1-alpha)^(dose-sizeD))
    norbeta<-prod(beta(FInc*c+1,(1-FInc)*c+1))*dbeta(alpha,alpha*c+1,(1-alpha)*c+1)
    coef1<-1
    coef2<-1
    if(dose>1){coef1<-prod(pbeta(alpha,FInc[1:(dose-1)]*c+1,(1-FInc[1:(dose-1)])*c+1))}
    if(dose<6){coef2<-prod(1-pbeta(alpha,FInc[(dose+1):sizeD]*c+1,(1-FInc[(dose+1):sizeD])*c+1))}
    resDose<-norunif*norbeta*coef1*coef2
  } else {
    FInc<-mat_prior_model[,dose]
    norunif<-((alpha-epsilon)^(1-dose))*((1-alpha-epsilon)^(dose-sizeD))*(2*epsilon)
    norbeta<-prod(beta(FInc*c+1,(1-FInc)*c+1))
    coef1<-1
    coef2<-1
    if(dose>1){coef1<-prod(pbeta(alpha-epsilon,FInc[1:(dose-1)]*c+1,(1-FInc[1:(dose-1)])*c+1))}
    if(dose<6){coef2<-prod(1-pbeta(alpha+epsilon,FInc[(dose+1):sizeD]*c+1,(1-FInc[(dose+1):sizeD])*c+1))}
    coef3<-pbeta(alpha+epsilon,FInc[dose]*c+1,(1-FInc[dose])*c+1)-pbeta(alpha-epsilon,FInc[dose]*c+1,(1-FInc[dose])*c+1)
    resDose<-norunif*norbeta*coef1*coef2*coef3
  }
  
  
  return(resDose)
  
}

Ini<-function(mat_prior_model,c,epsilon)
{
  Vec_norm<-sapply(D,Normalize_dose,mat_prior_model=mat_prior_model,c=c,epsilon=epsilon)
  return(Vec_norm)
}

#SPM_Ini(mat_prior_model,c,epsilon)



