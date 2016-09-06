#Generate a trial

Trial<-function(sizeTr,T,c,mat_prior_model,Vec_norm, Pi,epsilon,cohort,skip,caution,dose_exclusion)
{
  
  #removable
  Sum_Obs<-matrix(rep(0, 2*sizeD),2,sizeD)
  Seq_Obs<-matrix(rep(0, 2*sizeTr),2,sizeTr)
  
  #Variable
  matPi_n<-matrix(rep(0,sizeD*sizeTr),sizeTr,sizeD)
  Caution<-rep(caution,sizeD)
  dcur<-1
  i<-1
  max<-sizeD
  
  while(i<(sizeTr+1) && max>0)
      {
        #i=1
        ind<-dcur
        inf<-rbinom(1,1,T[ind])
        Seq_Obs[,i]<-c(ind,inf)
        Sum_Obs[,ind]<-Sum_Obs[,ind]+c(inf,1-inf)
        Caution[dcur]<-Caution[dcur]-1
        
        if((i%%cohort)==0 && Caution[dcur]<1){
          Pi_n<-A_Posteriori(Sum_Obs,mat_prior_model,c,Vec_norm,Pi,epsilon)
          matPi_n[i,]<-Pi_n
          temp<-which.max(Pi_n)
          if(skip==1){dcur<-temp}
          if(skip==0){dcur<-dcur+(temp-dcur>0)-(temp-dcur<0)}
        }
        if (dose_exclusion>0 && (i%%cohort)==0){
          elim<-1-pbeta(alpha,Sum_Obs[1,]+1,Sum_Obs[2,]+1)
          if(sum(elim>dose_exclusion)>0){max<-which.max(elim>0.95)-1}
          if(dcur>max && sum(Sum_Obs[,max+1])>2){dcur<-max}
          if(max==0 && sum(Sum_Obs[,max+1])<3){max<-1;dcur<-max}
        }
        i<-i+1
        
        #Pi_n
        #seq
        #Sum_Obs
      }
    
    
  
  
  nbrpat<-i-1
  mtd<-dcur
  Observations<-list(Sum_Obs=Sum_Obs,Seq_Obs=Seq_Obs)
  return(list(Observations=Observations,Pi_n=Pi_n,mtd=mtd,nbrpat=nbrpat,matPi_n=matPi_n))
  
}

#Trial(sizeTr,T,c,mat_prior_model,Vec_norm, Pi,epsilon,cohort,skip,caution,dose_exclusion)