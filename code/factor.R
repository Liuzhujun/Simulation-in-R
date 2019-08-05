get_residuals=function(Y,X,beta){
  res=Y
  for(k in 1:K){
    res=res-beta[k]*drop(X[,,k])
  }
  return(res)
}

principal_components=function(res,R){
  T=ncol(res)
  V=eigen(t(res)%*%res)$vectors;
  D=eigen(t(res)%*%res)$values;
  #e=cbind(t(V),D);
  #newe=e[order(e[,T+1]),]
  #newV=t(newe[,1:T])
  #f=newV[,1:R]
  f=V[,1:R]
  for(r in 1:R){
    f[,r]=f[,r]/sqrt(sum(f[,r]^2));
    if (mean(f[,r])<0){
      f[,r]=-f[,r]
    }
  }
  lambda=res%*%f
  w=list(lambda=lambda,f=f)
  return(w)
}

LS_obj=function(beta,Y,X,R){
  res=get_residuals(Y=Y,X=X,beta=beta);
  ev=sort((eigen(t(res)%*%res)$values));
  obj=sum(ev[1:(T-R)])/N/T;
  return(obj)
}

minimize_obj_method1=function(Y,X,R,st,precision_beta){
  N=nrow(Y);#????Y????????
  T=ncol(Y);
  K=2;
  SST=sum(diag(Y%*%t(Y)))/N/T;
  beta=st;  # beta??С????ֵ;
  beta_old=st+Inf;
  obj=Inf;
  diff_obj=-Inf;
  #cat("new iter: \n")
  while((max(abs(beta-beta_old))>precision_beta)&&(diff_obj<=SST*10^-10)){
    res=get_residuals(Y=Y,X=X,beta=beta)
    w=principal_components(res=res,R)
    lambda=w[[1]]
    f=w[[2]]
#print(  all.equal(res-lambda%*%t(f),res))
#print(beta)
    res=res-lambda%*%t(f);
    obj_old=obj;              #???ɵ?Ŀ?꺯??
    obj=sum(diag(t(res)%*%res))/N/T;  #??С????Ŀ?꺯??
    diff_obj=obj-obj_old;
 # print(diff_obj)
    if(diff_obj<=0){
      YY=as.vector(Y)
      XX=matrix(0,N*T,K)
      for(k in 1:K){
        xx=drop(X[,,k]);
      #  xx=xx-lambda %*% (t(lambda) %*% t(ginv(xx)));    #??lambdaͶӰ??x?????ĵط?
      #  xx=xx-t(f %*% (t(f) %*% ginv(xx)))
        xx=xx-lambda %*% (ginv(lambda) %*% ((xx)));    #??lambdaͶӰ??x?????ĵط?
        xx=xx-t(f %*% (ginv(f) %*% t(xx)))
        XX[,k]=as.vector(xx)  #ʹ XX ??ΪNTxK matrix
      }
      beta_old=beta;                 #???ɵ?beta
      beta=solve(t(XX)%*%XX)%*%t(XX)%*%YY;      #??????С???˹��?
     # print(beta)
     # cat(beta_old,diff_obj,beta,"\n")
          }
  }
  if(diff_obj<=0){
    ef=1;       #good solution found
  }else{
    ef=-1;
  }
  obj=LS_obj(beta=beta,Y=Y,X=X,R=R); #????Ŀ?꺯??
  h=list(beta=beta,obj=obj,ef=ef)
  return(h)
}

#?˺???
K0=function(roth,X,res){
  temp0=0
  for(i in 1:N){
    for(j in 1:N){
    xx1=outer(X[,,1][i,],X[,,1][j,],"-")
    xx2=outer(X[,,2][i,],X[,,2][j,],"-")
    k_1 = dnorm((xx1)/roth[1],0,1,log=FALSE)/roth[1]
    k_2 = dnorm((xx2)/roth[2],0,1,log=FALSE)/roth[2]
    kij= k_1*k_2
    temp0 = temp0 + t(res[i,])%*%kij%*%res[j,]
    }
  }
  return(temp0)
}


kii_kernel=function(roth,X,res){
  temp1=0;
  for(i in 1:N){
    xx1=outer(X[,,1][i,],X[,,1][i,],"-")
    xx2=outer(X[,,2][i,],X[,,2][i,],"-")
    k_value1 = dnorm((xx1)/roth[1],0,1,log=FALSE)/roth[1];
    k_value2 = dnorm((xx2)/roth[2],0,1,log=FALSE)/roth[2];
    kii = k_value1*k_value2
    temp1 = temp1 + t(res[i,])%*%kii%*%res[i,]
  }
  return(temp1)
}

Kh_R=function(roth,X,res){
  temp=0;
  for(i in 1:N){
    for(j in 1:N){
      if(i!=j){
      xx1=outer(X[,,1][i,],X[,,1][j,],"-")
      xx2=outer(X[,,2][i,],X[,,2][j,],"-")
      ka1 = dnorm((xx1)/roth[1],0,1,log=FALSE)/roth[1];
      ka2 = dnorm((xx2)/roth[2],0,1,log=FALSE)/roth[2];
      kj= ka1*ka2
      temp = temp + t(res[i,])^2%*%kj^2%*%res[j,]^2
      }
    }
  }
  return(temp)
}
