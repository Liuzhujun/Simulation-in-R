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
  e=cbind(t(V),D);
  newe=e[order(e[,T+1]),]
  newV=t(newe[,1:T])
  f=newV[,1:R]
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
  N=nrow(Y);#返回Y矩阵行数
  T=ncol(Y);
  K=2;
  SST=sum(diag(Y%*%t(Y)))/N/T;
  beta=st;  # beta最小化初值;
  beta_old=st+Inf;
  obj=Inf;
  diff_obj=-Inf;
  while((max(abs(beta-beta_old))>precision_beta)&&(diff_obj<=SST*10^-10)){
    res=get_residuals(Y=Y,X=X,beta=beta)
    w=principal_components(res=res,R)
    lambda=w[[1]]
    f=w[[2]]
    res=res-lambda%*%t(f);       
    obj_old=obj;              #存旧的目标函数
    obj=sum(diag(t(res)%*%res))/N/T;  #最小二乘目标函数
    diff_obj=obj-obj_old;     
    if(diff_obj<=0){  
      YY=as.vector(Y)
      XX=matrix(0,N*T,K)
      for(k in 1:K){
        xx=drop(X[,,k]);        
        xx=xx-lambda %*% (t(lambda) %*% t(ginv(xx)));    #将lambda投影到x以外的地方
        xx=xx-t(f %*% (t(f) %*% ginv(xx)))     
        XX[,k]=as.vector(xx)  #使 XX 成为NTxK matrix
      }
      beta_old=beta;                 #存旧的beta
      beta=ginv(t(XX)%*%XX)%*%t(XX)%*%YY;      #计算最小二乘估计
    }
  }
  if(diff_obj<=0){
    ef=1;       #good solution found
  }else{
    ef=-1; 
  }
  obj=LS_obj(beta=beta,Y=Y,X=X,R=R); #计算目标函数
  h=list(beta=beta,obj=obj,ef=ef)
  return(h)
}

#核函数
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
