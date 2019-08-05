  #??��????
  library('MASS')
  N=20;T=40;mu=c(0,0);Sigma=matrix(c(1,0,0,1),2,2)
  r0=mvrnorm(N,mu,Sigma)
  r0=t(r0)
  f0=mvrnorm(T,mu,Sigma)
  f0=t(f0)
  n=mvrnorm(N,rep(0,T*2),Sigma=diag(1,nrow=T*2))
  e=mvrnorm(N,rep(0,T),Sigma=diag(1,nrow=T))
  X_1=matrix(0.25,N,T)+0.25*(t(r0)%*%f0)+n[,1:T]
  X_2=matrix(0.5,N,T)+0.5*(t(r0)%*%f0)+n[,(T+1):(2*T)]
#  X_1=n[,1:T]
#  X_2=n[,(T+1):(2*T)]
  #GDP4
  #e=(mvrnorm(N,rep(0,T),Sigma=diag(1,nrow=T)))*(0.1+0.1*(X_1*X_1+X_2*X_2))^(1/2)
  #GDP7
  #e=mvrnorm(N,rep(0,T+1),Sigma=diag(1,nrow=T+1))
  #e=0.5*e[,1:T]+e[,2:(T+1)]
  Y=-11*X_1+25*X_2+t(r0)%*%f0+e
  X=array(c(as.array(X_1),as.array(X_2)))
  dim(X)=c(N,T,2)
  K=2
  start=matrix(-1,K,1)
  p=2
  R=2
  precision_beta = 10^-6
  repMIN=300
  repMAX=10*repMIN
  c0<-c(0.5,1,2)

  #?????��?
  trans=0;
  if(N<T){
    trans=1;
    NN=N; N=T; T=NN;
    Y=t(Y);
    X=aperm(X, c(2,1,3))
  }

  beta=start
  source("code/factor.R")

  #ѡ????ʼֵ
  beta=Inf*matrix(1,length(start),1);   #????beta
  obj0=Inf; #????beta?µ?Ŀ?꺯??
  count=0;
  exitflag=-1
  for (i in 1:repMAX ){
    if(count<repMIN){
      if(i==1){
        st=start
      }else{
        st=start+10*rnorm(length(start),mean=0,sd=1);
      }
    #  print(st)
      h=minimize_obj_method1(Y,X,R,st,precision_beta)
      para=h[['beta']]
      obj=h[['obj']]
      ef=h[['ef']]
      #cat('obj0 ',obj0,'\n')
      #cat('obj ',obj,'\n')

      if (obj < obj0){#objΪĿ?꺯??
        obj0=obj;
        beta=para;
        if(ef>0){
          exitflag=1;
        }else{
          exitflag=-1}
      }
      if(ef>0){
        count=count+1
      }
   #   print(beta)

    }
    #cat('beta',beta,'\n')
  }

  count
  exitflag
  beta

  #????lambda?Ĺ��ƺ?f?��Ƽ??в??��?:
  res1=Y;
  for(k in 1:K){
    res1=res1-beta[k]*drop(X[,,k])
    #return(res1)
  }
  res1
  V=eigen(t(res1)%*%res1)$vectors
  D=eigen(t(res1)%*%res1)$values
  e=cbind(t(V),D)
  newe=e[order(e[,T+1]),]
  newV=t(newe[,1:T])
  f=newV[,1:R]
  for(r in 1:R){
    f[,r]=f[,r]/sqrt(sum(f[,r]^2))
    if (mean(f[,r])<0){
      f[,r]=-f[,r]
    }
  }

  lambda=res1%*%f
  res=res1-lambda%*%t(f) #?в??Ĺ��?

  if(trans==1){
    save=lambda; lambda=f; f=save;
    res=t(res);
    NN=N; N=T; T=NN;
    Y=t(Y);
    X=aperm(X,c(2,1,3))
  }

  lambda
  f
  res

  #roth????
  roth <- NULL
  s1<-sd(X[,,1])
  roth[1]<-c0[1]*s1*(N*T)^(-1/(4+p))
  s2<-sd(X[,,2])
  roth[2]<-c0[1]*s2*(N*T)^(-1/(4+p))
  #roth=c(0.005,0.016)

  #HS????
  #huaN <- floor(log(N*T))+1
  #hmin <- 0.4*(N*T)^(-1/(2.1*p))
  #hmax <- 3*(N*T)^(-1/1000)
  #w <- (hmax/hmin)^(1/(huaN-1))

  #gamaNThs <- c()
  #for(s in 1:(huaN-1) ){
   #roth <- c(w*s*hmin,w^2*s*hmin)

  #??˹?˺???
   temp0=K0(roth,X,res)
   temp1=kii_kernel(roth,X,res)
   temp=Kh_R(roth,X,res)

  #????ͳ??��
   B1NT=1/N/T*(roth[1]*roth[2])^(1/2)*temp1
   Jnt=1/N^2/T^2*temp0
   VNT=1/N^2/T^2*(roth[1]*roth[2])*temp*2
   gamaNT=(N*T*(roth[1]*roth[2])^(1/2)*Jnt-B1NT)/VNT^(1/2)
   gamaNT
   #gamaNThs[s] <- gamaNT
  #}
  #supgamaNT <- max(gamaNThs)
