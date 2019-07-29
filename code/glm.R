myglm<-function(y,x,distribution,theta=0,w=1,max_iter=500,acc=1e-3)
{ # B(N,p) link=logit
  # poisson(lambda) link=log
  require(MASS)
  n=nrow(x)
  m=ncol(x)
  y=matrix(y,n,length(y)/n)
  if (length(w)==1) w=rep(1,n)

  #initial value,default start from (0,0,...,0)
  if(length(theta)==1)  theta.new=rep(0,m)  else theta.new=theta

  if((distribution=="binom") & (ncol(y)==1))# y ~ 0-1 distribution
  {for (i in 1:max_iter)
  {
    theta.old=theta.new
    G=cmp_G1(y,x,theta.old,w=w)$G
    dG=cmp_dG1(y,x,theta.old,w=w)
    theta.new=theta.old-solve(dG)%*%G
    if(hasNaN(theta.new))
    {
      cat("Nan occur!!!\n")
      break
    }
    error=(sum(theta.new-theta.old)^2)
    if((error < acc)) break
  }}else if((distribution=="binom") & (ncol(y)==2)){# y ~ B(N,p)
    for (i in 1:max_iter)
    {
      theta.old=theta.new
      G=cmp_G2(y,x,theta.old,w=w)$G
      dG=cmp_dG2(y,x,theta.old,w=w)
      theta.new=theta.old-solve(dG)%*%G

      if(hasNaN(theta.new))
      {
        cat("Nan occur!!!\n")
        break
      }

      error=(sum(theta.new-theta.old)^2)
      if((error < acc)) break
    }

  }
  else if (distribution=="poisson"){ # y ~ poisson(lambda)
    for (i in 1:max_iter)
    {
      theta.old=theta.new
      G=cmp_G3(y,x,theta.old,w=w)$G
      dG=cmp_dG3(y,x,theta.old,w=w)
      theta.new=theta.old-solve(dG)%*%G
      if(hasNaN(theta.new))
      {
        cat("Nan occur!!!\n")
        break
      }
      error=(sum(theta.new-theta.old)^2)
      if((error < acc)) break
    }

  }
  if(i==max_iter) cat("NOT converge\n")
  return (list(theta=theta.old,steps=i))
}

# y~B(1,mu)
cmp_G1<-function(y,x,theta,del=0,w=1,c=1.5)
{
  n=nrow(x)
  m=ncol(x)
  if (length(w)==1) w=rep(1,n)
  if (length(del)==1) del=rep(0,m)
  theta=theta+del
  mu=1/(1+exp(-x%*%theta))
  #E
  E=matrix(0,m,1)
  for (i in 1:n) {
    r0=(0-mu[i])/(sqrt(mu[i]*(1-mu[i])))
    r1=(1-mu[i])/sqrt(mu[i]*(1-mu[i]))
    E=E+t((min(c,max(-c,r0))*(1-mu[i])+min(c,max(-c,r1))*(mu[i]))*w[i]*(sqrt(mu[i]*(1-mu[i])))*t(x[i,]))
  }
  E=E/n
  res=matrix(0,m,1)
  for (i in 1:n) {
    r_i=(y[i]-mu[i])/sqrt(mu[i]*(1-mu[i]))
    res=res+t(min(c,max(-c,r_i))*w[i]*(sqrt(mu[i]*(1-mu[i])))*t(x[i,]))-E
  }
  return(list(G=res,E=E))
}


cmp_dG1<-function(y,x,theta,c=1.5,w=1)
{
  n=nrow(x)
  m=ncol(x)
  if (length(w)==1) w=rep(1,n)
  dG=matrix(0,m,m)
  mu=1/(1+exp(-x%*%theta))
  for (i in 1:n) {

    dmu=t(t(x[i,])*mu[i]*(1-mu[i]))
    t1=sqrt(mu[i]*(1-mu[i]))
    dt1=t1*(1/2)*(dmu/(mu[i]) - dmu/(1-mu[i]))

    t2=1/(sqrt(mu[i]*(1-mu[i])))
    dt2=t2*(-1/2)*(dmu/mu[i]-dmu/(1-mu[i]))

    t3=sqrt(mu[i]/(1-mu[i]))
    dt3=t3*(1/2*dmu/mu[i]+1/2*dmu/(1-mu[i]))

    t4=mu[i]^(3/2)*(1-mu[i])^(1/2)
    dt4=t4*((3/2)*1/mu[i]-(1/2)*1/(1-mu[i]))*dmu

    t5=mu[i]^(1/2)*(1-mu[i])^(3/2)
    dt5=t5*((1/2)*1/mu[i]-(3/2)*1/(1-mu[i]))*dmu


    r_i=(y[i]-mu[i])/sqrt(mu[i]*(1-mu[i]))
    r0=(0-mu[i])/(sqrt(mu[i]*(1-mu[i])))
    r1=(1-mu[i])/sqrt(mu[i]*(1-mu[i]))

    dG=dG+w[i]*min(c,max(-c,r_i))*dt1%*%t(x[i,])

    if((r_i< -c)|(r_i>c)){
      dG=dG+matrix(0,m,m)
    }
    else{
      dG=dG+w[i]*t1*(y[i]*dt2-dt3)%*%t(x[i,])
    }

    dG=dG+w[i]*(min(c,max(-c,r0))*dt5+min(c,max(-c,r1))*dt4)%*%t(x[i,])

    if((r0< -c)|(r0>c)){
      dG=dG+matrix(0,m,m)
    }
    else{
      dG=dG+w[i]*((-dt3)*t5)%*%t(x[i,])
    }
    if((r1< -c)|(r1>c)){
      dG=dG+matrix(0,m,m)
    }
    else{
      dG=dG+w[i]*((dt2-dt3)*t4)%*%t(x[i,])
    }
  }
  dG
}

#compute for y~B(N,p)
cmp_G2<-function(y,x,theta,del=0,w=1,c=1.5)
{
  n=nrow(x)
  m=ncol(x)
  N=y[1,1]+y[1,2]
  y1=y[,1]
  y0=y[,2]
  if (length(w)==1) w=rep(1,n)
  if (length(del)==1) del=rep(0,m)
  theta=theta+del
  p=1/(1+exp(-x%*%theta))
  mu=N*p
  #E
  E=matrix(0,m,1)
  for (i in 1:n) {
    mu_i=mu[i]
    dmu=N*t(t(x[i,])*p[i]*(1-p[i]))

    V_i=N*p[i]*(1-p[i])
    j1=floor(mu_i-c*sqrt(V_i))
    j2=floor(mu_i+c*sqrt(V_i))
    E=E+(c*(1-pbinom(j2,N,p[i])-pbinom(j1,N,p[i]))+
           mu_i/sqrt(V_i)*(
             pbinom(j2-1,N-1,p[i])-pbinom(j1-1,N-1,p[i])-(pbinom(j2,N,p[i])-pbinom(j1,N,p[i])) )
    )/sqrt(V_i)*dmu
  }
  E=E/n
  res=matrix(0,m,1)
  for (i in 1:n) {
    r_i=(y1[i]-mu_i)/sqrt(V_i)
    res=res+min(c,max(-c,r_i))*w[i]/(sqrt(V_i))*dmu-E
  }
  return(list(G=res,E=E))
}


cmp_dG2<-function(y,x,theta,c=1.5,w=1)
{
  n=nrow(x)
  m=ncol(x)
  N=y[1,1]+y[1,2]
  y1=y[,1]
  y0=y[,2]
  if (length(w)==1) w=rep(1,n)
  dG=matrix(0,m,m)
  p=1/(1+exp(-x%*%theta))
  mu=N*p
  for (i in 1:n) {
    dmu=N*t(t(x[i,])*p[i]*(1-p[i]))
    dp=t(t(x[i,])*p[i]*(1-p[i]))
    V_i=N*p[i]*(1-p[i])
    dV=N*t(t(x[i,])*p[i]*(1-p[i]))*(1-2*p[i])
    r_i=(y1[i]-mu[i])/sqrt(V_i)

    t1=N^(1/2)*p[i]^(1/2)*(1-p[i])^(-1/2)
    dt1=t1/2*(1/p[i]+1/(1-p[i]))*dp

    TEMP=matrix(0,m,m)
    j1=floor(mu[i]-c*sqrt(V_i))
    j2=floor(mu[i]+c*sqrt(V_i))
    dE=(c*(-cmp_dE2(p[i],j2,N)-cmp_dE2(p[i],j1,N))+
             t1*(cmp_dE2(p[i],j2-1,N-1)-cmp_dE2(p[i],j1-1,N-1)-
                                cmp_dE2(p[i],j2,N)+cmp_dE2(p[i],j1,N)))*dp+
      (pbinom(j2-1,N-1,p[i])-pbinom(j1-1,N-1,p[i])-(pbinom(j2,N,p[i])-pbinom(j1,N,p[i])))*dt1

    t0=sqrt(p[i]*(1-p[i]))
    dt0=t0/2*(1/p[i]-1/(1-p[i]))*dp
    if((r_i< -c)|(r_i>c)){
      TEMP=-dE
    }else{
      TEMP=((y1[i]*(-1/2)*(V_i)^(-3/2)*dV-dt1)-dE)
    }

    dG=dG+N^(1/2)*w[i]*t0*(TEMP)%*%(x[i,])

    E=(c*(1-pbinom(j2,N,p[i])-pbinom(j1,N,p[i]))+
         mu[i]/sqrt(V_i)*(
           pbinom(j2-1,N-1,p[i])-pbinom(j1-1,N-1,p[i])-(pbinom(j2,N,p[i])-pbinom(j1,N,p[i])) )
    )/sqrt(V_i)*dmu
    dG=dG+N^(1/2)*((min(c,max(-c,r_i))-E)*w[i]*dt0)%*%x[i,]

  }
  dG
}

cmp_dE2<-function(p,k,N){
  dE=0
  if(k>=0) {
    for (i in 0:k) {
      dE=dE+dbinom(i,N,p)*(i/p-(N-i)/(1-p))
    }
  }
  dE
}


#compute for y~poisson(lambda)
cmp_G3<-function(y,x,theta,del=0,w=1,c=1.5)
{
  n=nrow(x)
  m=ncol(x)
  if (length(w)==1) w=rep(1,n)
  if (length(del)==1) del=rep(0,m)
  theta=theta+del
  lambda=exp(x%*%theta)
  #E
  E=matrix(0,m,1)
  for (i in 1:n) {
    mu_i=V_i=lambda[i]
    dmu=t(t(x[i,])*lambda[i])
    j1=floor(mu_i-c*sqrt(V_i))
    j2=floor(mu_i+c*sqrt(V_i))
    E=E+(c*(1-ppois(j2,lambda[i])-ppois(j1,lambda[i]))+
           sqrt(V_i)*(dpois(j1,lambda[i]) -dpois(j2,lambda[i]))
    )/sqrt(V_i)*dmu
  }
  E=E/n
  res=matrix(0,m,1)
  for (i in 1:n) {
    r_i=(y[i]-mu_i)/sqrt(V_i)
    res=res+min(c,max(-c,r_i))*w[i]/(sqrt(V_i))*dmu-E
  }
  return(list(G=res,E=E))
}


cmp_dG3<-function(y,x,theta,c=1.5,w=1)
{
  n=nrow(x)
  m=ncol(x)
  if (length(w)==1) w=rep(1,n)
  dG=matrix(0,m,m)
  p=(exp(x%*%theta))

  for (i in 1:n) {
    dmu=t(t(x[i,])*p[i])
    r_i=(y[i]-p[i])/sqrt(p[i])
    TEMP=matrix(0,m,m)
    mu_i=p[i]
    V_i=mu[i]
    j1=floor(mu_i-c*sqrt(V_i))
    j2=floor(mu_i+c*sqrt(V_i))
    dE=(c*(-cmp_dE3(p[i],j2)-cmp_dE3(p[i],j1))+
             sqrt(p[i])*(dpois(j1,p[i])*(j1-p[i])-dpois(j2,p[i])*(j2-p[i])))*(x[i,])+
      (p[i]^(-1/2))/2*(dpois(j1,p[i])-dpois(j2,p[i]))*dmu
    if((r_i< -c)|(r_i>c)){
      TEMP=-dE
    }else{
      TEMP=(y[i]*(-1/2)*(mu[i])^(-3/2)*dmu-1/2*(mu[i])^(-1/2)*dmu)-(dE)
    }

    dG=dG+w[i]*sqrt(mu[i])*(TEMP)%*%(x[i,])


    E=(c*(1-ppois(j2,mu[i])-ppois(j1,mu[i]))+
         mu_i/sqrt(V_i)*(dpois(j1,p[i]) -dpois(j2,p[i]))
    )
    dG=dG+(min(c,max(-c,r_i))-E)*w[i]*1/2*sqrt(mu[i])*(x[i,])%*%t(x[i,])

  }
  dG
}

cmp_dE3<-function(lambda,k){
  dE=0
  if(k>=0) {
    for (i in 0:k) {
      dE=dE+dpois(i,lambda)*(i-lambda)
    }
  }
  dE
}

hasNaN<-function(vec){
  for (i in 1:length(vec)) {
    if (is.nan(vec[i])) return (TRUE)
  }
  if(i==length(vec)) return(FALSE)
}


#
# n=500
# m1=2
# m2=0
# m3=2
# m=m1+m2+m3
# beta=c(rep(-1,m1),rep(0,m2),rep(1,m3))
# X=runif(n*m,-1,1)
# X=matrix(X,n,m)
# eta=X%*%beta
# mu=1/(1+exp(-eta))
#
# # Y~B(1,p)
# Y=runif(n)
# Y[Y>=mu]=0
# Y[Y>0]=1
# glm(Y~X+0,family=binomial(link="logit"))
# myglm(Y,X,distribution  = "binom")
#
# # Y~B(n,p)
# N=5
# y1=rbinom(n,N,mu)
# y2=N-y1
# Y=cbind(y1,y2)
# glm(Y~X+0,family=binomial(link="logit"))
# myglm(Y,X,distribution = "binom")
#
# # Y~poisson(lambda)
# lambda=exp(X%*%beta)
# Y=rpois(n,lambda)
# glm(Y~X+0,family = poisson(link = "log"))
# myglm(Y,X,distribution = "poisson")
