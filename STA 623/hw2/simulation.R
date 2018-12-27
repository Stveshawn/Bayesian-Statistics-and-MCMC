pi=c(0.2,0.7,0.1) #proportion of good, Gil, bad
theta=c(0.2,0.5,0.9) #prob of test being positive
N=100 #the total sample size



likelihood=function(pi,theta,N){
  n=round(N*pi) #size for each status
  S=c(1,2,3) #set of status 1=good, 2=Gil, 3=bad
  s=c(rep(S[1],n[1]),rep(S[2],n[2]),rep(S[3],n[3]))
  test=rep(0,N)
  set.seed(1)
  for(i in 1:3){
    test[s==i]=rbinom(n[i],1,theta[i])
  }
  
  #Plug-in
  df=data.frame(test,s)
  theta.hat=aggregate(test~s,df,FUN='mean')$test
  #define the function to calculate likelihood
  L1=function(x,j){
    return(theta.hat[j]^x*(1-theta.hat[j])^(1-x))
  }
  
  #Full Bayes
  alpha=1
  beta=1
  L2=function(x,j){
    nj=sum(s==j)
    xjp=sum(test[s==j])
    if(x==1){
      r=(xjp+alpha)/(nj+alpha+beta)
    } else {
      r=(nj+beta-xjp)/(nj+alpha+beta)
    }
    return(r)
  }
  
  result=matrix(rep(0,12),nrow=6)
  it=0
  rname=NULL
  for(j in 1:3){
    for(x in 0:1){
      it=it+1
      result[it,]=c(L1(x,j),L2(x,j))
      rname=c(rname,paste0('x=',x,',','j=',j))
    }
  }
  colnames(result)=c('Plug-in','Full Bayes')
  rownames(result)=rname
  return(result)
}

##change prior mean and posterior mean to show the changes