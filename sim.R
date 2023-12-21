##############################################################
######      some useful functions      #######################
##############################################################

standvec = function(x){
  return((x - mean(x))/sd(x))}


matpower = function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp[[2]]%*%diag((tmp[[1]])^alpha)%*%
           t(tmp[[2]]))}

CalGam=function(A){
  if(is.vector(A)){
    n = length(A)}
  else n = dim(A)[1]
  tmp=rowSums(as.matrix(A*A))%*%t(rep(1,n))
  K=as.numeric(tmp+t(tmp)-2*A%*%t(A))
  K=K*(K>=0)
  tou=sum(sqrt(K))/(n*(n-1))
  gam=1/(2*tou^2)
  return(gam)
}

KGaussian=function(gamma,A,B){
  if(is.vector(A)){
    n = length(A)}
  else n = dim(A)[1]
  if(is.vector(B)){
    m = length(B)}
  else m = dim(B)[1]
  tmp_1=rowSums(as.matrix(A*A))%*%matrix(1,1,m)
  tmp_2=rowSums(as.matrix(B*B))%*%matrix(1,1,n)
  K=tmp_1+t(tmp_2)-2*A%*%t(B)
  K=exp(-K*gamma)
  return(K)
}


standard = function(x){
  xm = apply(x,2,mean)
  xc = t(t(x) - xm)
  xs = xc%*%matpower(var(x),-1/2)
  return(xs)
}



center = function(x){
  xm = apply(x,2,mean)
  xc = t(t(x) - xm)
  return(xc)
}



margstand  = function(x){
  xm = apply(x,2,mean)
  xc = t(t(x) - xm)
  xs = xc%*%matpower(diag(diag(var(x))),-1/2)
  return(xs)
}



hyper = function(x1,psi,rho){
  x2 = (-rho - psi[1]*x1)/psi[2]
  return(c(x1,x2))}


connect = function(a,b){
  lines(c(a[1],b[1]),c(a[2],b[2]))}



twoclass = function(x,y){
  plot(x[,1],x[,2],pch=" ")
  points(x[y==-1,1],x[y==-1,2])
  points(x[y==1,1],x[y==1,2],pch="+")
}


dis = function(v1,v2){
  if (dim(as.matrix(v1))[2] == 1)
  {p1 = v1%*%t(v1)/c(t(v1)%*%v1)
  p2 = v2%*%t(v2)/c(t(v2)%*%v2)
  d = sqrt(sum((p1-p2)*(p1-p2)))
  return(d)} else
    p1 <- v1%*%matpower(t(v1)%*%v1,-1)%*%t(v1)
  p2 <- v2%*%matpower(t(v2)%*%v2,-1)%*%t(v2)
  d <- sqrt(sum((p1-p2)*(p1-p2)))
  return(d)
}


bench = function(p,d,nmonte){
  tmp = 0
  for(i in 1:nmonte){
    v1 = matrix(rnorm(p*d),p,d)
    v2 = matrix(rnorm(p*d),p,d)
    tmp = tmp + dis(v1,v2)/nmonte}
  return(tmp)
}



mpstand = function(x,ignore){
  return(center(x)%*%mppower(var(x),-1/2,ignore))}



rigpower = function(matrix,power,epsilon){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  p = nrow(matrix)
  matrix1 = matrix + epsilon*eval[1]*diag(p)
  eig1 = eigen(matrix1)
  eval1 = eig1$values
  evec1 = eig1$vectors
  tmp = evec1%*%diag(eval1^power)%*%t(evec1)
  return(tmp)
}



rigstand = function(x,epsilon){
  return(center(x)%*%rigpower(var(x),-1/2,epsilon))}



mvar = function(x){
  return(diag(apply(x,2,var)))}



unevendiv = function(n,proportion){
  nmid = round(n*proportion)
  start = round(n*(1-proportion)/2)
  tmp = seq(from = start, to = start + nmid -1)
  return(tmp)
}



evendiv = function(n,ndiv){
  indiv=round(seq(from=1,to=n,length=ndiv+2)[c(-1,-(ndiv+2))])
  return(indiv)
}



linker = function(u,v){
  return(t(u)%*%v)}



lkermat = function(x){
  n = nrow(x)
  kerm = matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      kerm[i,j] = linker(x[i,],x[j,])}}
  return(kerm)
}
#


lonepsi = function(x,y,div,co){
  y1 = slice(y,div)
  out = svm(x, y1, kernel="linear", type="C", cost = co)
  psi=c(t(out$coefs)%*%margstand(x)[out$index,])
  return(psi)
}



polyker = function(u,v,deg,gamma,coef0){
  return((gamma*t(u)%*%v + coef0)^deg)}



linspear = function(truebeta,beta,x,whichmodel){
  u = center(x)%*%truebeta
  v = center(x)%*%beta
  vv = cbind(1,v[,1],v[,2],v[,1]*v[,2],v[,1]^2,v[,2]^2)
  truepred = modfun(x,whichmodel)
  proj = vv%*%solve(t(vv)%*%vv)%*%t(vv)%*%truepred
  return(cor(rank(proj),rank(truepred)))
}




rescale = function(x){
  p = ncol(x)
  xr = numeric()
  for(i in 1:p){
    num1 = max(x[,i])
    num2 = min(x[,i])
    xr = cbind(xr,(x[,i]-num2)/(num1 - num2))
  }
  return(xr)
}



bic = function(eval,n,constant){
  bicout = numeric()
  for(i in 1:p){
    bicout = c(bicout,sum(eval[1:i])-constant*log(n)*i/sqrt(n))}
  return(bicout)
}



geneigen = function(a,b,ignore,d){
  b1 = mppower(b,-1/2,ignore)
  a1 = b1 %*% a %*% t(b1)
  a1 = (a1 + t(a1))/2
  evec = eigen(a1)$vectors[,1:d]
  return(b1%*%evec)
}

sym = function(a){
  return((a + t(a))/2)}



madm = function(x){
  return(median(abs(x-median(x))))}


qmat = function(n){
  ii = diag(n)
  jj = rep(1,n)%*%t(rep(1,n))
  return(ii-jj/n)}




########################################################################
#                GSIR (This part is borrowed from Prof Kuang-Yao Lee)  #
########################################################################


gsir = function(kx,ky,epsx,epsy,numofdirs){
  n = dim(kx)[1]
  iden = diag(n)
  one = rep(1,n)
  qmat = diag(n) - one %*% t(one) / n
  gx = qmat %*% kx %*% qmat
  gy = qmat %*% ky %*% qmat
  gxinv = matpower(gx + epsx *  iden, -1)
  gyinv = matpower(gy + epsy *  iden, -1)
  eprob = gxinv %*% gx %*% gy %*% gyinv %*% gx %*% gxinv
  saveeig = eigen(eprob)
  return(
    list(eigvec = gxinv %*% saveeig[[2]][,1:numofdirs],
         eigval = saveeig[[1]][1:numofdirs]))
}


###############################################################
#           Moore-Penrose type power                          #
###############################################################


mppower = function(matrix,power,ignore){
  eig = eigen(matrix)
  eval = eig[[1]]
  evec = eig[[2]]
  m = length(eval[abs(eval)>ignore])
  tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
    t(evec[,1:m])
  return(tmp)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



###############################################################
#          Conjoined Conditional Covariance  Operator         #
###############################################################


ccco<-function(a,b,x,sgamx,sgamy,sepsx,sepsy,egamx,egamy,egamz,eepsx, numofdirs){
  
  res=x[,c(a,b)]
  pred=x[,-c(a,b)]
  

  gamy = CalGam(res)*sgamy
  gamx = CalGam(pred)*sgamx
  
  kres = KGaussian(gamy,res,res)
  kpred = KGaussian(gamx,pred,pred)
  
  ignore = 10^(-7) 
  numofdirs = 2
  gsirsave = gsir(kpred,kres,sepsx,sepsy,ignore,numofdirs)
  
  
  
  gsirpred = kpred %*% qmat %*% gsirsave[[1]][,c(1:numofdirs)]
  
 
  xi <- cbind(res[,1],gsirpred)
  xj <- cbind(res[,2],gsirpred)
  
  
  gamx1 = CalGam(xi)*egamx
  gamx2 = CalGam(xj)*egamy
  gamgsirpred = CalGam(gsirpred)*egamz
  
  kx1 = KGaussian(gamx1,xi,xi)
  kx2 = KGaussian(gamx2,xj,xj)
  kgsirpred = KGaussian(gamgsirpred,gsirpred,gsirpred)
  
  
  gramx1=qmat%*%kx1%*%qmat
  gramx2=qmat%*%kx2%*%qmat
  gramgsirpred=qmat%*%kgsirpred%*%qmat
  gpredinv = matpower(gramgsirpred + eepsx*qmat, -1)
  
  robjres <- norm(matpower(gramx1,1/2)%*%matower(gramx2,1/2)-matpower(gramx1,1/2)%*%gramgsirpred%*%gpredinv%*%matpower(gramx2,1/2)),type="F")
 
  return(robjres)
}




############################################################### 
# function: operator norm                                 
############################################################### 
onorm=function(a) return(eigen(round((a+t(a))/2,8))$values[1])


############################################################### 
# function: implementation                                 
############################################################### 
obj <- function(x, method){
  if(method==0){
    
    resmat<-matrix(0,dim(x)[[2]],dim(x)[[2]])
    for( i in 1:(dim(x)[[2]]-1)){
      for(j in (i+1):dim(x)[[2]]){
        resmat[i,j]=ccco(i,j,x,sgamx,sgamy,sepsx,sepsy,egamx,egamy,egamz, eepsx)
      }
    }
    
    
    final<-resmat
    return(final)
  }
  
  if(method==1){
    resmat<-matrix(0,dim(x)[[2]],dim(x)[[2]])
    
    for( i in 1:(dim(x)[[2]]-1)){
      for(j in (i+1):dim(x)[[2]]){
        resmat[i,j]=pccco(i,j,x,sgamx,sgamy,sepsx,sepsy,egamx,egamy,egamz,eepsx)
      }
    }
    
    
    final<-resmat
    return(final)
  }
}

