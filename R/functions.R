####################################################
#Function to calculate M, the dominant temporal dependence
#If n is small, the entire data set is used to compute M
#If n is large (greater than 150), the function takes 10 random samples with a minimum length of 50. The following procedure is repeated for each of the random samples and the mean M value is rounded and returned.
#h, the lag parameter, is incremented by one
#When the difference in the autocorrelation values for consecutive h values falls below 0.05, h and m are estiamted and the loop is broken
#Otherwise, the loop continues checking larger numbers up to a value of 10 (or the largest M possible for small data sets)
#A stopping point of 10 is used to ensure that the function does not become stuck in an infinite loop
####################################################

M.est <- function(data_M,M_threshold){
  n <- dim(data_M)[1]
  m <- NULL
  if (n>150){
    for (i in 1:10){
      left <- sample(1:(n-51),1)
      right <- sample((left+50):n,1)
      data <- data_M[left:right,]
      for(H in 1:10){
        h<- H-1
        cor.mat<-cor.est(H,(data))
        if(abs(cor.mat[1,(2*H+1)]/ cor.mat[(H+1),(H+1)] - (cor.mat[2,(2*H)]/ cor.mat[(H+1),(H+1)])) < M_threshold ){
          m <- c(m,max(0,(h-1)))
          break}
      }
    }
  }
  else
  {
    data <- data_M
    for(H in 1:min(10,floor(((n-4)/3)))){
      h<- H-1
      cor.mat<-cor.est(H,(data))
      if(abs(cor.mat[1,(2*H+1)]/ cor.mat[(H+1),(H+1)] - (cor.mat[2,(2*H)]/ cor.mat[(H+1),(H+1)])) < M_threshold ){
        m <- c(m,max(0,(h-1)))
        break}
    }
  }
  ifelse (is.null(m), m <- 0, m <- round(mean(m)))
  print(paste("dominant.temporal.dependence_M =",m))
  return(m)
}


####################################################
# function to estimate the variance of the test statistic
# under the null hypothesis
# M:   dominant dependence
# n:   the total number of time points
# cor: a (2M+1)x(2M+1) matrix with elements including estimate of
#      tr{C(h1)C(h2)} for 0<=|h1|<=M and 0<=|h2|<=M
####################################################

var.est<-function(M,n,cor){
  B.mat<-matrix(0,n,n)
  for (t in 1:(n-1)){
    E.mat<-matrix(0,n,n)
    E.mat[1:t,1:t]<-(n-t)/t
    E.mat[1:t,(t+1):n]<--1
    E.mat[(t+1):n,1:t]<--1
    E.mat[(t+1):n,(t+1):n]<-t/(n-t)
    B.mat<-B.mat+E.mat
  }
  M.n<-min(M,n-1)
  var<-0
  for(i in 1:(2*M.n+1)){
    for(j in 1:(2*M.n+1)){
      h1<-i-(M.n+1);h2<-j-(M.n+1)
      if(h1<=0){
        id1<-(abs(h1)+1):n;id3<-1:(n-abs(h1))
      }else{
        id1<-1:(n-abs(h1));id3<-(abs(h1)+1):n
      }
      if(h2<=0){
        id2<-1:(n-abs(h2)); id4<-(abs(h2)+1):n
      }else{
        id2<-(abs(h2)+1):n;id4<-1:(n-abs(h2))
      }
      var<-var+sum(B.mat[id1,id2]*(B.mat[id3,id4]+t(B.mat)[id3,id4])*cor[i,j])
    }
  }
  var<-var/(n^4)
  return(var)
}

####################################################
# function to estimate the variance of the test statistic
# under the null hypothesis
# M:   dominant dependence
# n:   the total number of time points
# cor: a (2M+1)x(2M+1) matrix with elements including estimate of
#      tr{C(h1)C(h2)} for 0<=|h1|<=M and 0<=|h2|<=M
####################################################

var.est.nt<-function(M,n,t,cor){
  B.mat<-matrix(0,n,n)
  B.mat[1:t,1:t]<-(n-t)/t
  B.mat[1:t,(t+1):n]<--1
  B.mat[(t+1):n,1:t]<--1
  B.mat[(t+1):n,(t+1):n]<-t/(n-t)
  M.n<-min(M,n-1)
  var<-0
  for(i in 1:(2*M.n+1)){
    for(j in 1:(2*M.n+1)){
      h1<-i-(M.n+1);h2<-j-(M.n+1)
      if(h1<=0){
        id1<-(abs(h1)+1):n;id3<-1:(n-abs(h1))
      }else{
        id1<-1:(n-abs(h1));id3<-(abs(h1)+1):n
      }
      if(h2<=0){
        id2<-1:(n-abs(h2)); id4<-(abs(h2)+1):n
      }else{
        id2<-(abs(h2)+1):n;id4<-1:(n-abs(h2))
      }
      var<-var+sum(B.mat[id1,id2]*(B.mat[id3,id4]+t(B.mat)[id3,id4])*cor[i,j])
    }
  }
  var<-var/(n^4)
  return(var)
}

############################################################
# function to estimate tr{C(h1)C(h2)}
# for 0<=|h1|<=M and 0<=|h2|<=M
# M: dominant dependence; data: entire sequence of data
############################################################

cor.est<-function(M,data){
  n<-dim(data)[1]
  a<-data%*%t(data)
  correlation.est<-matrix(0,(2*M+1),(2*M+1))
  out<-0
  storage.mode(out)<-"double"
  storage.mode(a)<-"double"
  n<-as.integer(n)
  M<-as.integer(M)
  result4<-.Fortran(code3,n,M,a,out=out)

  for(i in 1:(2*M+1)){
    for(j in 1:(2*M+1)){
      h1<-i-(M+1);h2<-j-(M+1)
      h1<-as.integer(h1)
      h2<-as.integer(h2)
      h3<-abs(h1); h4<-abs(h2)
      h3<-as.integer(h3)
      h4<-as.integer(h4)
      result1<-.Fortran(code1,n,M,h1,h2,a,out=out)
      result2<-.Fortran(code2,n,M,h3,a,out=out)
      result3<-.Fortran(code2,n,M,h4,a,out=out)
      correlation.est[i,j]<-result1[[6]]-result2[[5]]-result3[[5]]+result4[[4]]
    }
  }
  return(correlation.est)
}

####################################################
# function to compute the test statistic
# bias: the bias-correction in eqn(2.3) of our manuscript
####################################################

test.stat<-function(M,t,data,bias){
  n<-dim(data)[1]; p<-dim(data)[2]
  M.n<-min(M,n-1)
  if(M==0){
    f<-1/t+1/(n-t)
  }else{
    f<-rep(0,(M.n+1))
    f[1]<-1/t+1/(n-t)
    for (i in 2:(M.n+1)){
      a<-0;b<-0;c<-0
      if(i<=t){
        a<-2*(t+1-i)/(t^2)
      }
      if(i<=(n-t)){
        b<-2*(n-t-i+1)/((n-t)^2)
      }
      for(j in 1:(i-1)){
        d<-0
        if(j<=t & j>=i-n+t){
          d<--2/(t*(n-t))
        }
        c<-c+d
      }
      f[i]<-a+b+c
    }
  }

  mat1<-matrix(data[1:t,],t,p)
  mat2<-matrix(data[(t+1):n,],(n-t),p)
  mean1<-colMeans(mat1)
  mean2<-colMeans(mat2)
  stat1<-t(mean1-mean2)%*%(mean1-mean2)
  stat2<-t(f)%*%bias[1:(M.n+1)]
  stat<-t*(n-t)/(n^2)*(stat1-stat2)
  return(stat)
}

####################################################
#Function to estimate f_t
#Called by BS.CPdetection
#Returns the transpose of f_t, a row vector with (M+1) components
#Used in estimating the standard deviation sigma_t,j for the power enhancement term L0
####################################################

ft <-function(M,t,data){
  n<-dim(data)[1]; p<-dim(data)[2]
  M.n<-min(M,n-1)
  if(M==0){
    f<-1/t+1/(n-t)
  }else{
    f<-rep(0,(M.n+1))
    f[1]<-1/t+1/(n-t)
    for (i in 2:(M.n+1)){
      a<-0;b<-0;c<-0
      if(i<=t){
        a<-2*(t+1-i)/(t^2)
      }
      if(i<=(n-t)){
        b<-2*(n-t-i+1)/((n-t)^2)
      }
      for(j in 1:(i-1)){
        d<-0
        if(j<=t & j>=i-n+t){
          d<--2/(t*(n-t))
        }
        c<-c+d
      }
      f[i]<-a+b+c
    }
  }
  return(t(f))
}


####################################################
# function for the bias-correction term in eqn(2.3)
####################################################

Bias<-function(M,data){
  n<-dim(data)[1]; p<-dim(data)[2]
  if(M==0){
    Mat<-1-1/n
  }else{
    Mat1<-diag(1,M+1)-diag(rep(0:M),M+1)*(1/n)
    Mat<-Mat2<-matrix(0,M+1,M+1)
    for(i in 1:(M+1)){
      for(j in 1:(M+1)){
        if(j==1){
          Mat2[i,j]=-(1-(i-1)/n)*(1/n)
        }else{
          Mat2[i,j]=2/n*(-1-1/n+max(i,j)/n+(i-1)*(j-1)/(n^2))
        }
      }

    }
    Mat<-Mat1+Mat2
  }
  mean.o<-colMeans(data)
  data.cen<-data-matrix(mean.o,n,p,byrow=T)
  Mat3<-data.cen%*%t(data.cen)
  if(M==0){
    vec<-1/n*sum(diag(Mat3))
  }else{
    vec<-NULL
    vec[1]<-1/n*sum(diag(Mat3))
    for(i in 2:(M+1)){
      if((i-1)<n){
        vec[i]<-1/n*sum(diag(Mat3[-(1:(i-1)),-((n-i+2):n)]))
      }else{
        vec[i]<-0
      }
    }
  }

  bias.mat<-solve(Mat)%*%vec
  return(bias.mat)

}



####################################################
#Function to estimate F^{-1}_{n,m}
#Called by BS.CPdetection
#Returns an (M+1)x(M+1) matrix that's used in estimating the standard deviation sigma_t,j for the power enhancement term L0
####################################################

Finversenm<-function(M,data){
  n<-dim(data)[1]; p<-dim(data)[2]
  if(M==0){
    Mat<-1-1/n
  }else{
    Mat1<-diag(1,M+1)-diag(rep(0:M),M+1)*(1/n)
    Mat<-Mat2<-matrix(0,M+1,M+1)
    for(i in 1:(M+1)){
      for(j in 1:(M+1)){
        if(j==1){
          Mat2[i,j]=-(1-(i-1)/n)*(1/n)
        }else{
          Mat2[i,j]=2/n*(-1-1/n+max(i,j)/n+(i-1)*(j-1)/(n^2))
        }
      }

    }
    Mat<-Mat1+Mat2
  }
  return(solve(Mat))
}

####################################################
#Function to estimate Vj
#Called by BS.CPdetection
#Returns Vj, which is an (M+1) component vector that's used in estimating the standard deviation sigma_t,j for the power enhancement term L0
####################################################

Vj<-function(M,data,j){
  n<-dim(data)[1]; p<-dim(data)[2]
  mean.o<-colMeans(data)
  data.cen<-data-matrix(mean.o,n,p,byrow=T)
  data.cen <- data.cen[,j]
  Mat3<-data.cen%*%t(data.cen)
  if(M==0){
    vec<-1/n*sum(diag(Mat3))
  }else{
    vec<-NULL
    vec[1]<-1/n*sum(diag(Mat3))
    for(i in 2:(M+1)){
      if((i-1)<n){
        vec[i]<-1/n*sum(diag(Mat3[-(1:(i-1)),-((n-i+2):n)]))
      }else{
        vec[i]<-0
      }
    }
  }
  return(vec)
}

####################################################
#Function to determine if a statistically significant change has occurred using the power enhancement method
#Returns the value of the variable 'Decision' (which indicates whether or not a statistically significant change was detected); the index of the change point; the value of the test statistic with the added power enhancement term, L0; and the p-value the test statistic with the added power enhancement term
#Called by the BS.CPtracker function
####################################################

BS.CPdetection<-function(M,data,bias,cor,ti,tf,a){
  data.M<-data[ti:tf,]
  n<-dim(data.M)[1]; p<-dim(data.M)[2]
  Decision <- 0; p_value<- 0; stat.test<-0; max.ind<-0; L0 <- 0
  stat.margin<- rep(0,(n-1))
  Fnm <- Finversenm(M,data.M)
  for (t in 1:(n-1)){
    stat.margin[t] <- test.stat(M,t,data.M,bias)
    sdnt<-abs((var.est.nt(M,n,t,cor)))^0.5
    if (sdnt==0){
      next
    }
    if (stat.margin[t]/sdnt < sqrt(2*log(n)-log(log(n)+log(log(n^(1/pi)))))){
      next
    }
    for (j in 1:(p)){
      mat1<-matrix(data.M[1:t,j],t,1)
      mat2<-matrix(data.M[(t+1):n,j],(n-t),1)
      mean1<-mean(mat1)
      mean2<-mean(mat2)
      meandiff<-(mean1-mean2)^2
      biasMat <- Fnm %*% Vj(M,data.M,j)
      vartj <- ft(M,t,data.M)%*%biasMat[1:(min(M,n-1)+1)]
      sdtj <- abs(((n/(t*(n-t)))*vartj))^0.05
      if (sdtj==0){
        next
      }
      if ( abs(meandiff)/sdtj < (log(log(n))*sqrt(log(p)))){
        next
      }
      L0<-c(L0,meandiff)
    }
  }
  sd<-abs((var.est(M,n,cor)))^0.5
  ind<-which.max(stat.margin)
  ifelse(sd==0,stat.test<-0,stat.test<-sum(stat.margin)/sd)
  stat.test<-stat.test+(sum(L0)*sqrt(p))
  if(stat.test > qnorm(1-a)){Decision<-1;max.ind<-ind;p_value<-(1-pnorm(stat.test))}
  result<-c(Decision,ti+max.ind-1,stat.test,p_value)
  return(result)
}

####################################################
#Function to find change points without the power enhancement method
#Called by BS.CPtracker
#Computes the value of the test statistic for the interval, i.e. the sum of the marginal L_t's for each t in the interval
#Identifies which index t corresponds to the largest marginal L_t in the interval
#Returns this index as well as the value of the test statistic for the given interval and the corresponding p-value.
####################################################

NoPE.BS.CPdetection<-function(M,data,bias,cor,ti,tf,a){
  data.M<-data[ti:tf,]
  n<-dim(data.M)[1]; p<-dim(data.M)[2]
  Decision <- 0; stat.test<-0; p_value <- 0; max.ind<-0
  stat.margin<-NULL
  for (tau in 1:(n-1)){
    stat.margin[tau]<-test.stat(M,tau,data.M,bias)
  }
  sd<-abs((var.est(M,n,cor)))^0.5
  ind<-which.max(stat.margin)
  ifelse(sd==0,stat.test<-0,stat.test<-sum(stat.margin)/sd)
  if(stat.test > qnorm(1-a)){Decision<-1;max.ind<-ind;p_value<-(1-pnorm(stat.test))}
  result<-c(Decision,ti+max.ind-1,stat.test,p_value)
  return(result)
}
####################################################
#Function to find change points without the power enhancement method
#Called by WBS
#Computes the value of the test statistic for the interval, i.e. the sum of the marginal L_t's for each t in the interval
#Identifies which index t corresponds to the largest marginal L_t in the interval
#Returns this index as well as the value of the test statistic for the given interval
####################################################

WBS.CPdetection<-function(M,data,bias,cor,ti,tf){
  data.M<-data[ti:tf,]
  n<-dim(data.M)[1]; p<-dim(data.M)[2]
  stat.test<-0; max.ind<-0
  stat.margin<-NULL
  for (tau in 1:(n-1)){
    stat.margin[tau]<-test.stat(M,tau,data.M,bias)
  }
  sd<-abs((var.est(M,n,cor)))^0.5
  max.ind<-which.max(stat.margin)
  ifelse(sd==0,stat.test<-0,stat.test<-sum(stat.margin)/sd)
  result<-c(ti+max.ind-1,stat.test)
  return(result)
}

####################################################
#Function to Track and return the location of any change point(s) under the regular binary segmentation algorithm
#Called by the binary.segmentation function
#Starts by testing the entire interval
#Calls the BS.CPdetection function to determine if there is a statistically significant change in the interval
#If a change point is detected, two new intervals are added to SegmentList by bisecting the original interval at the location of the change point
#Updates the correlation matrix
#Removes the original interval from the list of SegmentList and proceeds to test the next interval until SegmentList is empty
####################################################

BS.CPtracker<-function(M,data,bias,cor,size,dim,ti,tf,a,power_enhancement){
  SegmentList<-c(ti,tf)
  FoundList<-NULL
  pvalues <- NULL
  while (length(SegmentList)>0)
  {
    s1<-SegmentList[1]; e1<-SegmentList[2]
    Tlen<-e1-s1+1
    if(Tlen>2+M){
      ifelse(power_enhancement==TRUE,Cpvec<-BS.CPdetection(M,data,bias,cor,s1,e1,a),Cpvec <-NoPE.BS.CPdetection(M,data,bias,cor,s1,e1,a))
      if(Cpvec[1]==1){
        k1<-Cpvec[2]
        SegmentList<-c(SegmentList,s1,k1,(k1+1),e1)
        FoundList<-c(FoundList,k1)
        pvalues <- c(pvalues,Cpvec[4])
        num<-length(FoundList)+1
        div<-c(0,sort(FoundList),size)
        cen.mat<-att.mat<-NULL
        for(i in 1:num){
          att.mat<-matrix(rep(colMeans(as.matrix(data[(div[i]+1):div[i+1],])),(div[i+1]-div[i])),(div[i+1]-div[i]),dim,byrow=T)
          cen.mat<-rbind(cen.mat,att.mat)
        }
        cor<-cor.est(M,(data-cen.mat))
      }

    }
    SegmentList<-SegmentList[-c(1,2)]
  }
  Found_with_pvalues <- cbind(FoundList,pvalues)
  if(length(FoundList)>1){Found_with_pvalues <- Found_with_pvalues[order(Found_with_pvalues[,1]),]}
  return(Found_with_pvalues)
}


####################################################
#Function to identify the best candidate change point in each successive interval by generating random sub-intervals, find the largest normalized change among these intervals, and then check that the test statistic at the proposed change point exceeds the threshold zeta
#Called by the WBS.CPtracker function
#Generates 2500 random intervals by choosing 5000 random numbers between two endpoints
#Of the 2500, any intervals with a negative length or a positive length that's less than the user's specified minimum size are eliminated
#Either CP.detection1 or CP.detection2 is called to test for change points in each of the approximately 1,000 remaining intervals
#Changemax and changepoint are used as trackers to record the size of the normalized change in each interval and the index at which the largest marginal change occurs. They are updated when CP.detection1 or CP.detection2 identify a larger normalized change than has been recorded so far.
####################################################

WBS<-function(M,data,bias,cor,size,dim,ti,tf,minsize,num_intervals){
  tot_intervals=2*num_intervals
  num_endpts=2*tot_intervals
  n <- dim(data)[1]
  intervals <- sample(ti:tf,num_endpts,replace=T)
  zeta <- sqrt(2*log(n))
  changemax <- 0
  changepoint <- 0
  while (length(intervals)>0)
  {
    sm<-intervals[1]; tm<-intervals[2]
    Tlen<-tm-sm+1
    if(Tlen>max(minsize,2+M)){
      Changecheck<-WBS.CPdetection(M,data,bias,cor,sm,tm)
      if(Changecheck[2]>changemax){
        changemax <- Changecheck[2]
        changepoint <- Changecheck[1]
      }
    }
    intervals <- intervals[-c(1,2)]
  }
  ifelse (changemax>zeta, return (changepoint), return (0))
}

####################################################
#Called by wild.binary.segmentation function
#The analogue of BS.CPdetection for wild binary segmentation
#Tracks the location of any change point(s)
#Uses SegmentList to determine which intervals to test for change points
#Starts with the endpoints of the original interval as the only elements of SegmentList
#Calls the WBS function to perform wild binary segmentation of the interval
#If a change point is detected by WBS, the original interval is bisected at the change point, and the endpoints of the two new intervals are added to SegmentList
#Updates the correlation matrix
#The original interval is removed from SegmentList and the next interval is tested, continuing until SegmentList is empty
####################################################

WBS.CPtracker<-function(M,data,bias,cor,size,dim,ti,tf,minsize,num_intervals){
  SegmentList<-c(ti,tf)
  FoundList<-NULL
  while (length(SegmentList)>0){
    s1<-SegmentList[1]; e1<-SegmentList[2]
    Tlen<-e1-s1+1
    if(Tlen>10){
      Cpvec<- WBS(M,data,bias,cor,size,dim,s1,e1,minsize,num_intervals)
      if(Cpvec[1]!=0){
        k1<-Cpvec[1]
        SegmentList<-c(SegmentList,s1,k1,(k1+1),e1)
        FoundList<-c(FoundList,k1)
        num<-length(FoundList)+1
        div<-c(0,sort(FoundList),size)
        cen.mat<-att.mat<-NULL
        for(i in 1:num){
          att.mat<-matrix(rep(colMeans(as.matrix(data[(div[i]+1):div[i+1],])),(div[i+1]-div[i])),(div[i+1]-div[i]),dim,byrow=T)
          cen.mat<-rbind(cen.mat,att.mat)
        }
        cor<-cor.est(M,(data-cen.mat))
      }

    }
    SegmentList<-SegmentList[-c(1,2)]
  }
  FoundList <- sort(FoundList)
  return(FoundList)
}

####################################################
#Function to perform regular binary segmentation on the data matrix data_M
#Takes as input the data and the critical value alpha
#Calls the BS.CPtracker function to determine the existence and location of any change point(s)
#Returns "No Change Points Found" if there is no statistically significant change in the mean anywhere in the data
#Otherwise returns a list of the change points and the corresponding p-values
####################################################

binary.segmentation <-function(data_M,alpha=.05,power_enhancement=TRUE,M_threshold=0.05){
  data_M <- as.matrix(data_M)
  M <- M.est(data_M,M_threshold)
  bias.est<-Bias(M,data_M)
  cor.mat<-cor.est(M,(data_M))
  n <- nrow(data_M)
  p <- ncol(data_M)
  cpts<-BS.CPtracker(M,data_M,bias.est,cor.mat,n,p,1,n,alpha,power_enhancement)
  ifelse(is.null(cpts),return("No Change Points Found"),return(cpts))
}

####################################################
#Function to perform wild binary segmentation on the data matrix data_M
#Takes as input the data matrix data_M, the minimum length of the random intervals, and the threshold that the test statistic must exceed
#Returns "No Change Points Found" if there is no statistically significant change in the mean for any randomly generated sub-interval and the location of the change points otherwise
####################################################

wild.binary.segmentation <-function(data_M,minsize=15,num_intervals=1250,M_threshold=0.05){
  data_M <- as.matrix(data_M)
  M <- M.est(data_M,M_threshold)
  bias.est<-Bias(M,data_M)
  cor.mat<-cor.est(M,(data_M))
  n <- nrow(data_M)
  p <- ncol(data_M)
  cpts<-WBS.CPtracker(M,data_M,bias.est,cor.mat,n,p,1,n,minsize,num_intervals)
  if(is.null(cpts))
    {return("No Change Points Found")}
  else{
    print("FoundList")
    return(cpts)
  }
}

####################################################
# function to generate data
####################################################

data_Ger<-function(n,p,M,k,Gam,mu){
  data_Mat<-matrix(0,n,p)
  L<-M+1
  Z<-matrix(rnorm(p*(n+L-1)),p*(n+L-1),1) #Gaussian innovation
  #    Z<-matrix(rt(p*(n+L-1),8),p*(n+L-1),1)*sqrt(3/4) # t innovation
  vec.coef<-1/rep(c(L:1),each=p)
  Gam.mat<-t(apply(Gam,1,rep,L))*matrix(vec.coef,ncol=L*p,nrow=p,byrow=T)
  iter<-length(k)+1
  k<-c(0,k,n)
  for(i in 1:iter){
    for (t in (k[i]+1):k[i+1])
    {
      data_Mat[t,]<-matrix((matrix(mu[,i],ncol=1,nrow=p,byrow=F)-Gam.mat%*%Z[((t-1)*p+1):((t+L-1)*p),]),1,p,byrow=F)
    }
  }
  return(data_Mat)
}
