<!-- README.md is generated from README.Rmd. Please edit that file -->
HDcpDetect
==========

The package contains two functions - binary.segmentation and wild.binary.segmentation. The first iteratively uses a bisection algorithm that first tests for the existence of change point(s) and then estimates their location(s). The algorithm starts over the entire interval and continues on sub-intervals until no more change points can be identified. The second (wild binary segmentation) generates thousands of random sub-intervals within the data and computes the test statistic for each of these intervals. The largest of these values is then compared to a threshold, and the original interval is divided in half if a change point is found. This algorithm continues until no more change points can be identified. The binary.segmentation function returns an estimate of M (the dominant temporal dependence), a list of the detected change points, and the corresponding p-values for each change. The wild.binary.segmentation function returns an estimate of M and a list of detected change points.

References: Li, J., Li, L., Xu, M., Zhong, P (2018). Change Point Detection in the Mean of High-Dimensional Time Series Data under Dependence. Manuscript. Fryzlewicz, P. (2014). Wild Binary Segmentation for Multiple Change-point Detection. The Annals of Statistics.

Example
-------

The code below simulates a 150x200 time series with temporal dependence and change points located at 70 and 80. Wild binary segmentation is called on this simulated data matrix to assess its accuracy. The results are tracked over several iterations, and the probability of detecting a change point at each point in the time series is plotted. The number of true positives, false positives, and false negatives is recorded for each iteration, and the mean and standard deviation of each is returned at the end:

``` r
library(ggplot2)
library(ggpubr)
#> Loading required package: magrittr
library(tidyverse)
#> -- Attaching packages ---------------------------------------------------------------------------------------------------------- tidyverse 1.2.1 --
#> v tibble  1.4.2     v purrr   0.2.5
#> v tidyr   0.8.1     v dplyr   0.7.6
#> v readr   1.1.1     v stringr 1.3.1
#> v tibble  1.4.2     v forcats 0.3.0
#> -- Conflicts ------------------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
#> x tidyr::extract()   masks magrittr::extract()
#> x dplyr::filter()    masks stats::filter()
#> x dplyr::lag()       masks stats::lag()
#> x purrr::set_names() masks magrittr::set_names()
library(HDcpDetect)

#################################################
#Compute the mean/sd of true positives and false positives
#################################################
means_n_sds <- function(true_cpts,data_mat){
  tp_vec<-NULL
  fp_vec<-NULL
  for (i in 1:nrow(data_mat)){
    tp_count<-0
    for (j in 1:ncol(data_mat)){
      if (data_mat[i,j] %in% true_cpts){
        tp_count<-tp_count+1
      }
    }
    tp_vec<-c(tp_vec,tp_count)
  }
  for (k in 1:nrow(data_mat)){
    fp_count<-0
    for (l in 1:ncol(data_mat)){
      if (data_mat[k,l] %in% true_cpts == FALSE & data_mat[k,l] != 0){
        fp_count <- fp_count + 1
      }
    }
    fp_vec <- c(fp_vec,fp_count)
  }
  tp.mean<-mean(tp_vec)
  tp.sd<-sd(tp_vec)
  fp.mean<-mean(fp_vec)
  fp.sd<-sd(fp_vec)
  fn.mean<-length(true_cpts)-tp.mean
  fn.sd<-tp.sd
  results <- cbind(tp.mean, tp.sd, fp.mean, fp.sd, fn.mean, fn.sd)
  colnames(results) <- c("tp.mean", "tp.sd","fp.mean","fp.sd","fn.mean","fn.sd")
  return(results)
}

#################################################
#Data generation function
#################################################

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

#################################################
#Wild Binary Segmentation simulation
#################################################
p<-200; n<-150
rat1<-0.4667
rat2 <- 0.5333
k<-c(round(n*rat1),round(n*rat2)) # location of two change points
M<-1
alpha<-0.05
###################################################
# generate different population mean vectors
###################################################
s_beta<-0.3;delta1<-0;delta2<-2;delta3 <- 0
mm<-round(p^(1-s_beta)) # number of nonzero signals
mu1<-mu2<-mu3<-rep(0,p)
signal.sign1<-sample(c(-1,1),mm,replace=T)
signal.sign2<-sample(c(-1,1),mm,replace=T)
signal.sign3<-sample(c(-1,1),mm,replace=T)
post1<-sort(sample(1:p,mm,replace=F))
post2<-sort(sample(1:p,mm,replace=F))
post3<-sort(sample(1:p,mm,replace=F))
mu1[post1]<-signal.sign1*delta1
mu2[post2]<-signal.sign2*delta2
mu3[post3]<-signal.sign3*delta3
mu<-cbind(mu1,mu2,mu3)
#####################################################
# generate covariance matrices
#####################################################

Gam<-matrix(0,p,p)
w<-0.5
for(i in 1:p)
  for(j in 1:p)
  {
    dij<-abs(i-j)
    if(dij<(p*w))
    {
      Gam[i,j]<-w^(dij)
    }
  }
times<-1:p; rho<-0.6
H<-abs(outer(times, times, "-"))
Gam<-rho^H

iter<-10 # the number of iterations

Allcpts<-matrix(0,iter,30)
list <- NULL
for (i in 1:iter){
  data_M<-data_Ger(n,p,M,k,Gam,mu)
  cpts<-wild.binary.segmentation(data_M, minsize = 15)
  print(cpts)
  if (!is.character(cpts))
  {Allcpts[i,1:length(cpts)]<-cpts
    list <- c(list,cpts)
  }
  
}
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 48 70 80
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1]  40  67  74  93 117
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 16 70 80
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 70 80
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 70 78
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 70 80
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 70 80
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 69 79
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 27 80
#> [1] "dominant.temporal.dependence_M = 1"
#> [1] "FoundList"
#> [1] 56 70
if (is.null(list)){list<-0}

df <- data.frame(list)
colnames(df) <- "cpts"
probs <- df%>%dplyr::count(cpts)%>%dplyr::mutate(probability = n/ iter)
all <- data.frame(1:150,0); colnames(all) <- c("point","probability")
for (i in 1:150)
{
  ifelse(i %in% probs$cpts, all[all$point == i,]$probability <- probs[probs$cpts == i,]$probability, all[all$point == i,]$probability <- 0)
  
}
wildbinaryseg <- ggplot(all,aes(x=point,y=probability))+geom_point()+geom_line()+labs(title="Wild Binary Segmentation Simulation (cpts at 70,80)", y= "Probability of Detection")+theme_bw()+theme(plot.title = element_text(hjust = 0.5))+xlim(0,150)+ylim(0,1)

numeric.wildbinaryseg <- means_n_sds(c(70,80),Allcpts)

#wildbinaryseg
numeric.wildbinaryseg
#>      tp.mean     tp.sd fp.mean   fp.sd fn.mean     fn.sd
#> [1,]     1.3 0.8232726     1.2 1.47573     0.7 0.8232726
```
