setwd("D:/experiment/Conference Paper/ECML/ECML2022")
rm(list = ls())

##############################################################################
#
#            sigma : gaussian kernel parameter
#            gamma : maximum coefficient
#            lambda: regularized  parameter
#            eta   : gradient descent step size
#            B     : budget 
#            loss  : hinge loss 
#            SV    : set of support vector 
#
##############################################################################
d_index <- 2

dpath          <- file.path("D:/experiment/online learning dataset/binary C/") 

Dataset        <- c("a9a", "phishing","magic04-N","SUSY20000-N") 

savepath1      <- paste0("D:/experiment/Conference Paper/ECML/ECML2022/ECML2022 Result/",
                         paste0("OKS++-logistic-",Dataset[d_index],".txt"))
savepath2      <- paste0("D:/experiment/Conference Paper/ECML/ECML2022/ECML2022 Result/",
                         paste0("OKS++-logistic-all-",Dataset[d_index],".txt"))

traindatapath    <- file.path(dpath, paste0(Dataset[d_index], ".train"))
traindatamatrix  <- as.matrix(read.table(traindatapath))
trdata           <- traindatamatrix[ ,-1]
ylabel           <- traindatamatrix[ ,1]                                             

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)


x         <- seq(-2,3,1)
sigma     <- 2^(x)
len_sigma <- length(sigma)

U   <- 15
C_0 <- 1
G   <- 1

reptimes <- 10
runtime  <- c(rep(0, reptimes))
errorrate<- c(rep(0, reptimes))
all_p    <- matrix(0,nrow = reptimes, ncol = len_sigma)
all_N    <- matrix(0,nrow = reptimes, ncol = len_sigma)
all_error<- matrix(0,nrow = reptimes, ncol = len_sigma)

all_infor<- matrix(0,nrow = reptimes, ncol = 3*len_sigma)

for(re in 1:reptimes)
{
  order    <- sample(1:length_tr,length_tr,replace = F)   
  
  N        <- c(rep(0, len_sigma))                                 # store the selected times of each kernel parameter
  error    <- c(rep(0, len_sigma))
  L        <- c(rep(0, len_sigma))
  k        <- c(rep(0, len_sigma))
  lambda   <- c(rep(0, len_sigma))
  nabla    <- c(rep(0, len_sigma))
  norm     <- c(rep(0, len_sigma))
  p        <- c(rep(1/len_sigma,len_sigma))
  q        <- p
  
  sv_coe_list <- list(
    svpara1   = array(0,1),
    svpara2   = array(0,1),
    svpara3   = array(0,1),
    svpara4   = array(0,1),
    svpara5   = array(0,1),
    svpara6   = array(0,1)
  )
  sv_max_list <- list(
    svmat1  = matrix(0,nrow = feature_tr,ncol=1),
    svmat2  = matrix(0,nrow = feature_tr,ncol=1),
    svmat3  = matrix(0,nrow = feature_tr,ncol=1),
    svmat4  = matrix(0,nrow = feature_tr,ncol=1),
    svmat5  = matrix(0,nrow = feature_tr,ncol=1),
    svmat6  = matrix(0,nrow = feature_tr,ncol=1)
  )
  
  t1       <- proc.time()  #proc.time()
  tilde_CK <- 0
  tilde_CK_ <- 0
  tem1     <- (G*C_0)^{1/3}*(U*len_sigma)^{2/3}
  for (i in 1:length_tr)
  {
    err   <- 0
    It    <- sample(1:len_sigma, 1, replace = T, prob=p)
    
    tem_svmat    <- sv_max_list[[It]]
    tem_svpara   <- sv_coe_list[[It]]
    
    diff_S_i <- tem_svmat - trdata[order[i], ]
    tem      <- colSums(diff_S_i*diff_S_i)
    fx       <- tem_svpara %*% exp(tem/(-2*(sigma[It])^2))
    fx <- fx[1,1]
    haty = 1
    if(fx < 0)
      haty = -1
    if(haty != ylabel[order[i]])
      err = 1
    tem1 <- exp(-ylabel[order[i]]*fx)
    loss <- log(1+tem1)
    
    k[It] <- k[It]+1
    if(k[It] == 1)
    {
      tem_svmat[,1] <- trdata[order[i], ]
    }else{
      tem_svmat <- cbind(tem_svmat,trdata[order[i], ])
    }
    
    tilde_l     <- loss/p[It]
    tilde_CK    <- tilde_CK+tilde_l
    tilde_CK_   <- tilde_CK_+(tilde_l)^2*q[It]
    nabla[It]   <- nabla[It] + tilde_l
    tem2        <- len_sigma^{1/6}*sqrt(1+nabla[It])
    lambda[It]  <- sqrt(3/4)*U^{4/3}*max(tem1^3,8*tilde_CK)^{-1/6}/tem2
    
    tilde_nabla         <- -1*tem1/(1+tem1)*ylabel[order[i]]/p[It]
    tem_svpara[k[It]]   <- -1*lambda[It]*tilde_nabla
    norm[It]    <- sqrt(norm[It]^2+(lambda[It]*tilde_nabla)^2-2*tilde_nabla*lambda[It]*fx)
    if(norm[It]>U)
    {
      tem_svpara <- tem_svpara*U/norm[It]
      norm[It]   <- U
    }
    
    error[It] <- error[It] + err
    L[It]     <- L[It] + tilde_l

    delta     <- 0.5*tem1/(max(tem1,2*tilde_CK^{1/3}))
    p1        <- (1-delta)
    p2        <- delta/len_sigma
    eta       <- sqrt(2*log(len_sigma))/sqrt(1+tilde_CK_)
    q         <- exp(-1*eta*L)/sum(exp(-1*eta*L))
    p         <- p1*q + p2
    
    sv_max_list[[It]]   <- tem_svmat
    sv_coe_list[[It]]   <- tem_svpara
  }
  t2 <- proc.time()
  runtime[re]    <- (t2 - t1)[3]
  errorrate[re]  <- sum(error)/length_tr
  all_error[re,] <- error
  all_p[re,]     <- p
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_error"),
  alg_name = c("OKS++-logistic"),
  dataname = paste0(Dataset[d_index], ".train"),
  run_time = as.character(runtime),
  err_num  = as.character(errorrate), 
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)

all_infor[,1:len_sigma] = all_N[,1:len_sigma]
all_infor[,(len_sigma+1):(2*len_sigma)] = all_p[,1:len_sigma]
all_infor[,(2*len_sigma+1):(3*len_sigma)] = all_error[,1:len_sigma]

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = T) 
write.table(all_infor,file=savepath2,row.names =TRUE, col.names =FALSE, quote = T) 

sprintf("the candidate kernel parameter are :")
sprintf("%.5f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average error rate is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of error is %.5f in dataset", sd(errorrate))
