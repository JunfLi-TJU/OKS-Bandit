setwd("D:/experiment/Conference Paper/ECML/ECML2022")
rm(list = ls())

#source("F:/experiment/R/exercise/kernel_xy.R")

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
d_index <- 10

dpath          <- file.path("D:/experiment/online learning dataset/regression/") 

Dataset        <- c("elevators_all","bank_all","slice_all","Year_test","Nor-twitter100000",
                   "ailerons_all","puma32h-all","calhousing","N-twitter100000","N-TomsHardware") 

savepath1      <- paste0("D:/experiment/Conference Paper/ECML/ECML2022/ECML2022 Result/",
                         paste0("IOKS-square-",Dataset[d_index],".txt"))
savepath2      <- paste0("D:/experiment/Conference Paper/ECML/ECML2022/ECML2022 Result/",
                         paste0("IOKS-square-all-",Dataset[d_index],".txt"))

traindatapath    <- file.path(dpath, paste0(Dataset[d_index], ".train"))
traindatamatrix  <- as.matrix(read.table(traindatapath))
trdata           <- traindatamatrix[ ,-1]
ylabel           <- traindatamatrix[ ,1]                                             

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)


x         <- seq(-2,3,1)
sigma     <- 2^(x)
len_sigma <- length(sigma)
kap       <- exp(2/(3*log(length_tr)))
U         <- 1
G_1       <- 2*(U+1)
l_max     <- 1

reptimes  <- 10
runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
all_p     <- matrix(0,nrow = reptimes, ncol = len_sigma)
all_error <- matrix(0,nrow = reptimes, ncol = len_sigma)

all_infor <- matrix(0,nrow = reptimes, ncol = 3*len_sigma)

for(re in 1:reptimes)
{
  order    <- sample(1:length_tr,length_tr,replace = F)   
  
  error    <- c(rep(0, len_sigma))
  k        <- c(rep(0, len_sigma))
  p        <- c(rep(1/len_sigma,len_sigma))
  q        <- p
  eta0     <- 8*l_max*len_sigma^(3/8)/(G_1*U*sqrt(length_tr)*(log(length_tr))^(1/2)) 
  eta      <- c(rep(eta0,len_sigma))
  rho      <- c(rep(2*len_sigma,len_sigma))
  norm     <- c(rep(0, len_sigma))
  lambda   <- c(rep(0, len_sigma))
  norm_tilde_nabla   <- c(rep(0, len_sigma))
  
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
  
  t1    <- proc.time()  #proc.time()
  for (i in 1:length_tr)
  {
    err        <- 0
    tilde_ell  <- c(rep(0, len_sigma))
    It         <- sample(1:len_sigma, 1, replace = T, prob=p)
    
    tem_svmat  <- sv_max_list[[It]]
    tem_svpara <- sv_coe_list[[It]]

    diff_S_i <- tem_svmat - trdata[order[i], ]
    tem      <- colSums(diff_S_i*diff_S_i)
    fx       <- tem_svpara %*% exp(tem/(-2*(sigma[It])^2))
    fx <- fx[1,1] 
    
    err = (fx-ylabel[order[i]])^2
    k[It] <- k[It]+1
    if(k[It] == 1)
    {
      tem_svmat[,1]   <- trdata[order[i], ]
    }else{
      tem_svmat       <- cbind(tem_svmat,trdata[order[i], ])
    }
    
    tilde_nabla       <- 2*(fx-ylabel[order[i]])/p[It]    
    norm_tilde_nabla[It] <- norm_tilde_nabla[It]+tilde_nabla^2
    lambda[It]        <-sqrt(0.5)*U/sqrt(1+norm_tilde_nabla[It])

    tem_svpara[k[It]] <- -1*lambda[It]*tilde_nabla
    norm[It] <- sqrt(norm[It]^2 -2*lambda[It]*tilde_nabla*fx+(lambda[It]*tilde_nabla)^2)
    if(norm[It]>U)
    {
      tem_svpara <- tem_svpara*U/norm[It]
      norm[It]   <- U
    }
    
    error[It]     <- error[It] + err
    if(p[It]>=max(eta))
      tilde_ell[It] <- err/l_max/p[It]
    else
      tilde_ell[It] <- err/l_max/(p[It]+max(eta))
    
    sv_max_list[[It]]   <- tem_svmat
    sv_coe_list[[It]]   <- tem_svpara
    
    ###### Binary search ##################
    upper_lambda <- tilde_ell[It]    #### lambda^\ast < 0
    lambda_ast   <- upper_lambda
    tem          <- q^(-7/8)+eta*(tilde_ell-lambda_ast)
    sum_barpt    <- sum(tem^(-8/7))
    if(abs(sum_barpt - 1)>1e-5)
    {
      lower_lambda <- 0
      lambda_ast   <- lower_lambda
      tem          <- q^(-7/8)+eta*(tilde_ell-lambda_ast)
      sum_barpt    <- sum(tem^(-8/7))
      if(abs(sum_barpt - 1)>1e-5)
      {
        lambda_ast   <- (lower_lambda+upper_lambda)/2 
        tem          <- q^(-7/8)+eta*(tilde_ell-lambda_ast)
        sum_barpt    <- sum(tem^(-8/7))
      }
      while(abs(sum_barpt - 1)>1e-5)
      {
        if(sum_barpt<1)
        {
          lower_lambda  <- lambda_ast
        }else{
          upper_lambda  <- lambda_ast
        }
        lambda_ast      <- (lower_lambda+upper_lambda)/2
        tem             <- q^(-7/8)+eta*(tilde_ell-lambda_ast)
        sum_barpt       <- sum(tem^(-8/7))
      }
    }
    q <- tem^(-8/7)
    delta_ <- 1/length_tr^(3/4)
    p <- (1-delta_)*q+delta_/len_sigma  
    for(j in 1:len_sigma)
    {
      if((1/p[j])>rho[j])
      {
        rho[j] <- 2/p[j]
        eta[j] <- kap*eta[j]
      }
    }
  }
  t2 <- proc.time()
  runtime[re]    <- (t2 - t1)[3]
  errorrate[re]  <- sum(error)/length_tr
  all_error[re,] <- error
  all_p[re,]     <- p
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_error"),
  alg_name = c("IOKS-square"),
  dataname = paste0(Dataset[d_index], ".train"),
  run_time = as.character(runtime),
  err_num  = as.character(errorrate), 
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)

all_infor[,(len_sigma+1):(2*len_sigma)] = all_p[,1:len_sigma]
all_infor[,(2*len_sigma+1):(3*len_sigma)] = all_error[,1:len_sigma]

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = T) 
write.table(all_infor,file=savepath2,row.names =TRUE, col.names =FALSE, quote = T) 

sprintf("the candidate kernel parameter are :")
sprintf("%.5f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average MSE is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of error is %.5f in dataset", sd(errorrate))
