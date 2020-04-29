# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

##########################################################
# Input a target DV(vector) and a list of IVS            #
# Find studied lags and Kernelize the values             #
#                                                        #
# PARAMETERS:                                            #
# L=3      # Lag                                         #
# T=3      # Studied total period from now. NOTE: Dt=T/L #
# r=1.2    # Varied window inflation factor              #
# BW=2    # Kernel Bindwidth                            #
# alpha=1  # Pure LASSO; 0 for Ridge                     #
##########################################################

wats <- function(DV=DV, IV_list=IV_list, L=L, Total_period=Total_period,
                 r=r, BW=BW, alpha=alpha){

    P = length(IV_list)
    colnames(DV) <- c("value", "time_stamp")
    ##################################################

    ##################################################

    l1 = Total_period*(r-1)/(r^L-1) # First/shortest window size
    lags <- as.data.frame(matrix(NA,length(DV$time_stamp),L))
    for (z in 1:L){
      for (w in 1: length(DV$time_stamp)){
        lags[w,z] = c(DV$time_stamp[w] - l1*(r^(L+1-z)-1)/(r-1))
      }
    }
    np.index <- max(which(lags$V1<=0),0)
    lags <- subset(lags, V1>0) # Remove if any non-positive lags

    bm <- subset(DV, time_stamp > np.index ,select=value)

    Amatrix <- function(dv,iv){
      Am <- as.data.frame(matrix(NA,nrow(lags),L))
      colnames(iv) <- c("value", "time_stamp")
      tSelect <- cbind(replicate(L,iv$time_stamp))
      ySelect <- cbind(replicate(L,iv$value))
      for (j in 1:nrow(lags)){
        ti <- lags[j,][rep(seq_len(1), each = nrow(iv)), ]
        K = exp(-((ti-tSelect)^2)/BW)
        Am[j,] <- mapply(sum,ySelect*K)/mapply(sum,K)
      }
      return(Am)
    }

    Am.y <- Amatrix(DV,DV)
    Am.yx <- Am.y
    for (i in 1:P){
      Am.yx <- cbind(Am.yx, Amatrix(DV,IV_list[[i]])) #All time lags
    }



  #####################################
  # Use 'glmnet' to solve Elastic Net #
  # Lambda is deterined automatically #
  #####################################
  Am.YX <- as.data.frame(scale(Am.yx, center = TRUE, scale = TRUE));
  # STANDARDIZE X
  Bm <- scale(bm, center = TRUE, scale = FALSE);	# AT LEAST CENTER Y
  dat1 <- as.data.frame(cbind(Am.YX,Bm))

  elas <- glmnet::cv.glmnet(as.matrix(Am.yx), Bm,
                    type.measure="mse",
                    intercept = TRUE,
                    alpha=alpha,family="gaussian")
  beta.hat <- coef(elas, s=elas$lambda.min)
  beta.hat <- coef(elas, s=1)
  return(beta.hat)
}
