##################################################################################
#                                                                                #
#   R Code to analyze dietary data by usging CK_Normal method                    #
#   corresponding to the manuscript titled,                                      #
#   "Linear Mode Regression with Covariate Measurement Error"		               #
#                                                                                #
#                                                 Date: 10/02/2017.	          #
##################################################################################

#################################################################################
# Set the path to read the dietary data
# You may need specify your path to read the data and 'FFT.cpp' file
##################################################################################
setwd('Your path to Load Data')

##################################################################################
# required packages
##################################################################################
library(Rcpp)
library(RcppArmadillo)
library(nloptr)


##################################################################################
# Rcpp was used to compute the value of objective function in mode
# regression for a given beta
##################################################################################
sourceCpp('FFT.cpp')

##################################################################################
# Arguments in FFT_AP function:
# input  - m input values of K^{*}(t), m = 2^16 in our simulation.
# mconst - positive and negative sign corresponding to m input values in FFT.
# beta   - the interval in t for the m input values of K^{*}(t).
# m      - the numbe of input values of K^{*}(t).
##################################################################################

m = 2^16 
beta <- sqrt((2*pi)/m) #in order to have the same resolution both in t and x
beta2  = (2*pi)/(m*beta); 
input  = seq(-m*beta/2, m*beta/2-beta, by=beta)  ###this is for t
mconst = (-1)^(0:(m-1))


###################################################################################
# Part I: Bandwidth selection by using SIMEX method.
# both cV1 and cV2 are trying to minimize the MSE based on (W,Y).
# This bandwidth selection method was proposed in Delaigle and Hall (2008).
###################################################################################
# CV_1 function computes the value of h_{1} in SM_2 step in the manuscript.
# Arguments:
# start = staring point of estimtors
# x_con = contaminated covariate W
# y = response varialbe
# B = the number of estimators
# n_o = sample size
# Outputs:
# M1: MSE
cV1 <- function(start, lambda, x_con, y, B, n_o)
{
	h <- lambda
	beta_1 <- beta_2 <- matrix(0,B,2)
	d_st_1 <- array(0,B)
	for (i in 1:B)
	{
		w_st <- x_con + rnorm(n_o,0,0.3390)
	    beta_1[i,] <- bobyqa(start,fn=FFT_AP,x_cont=x_con,y_o=y, input=input,Lambda=h,beta=beta, mconst=mconst)$par
	    beta_2[i,] <- bobyqa(start,fn=FFT_AP,x_cont=w_st,y_o=y, input=input,Lambda=h,beta=beta, mconst=mconst)$par
    }
    
    S_star <- cov(beta_2-beta_1)

    d_st_1 <- diag( (beta_2-beta_1)%*%solve(S_star)%*%t(beta_2-beta_1) )

    M1 <- mean(d_st_1)
    
    return(M1)
	
}

# CV_2 function computes the value of h_{2} in SM_4 step
# Arguments:
# start = staring point of estimtors
# x_con = contaminated covariate W
# y = response varialbe
# B = the number of estimators
# n_o = sample size
# Outputs:
# M2: MSE
cV2 <- function(start, lambda, x_con, y, B, n_o)
{
	h <- lambda
    beta_2 <- beta_3 <- matrix(0,B,2)
	d_st_2 <- array(0,B)
	for (i in 1:B)
	{
		w_st <- x_con + rnorm(n_o,0,0.3390)
	    w_st2 <- w_st + rnorm(n_o,0,0.3390)
	    beta_2[i,] <- bobyqa(start,fn=FFT_AP,x_cont=w_st,y_o=y, input=input,Lambda=h,beta=beta, mconst=mconst)$par
        beta_3[i,] <- bobyqa(start,fn=FFT_AP,x_cont=w_st2,y_o=y, input=input,Lambda=h,beta=beta, mconst=mconst)$par
    }
    
        S_star2 <- cov(beta_3-beta_2)

	d_st_2 <- diag( (beta_3-beta_2)%*%solve(S_star2)%*%t(beta_3-beta_2) )

    M2 <-  mean(d_st_2) 
    
    
    return(M2)
	
}

##############################################################################
# Part II: Data cleaning
##############################################################################
# Data claim:
# y_total: FFQ intake measured as the percent calories from fat, Response Variable Y.
# x_total: average of these recalls from each subject as a surrogate of this subject's long-term intake, covariate W.
# (fr1_n, fr2_n, fr3_n, fr4_n, fr5_n, fr6_n): six 24-hour food recalls on randomly selected days.

# read data from 'wishreg.csv'
data = read.csv('wishreg.csv',header=TRUE,sep=",")
n = dim(data)[1]
# scale and center all variables
y_total = (data$ffq-mean(data$ffq)) / sqrt(var(data$ffq))
fr1_n = (data$fr1 - mean(data$fr1)) / sqrt(var(data$fr1));
fr2_n = (data$fr2 - mean(data$fr2)) / sqrt(var(data$fr2));
fr3_n = (data$fr3 - mean(data$fr3)) / sqrt(var(data$fr3));
fr4_n = (data$fr4 - mean(data$fr4)) / sqrt(var(data$fr4));
fr5_n = (data$fr5 - mean(data$fr5)) / sqrt(var(data$fr5));
fr6_n = (data$fr6 - mean(data$fr6)) / sqrt(var(data$fr6));

# obtain contaminated covariate W
x_total = ( fr1_n + fr2_n + fr3_n + fr4_n + fr5_n + fr6_n ) / 6 ;

##############################################################################
# Part III: Run CK_Normal method
##############################################################################

# Bandwidth Selection by using SIMEX method, (-0.2672,0.3627) is Naive estimator as starting point
y <- y_total
x_con <- x_total
h_1 <- optimize(cV1,c(0.5,1.5),start=c(-0.2672,0.3627),x_con = x_con, y = y, B=10, n)$minimum
h_2 <- optimize(cV2,c(0.5,1.7),start=c(-0.2672,0.3627),x_con = x_con, y = y, B=10, n)$minimum
h_select_n <- h_1^2/h_2
# obtain final estimator by using CK_Normal method
beta_e <- bobyqa(c(-0.2672,0.3627),fn=FFT_AP,x_cont=x_con,y_o=y, input=input, Lambda=h_select_n,beta=beta, mconst=mconst)$par

#Compute estimator beta_0 and beta_1
beta_e <- matrix(beta_e,1,2)
colnames(beta_e) = c("beta_0","beta_1")
rownames(beta_e) = c("")
round(beta_e,4)
