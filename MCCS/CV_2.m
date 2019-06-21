%****************************************************************%
% CV_2 function computes the value of h_{2} in SM_4 step 
% in the manuscript,
% "Linear Mode Regression with Covariate Measurement Error"
%****************************************************************%
% Arguments:
% parms = starting point of both intercept and slope paramters
% Lambda = Bandwidth 
% x_cont = Contaminated covariate W
% y_o = Response variable
% B = the number of further contaminated covariate data in Section 3.3
% n_o = sample size
% options = control do not show the result from "fslove" function
% Outputs:
% out: MSE
function [ out ] = CV_2( parms, x_cont, y_o,Lambda, B, n_o, options  )
  
 beta_2 = zeros(B,2);
 beta_3 = zeros(B,2);

 for i=1:B
        u_b =  normrnd(0,0.3390,1000,n_o);
        u_b1 =  normrnd(0,0.3390,1000,n_o);
        w_st = x_cont +  normrnd(0,0.3390,1,n_o);
        w_st2 = w_st +  normrnd(0,0.3390,1,n_o);
        f_1 = @(beta)MCCS(beta,Lambda,w_st,y_o,u_b,n_o);
        beta_2(i,:) = fsolve(f_1,parms,options);
        f_2 = @(beta)MCCS(beta,Lambda,w_st2,y_o,u_b1,n_o);
        beta_3(i,:) = fsolve(f_2,parms,options);
 end
 
   beta_f = beta_3-beta_2;
   beta_f( (beta_2(:,1)==2 & beta_2(:,2)==5) | (beta_3(:,1)==2 & beta_3(:,2)==5) , :)=[];
   dim_beta = size(beta_f);
   if dim_beta(1) <= 1
       d_st_2 = NaN(1);
   else
       S_star2 = cov(beta_f);
       d_st_2 = diag((beta_f) / S_star2 * (beta_f)');
   end
   out = mean(d_st_2);

end

