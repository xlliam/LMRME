%****************************************************************%
% MCCS function computes the sum of \Psi_{MC,B} in MC_4 step
% in the manuscript,
% "Linear Mode Regression with Covariate Measurement Error"
%****************************************************************%
% Arguments:
% beta = (\beta_{0},\beta_{1}) intercept and slope paramters
% Lambda = Bandwidth 
% x_cont = Contaminated covariate W
% y_o = Response variable
% u_b = Indepedent random errors generated from Normal distribution
% n_o = sample size
% Outputs:
% out: sum of \Psi_{MC,B} in MC_4 step

function [ out ] = MCCS( beta,Lambda, x_cont, y_o,u_b,n_o)
f = zeros(n_o,2);
for i=1:n_o
% phi_mcb function is used to compute \Psi_{MC,B} in MC_3 Step
  f(i,:) =  phi_mcb(beta,Lambda,x_cont(i),y_o(i),u_b(:,i));
end
out = sum(f,1);
end




