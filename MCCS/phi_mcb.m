%****************************************************************%
% phi_mcb function computes the value of \Psi_{MC,B} in MC_3 step
% in the manuscript,
% "Linear Mode Regression with Covariate Measurement Error"
%****************************************************************%
% Arguments:
% beta = (\beta_{0},\beta_{1}) intercept and slope paramters
% Lambda = Bandwidth 
% x_cont = Contaminated covariate W
% y_o = Response variable
% u_b = Indepedent random errors generated from Normal distribution
% Outputs:
% out: the value of \Psi_{MC,B} in MC_3 step

function [out]=phi_mcb(beta,Lambda,x_cont,y_o,u_b)
w_con = complex(x_cont,u_b);
I_one = ones(1000,1);
x_new = [I_one,w_con];
const = 1/sqrt(2*pi)*exp( - ( y_o - x_new*beta).^2/2/Lambda^2 ).* (( y_o - x_new*beta)/Lambda^2 );
const_1 = [const,const];
re = real(x_new.*const_1);
out = mean(re,1);
end

