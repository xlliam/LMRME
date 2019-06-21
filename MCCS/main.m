%****************************************************************%
%                                                                %
%  MATLAB code to analyze dietary data by using Monte Carlo      %
%  Corrected Score method (MCCS).                                %
%  corresponding to the manuscript title,                        %
%  "Linear Mode Regression with Covariate Measurement Error"                 %
%                                                                %
%                                    Date: 10/02/2017.             %
%                                                                %
%****************************************************************%
clear;clc;
%****************************************************************%
% Set path to read the dietary data
% You may need specify your path to read the data
%****************************************************************%
addpath('Your Path to Load Data')

%****************************************************************%
% Part I: Data cleaning 
%****************************************************************%
% Data claim:
% y_up: FFQ intake measured as the percent calories from fat, 
%       Response Variable Y.
% w_n1: average of these recalls from each subject as a surrogate 
%       of this subject's long-term intake, covariate W.
% fri_n,i=1,...,6: six 24-hour food recalls on randomly selected days.

% Read data from 'wishreg.csv'
data = readtable('wishreg.csv');
s_ize = size(data);
% Scale and center all variables
fr1_n = (data.fr1 - mean(data.fr1)) / std(data.fr1);
fr2_n = (data.fr2 - mean(data.fr2)) / std(data.fr2);
fr3_n = (data.fr3 - mean(data.fr3)) / std(data.fr3);
fr4_n = (data.fr4 - mean(data.fr4)) / std(data.fr4);
fr5_n = (data.fr5 - mean(data.fr5)) / std(data.fr5);
fr6_n = (data.fr6 - mean(data.fr6)) / std(data.fr6);
y=data.ffq;
y_up = (y-mean(y)) / std(y);
% Obtain contaminated covariate W
w_n1 = ( fr1_n + fr2_n + fr3_n + fr4_n + fr5_n + fr6_n ) / 6 ;

%****************************************************************%
% Part II: Estimate the variance of measurement errros 
%****************************************************************%
fr_t = cat(2,fr1_n,fr2_n,fr3_n,fr4_n,fr5_n,fr6_n);
gama_u2 = sum( sum((fr_t - repmat(w_n1,[1,6]) ).^2,2)) / 5 / s_ize(1);
sigma_u2 = gama_u2 / 6;
%Once run the three lines above, one can have 
%***************% 
% std = 0.3390; %
%***************%

%****************************************************************%
% Part III: Run MCCS method
%****************************************************************%

x_total = w_n1';
y_total = y_up';
% B and options are parameters in CV_1 and CV_2 function.
B=10;
options = optimset('Display','off');

x_con = x_total;
y = y_total; 
% Bandwidth Selection by using SIMEX method, (-0.2672,0.3627) is Naive estimator as starting point
parms = [-0.2672;0.3627];
f1 = @(Lambda)CV_1(parms, x_con, y,Lambda, B, s_ize(1) ,options );
f2 = @(Lambda)CV_2(parms, x_con, y,Lambda, B, s_ize(1) ,options );
h_1 = fminbnd(f1,0.6,1.6);
h_2 = fminbnd(f2,0.6,1.7);
h = h_1^2/h_2;

% Obtain final estimator by using MCCS method
if isnan(h)
    beta = NaN(1,2);
else
    u_b = normrnd(0,0.3390,1000,s_ize(1));
    f = @(beta)MCCS(beta,h,x_con,y,u_b,s_ize(1));
    beta = fsolve(f,[-0.2672;0.3627],options);
end

%compute estimator beta_0 and beta_1
table(beta(1),beta(2),'VariableNames',{'beta_0','beta_1'})
