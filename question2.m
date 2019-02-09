%Problem 2
%This code is written by Zhongming Shi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;
rng(2019)

global beta_all kappa nu

% Set the parameters
beta = .99^(1/12); % "The discount factor ... its annual value is 0.99."
N = 1000; % number of households
T = 12*40; % each month is a time period
sig2_u = .2; sig2_eps = .2;
nu = 1;

% Calibrating kappa
theta = .6; % labor share of output
hour_month = 28.5/7*30;
% mean monthly hours worked per adult in poor countries (Bick et al (2008))
out_con = 1/.25; % output over consumption
kappa = theta*hour_month*out_con^(-1-1/nu);

% Deterministic seasonal component, g(m)
gm_low = [-.073 -.185 .071 .066 .045 .029 .018 .018 .018 .001 -.017 -.041];
gm_mid = [-.147 -.37 .141 .131 .09 .058 .036 .036 .036 .002 -.033 -.082];
gm_hig = [-.293 -.739 .282 .262 .18 .116 .072 .072 .072 .004 -.066 -.164];

gm_low = exp(repmat(gm_low,N,T/12));
gm_mid = exp(repmat(gm_mid,N,T/12));
gm_hig = exp(repmat(gm_hig,N,T/12));
% exponnetial of deterministic seasonal component of all households across
% all periods; will drop "_all" to keep names of variables simple

gm_low_neg = 1./gm_low;
gm_mid_neg = 1./gm_mid;
gm_hig_neg = 1./gm_hig;
% exponential of negative deterministic seasonal component of all
% households across all periods; using simple algebra: e^(x)*e^(-x) = 1

% Stochastic seasonal component, sig2_m
sig2_m_low = [.043 .034 .145 .142 .137 .137 .119 .102 .094 .094 .085 .068];
sig2_m_mid = [.085 .068 .29 .283 .273 .273 .239 .205 .188 .188 .171 .137];
sig2_m_hig = [.171 .137 .58 .567 .546 .546 .478 .41 .376 .376 .341 .273];

beta_all = beta.^(12:T+11); % a sequence, beta^12, bata^13 ... beta^(T+11)

% Pos. correlated stochastic seasonal components of consumption and labor
for i = 1:N
    for k = 1:12
        cov_low = [sig2_m_low(k),.03;.03,sig2_m_low(k)];
        cov_mid = [sig2_m_mid(k),.03;.03,sig2_m_mid(k)];
        cov_hig = [sig2_m_hig(k),.03;.03,sig2_m_hig(k)];
        
        tmp_low = mvnrnd([0 0],cov_low,T/12);
        tmp_mid = mvnrnd([0 0],cov_mid,T/12);
        tmp_hig = mvnrnd([0 0],cov_hig,T/12);
        
        for j = 1:(T/12)
            con_ssc_low(i,k+(j-1)*12)=exp(-sig2_m_low(k)/2+tmp_low(j,1));
            con_ssc_mid(i,k+(j-1)*12)=exp(-sig2_m_mid(k)/2+tmp_mid(j,1));
            con_ssc_hig(i,k+(j-1)*12)=exp(-sig2_m_hig(k)/2+tmp_hig(j,1));
            
            lab_ssc_low(i,k+(j-1)*12)=exp(-sig2_m_low(k)/2+tmp_low(j,2));
            lab_ssc_mid(i,k+(j-1)*12)=exp(-sig2_m_mid(k)/2+tmp_mid(j,2));
            lab_ssc_hig(i,k+(j-1)*12)=exp(-sig2_m_hig(k)/2+tmp_hig(j,2));
        end
    end
end

% Neg. correlated stochastic seasonal components of consumption and labor
for i = 1:N
    for k = 1:12
        cov_low = [sig2_m_low(k),-.03;-.03,sig2_m_low(k)];
        cov_mid = [sig2_m_mid(k),-.03;-.03,sig2_m_mid(k)];
        cov_hig = [sig2_m_hig(k),-.03;-.03,sig2_m_hig(k)];
        
        tmp_low = mvnrnd([0 0],cov_low,T/12);
        tmp_mid = mvnrnd([0 0],cov_mid,T/12);
        tmp_hig = mvnrnd([0 0],cov_hig,T/12);
        
        for j = 1:(T/12)
            con_ssc_low_neg(i,k+(j-1)*12)=...
                exp(-sig2_m_low(k)/2+tmp_low(j,1));
            con_ssc_mid_neg(i,k+(j-1)*12)=...
                exp(-sig2_m_mid(k)/2+tmp_mid(j,1));
            con_ssc_hig_neg(i,k+(j-1)*12)=...
                exp(-sig2_m_hig(k)/2+tmp_hig(j,1));
            
            lab_ssc_low_neg(i,k+(j-1)*12)=...
                exp(-sig2_m_low(k)/2+tmp_low(j,2));
            lab_ssc_mid_neg(i,k+(j-1)*12)=...
                exp(-sig2_m_mid(k)/2+tmp_mid(j,2));
            lab_ssc_hig_neg(i,k+(j-1)*12)=...
                exp(-sig2_m_hig(k)/2+tmp_hig(j,2));
        end
    end
end

% The following sections of code have similar counterparts in Question 1;
% for similicity, I rarely explain what they mean

con_z = repmat(exp(-sig2_u/2+sqrt(sig2_u).*randn(N,1)),1,T);
lab_z = repmat(exp(-sig2_u/2+sqrt(sig2_u).*randn(N,1)),1,T);

con_idio = kron(exp(-sig2_eps/2+sqrt(sig2_eps).*randn(N,T/12)),ones(1,12));
lab_idio = kron(exp(-sig2_eps/2+sqrt(sig2_eps).*randn(N,T/12)),ones(1,12));

% Household consumptions
con_low = con_z.*gm_low.*con_ssc_low.*con_idio;
con_mid = con_z.*gm_mid.*con_ssc_mid.*con_idio;
con_hig = con_z.*gm_hig.*con_ssc_hig.*con_idio;
con_low_0idio = con_z.*gm_low.*con_ssc_low;
con_mid_0idio = con_z.*gm_mid.*con_ssc_mid;
con_hig_0idio = con_z.*gm_hig.*con_ssc_hig;
con_low_0ssc = con_z.*gm_low.*con_idio;
con_mid_0ssc = con_z.*gm_mid.*con_idio;
con_hig_0ssc = con_z.*gm_hig.*con_idio;
con_low_0ssc0idio = con_z.*gm_low;
con_mid_0ssc0idio = con_z.*gm_mid;
con_hig_0ssc0idio = con_z.*gm_hig;
con_0seas = con_z.*con_idio;

con_low_neg1 = con_z.*gm_low.*con_ssc_low_neg.*con_idio;
con_mid_neg1 = con_z.*gm_mid.*con_ssc_mid_neg.*con_idio;
con_hig_neg1 = con_z.*gm_hig.*con_ssc_hig_neg.*con_idio;
con_low_neg2 = con_z.*gm_low_neg.*con_ssc_low_neg.*con_idio;
con_mid_neg2 = con_z.*gm_mid_neg.*con_ssc_mid_neg.*con_idio;
con_hig_neg2 = con_z.*gm_hig_neg.*con_ssc_hig_neg.*con_idio;

% Household labor supplies
lab_low = lab_z.*gm_low.*lab_ssc_low.*lab_idio;
lab_mid = lab_z.*gm_mid.*lab_ssc_mid.*lab_idio;
lab_hig = lab_z.*gm_hig.*lab_ssc_hig.*lab_idio;
lab_low_0idio = lab_z.*gm_low.*lab_ssc_low;
lab_mid_0idio = lab_z.*gm_mid.*lab_ssc_mid;
lab_hig_0idio = lab_z.*gm_hig.*lab_ssc_hig;
lab_low_0ssc = lab_z.*gm_low.*lab_idio;
lab_mid_0ssc = lab_z.*gm_mid.*lab_idio;
lab_hig_0ssc = lab_z.*gm_hig.*lab_idio;
lab_low_0ssc0idio = lab_z.*gm_low;
lab_mid_0ssc0idio = lab_z.*gm_mid;
lab_hig_0ssc0idio = lab_z.*gm_hig;
lab_0seas = lab_z.*lab_idio;

lab_low_neg1 = lab_z.*gm_low_neg.*lab_ssc_low_neg.*lab_idio;
lab_mid_neg1 = lab_z.*gm_mid_neg.*lab_ssc_mid_neg.*lab_idio;
lab_hig_neg1 = lab_z.*gm_hig_neg.*lab_ssc_hig_neg.*lab_idio;
lab_low_neg2 = lab_z.*gm_low.*lab_ssc_low_neg.*lab_idio;
lab_mid_neg2 = lab_z.*gm_mid.*lab_ssc_mid_neg.*lab_idio;
lab_hig_neg2 = lab_z.*gm_hig.*lab_ssc_hig_neg.*lab_idio;

% (a). Assume a deterministic seasonal component and a stochastic seasonal
% component for labor supply both of which are higly positively correlated
% with their consumption counterparts. Then, compute the welfare gains of 
% removing seasons isolating the effects of consumption and leisure.

% Effects of both
for i = 1:N
    total_low(i) = welfaregain(con_low,lab_low,con_0seas,lab_0seas,i);
    total_mid(i) = welfaregain(con_mid,lab_mid,con_0seas,lab_0seas,i);
    total_hig(i) = welfaregain(con_hig,lab_hig,con_0seas,lab_0seas,i);
end

% Effects of consumption
for i = 1:N
    coneff_low(i) = welfaregain(con_low,lab_low,con_0seas,lab_low,i);
    coneff_mid(i) = welfaregain(con_mid,lab_mid,con_0seas,lab_mid,i);
    coneff_hig(i) = welfaregain(con_hig,lab_hig,con_0seas,lab_hig,i);
end

% Effects of labor
for i = 1:N
    labeff_low(i) = welfaregain(con_0seas,lab_low,con_0seas,lab_0seas,i);
    labeff_mid(i) = welfaregain(con_0seas,lab_mid,con_0seas,lab_0seas,i);
    labeff_hig(i) = welfaregain(con_0seas,lab_hig,con_0seas,lab_0seas,i);
end

disp('Part (a) Result')
Item = {'Effects of both';'Effects of consumption';'Effects of labor'};
Low = [mean(total_low);mean(coneff_low);mean(labeff_low)];
Middle = [mean(total_mid);mean(coneff_mid);mean(labeff_mid)];
High = [mean(total_hig);mean(coneff_hig);mean(labeff_hig)];
Parta_Result = table(Item,Low,Middle,High)

% (b). Assume a deterministic seasonal component and a stochastic seasonal 
% component for labor supply both of which are higly negatively correlated
% with their consumption counterparts. Then, compute the welfare gains of 
% removing seasons isolating the effects of consumption and leisure.

% Consider two scenarios: one is when the deterministic component of
% consumption is as original and that of labor is the other way around; the
% other scenario is the opposite situation

% Effects of both
for i = 1:N
    total_low1(i) = ...
        welfaregain(con_low_neg1,lab_low_neg1,con_0seas,lab_0seas,i);
    total_mid1(i) = ...
        welfaregain(con_mid_neg1,lab_mid_neg1,con_0seas,lab_0seas,i);
    total_hig1(i) = ...
        welfaregain(con_hig_neg1,lab_hig_neg1,con_0seas,lab_0seas,i);
    total_low2(i) = ...
        welfaregain(con_low_neg2,lab_low_neg2,con_0seas,lab_0seas,i);
    total_mid2(i) = ...
        welfaregain(con_mid_neg2,lab_mid_neg2,con_0seas,lab_0seas,i);
    total_hig2(i) = ...
        welfaregain(con_hig_neg2,lab_hig_neg2,con_0seas,lab_0seas,i);
end

% Effects of consumption
for i = 1:N
    coneff_low1(i) = ...
        welfaregain(con_low_neg1,lab_low_neg1,con_0seas,lab_low_neg1,i);
    coneff_mid1(i) = ...
        welfaregain(con_mid_neg1,lab_mid_neg1,con_0seas,lab_mid_neg1,i);
    coneff_hig1(i) = ...
        welfaregain(con_hig_neg1,lab_hig_neg1,con_0seas,lab_hig_neg1,i);
    coneff_low2(i) = ...
        welfaregain(con_low_neg2,lab_low_neg2,con_0seas,lab_low_neg2,i);
    coneff_mid2(i) = ...
        welfaregain(con_mid_neg2,lab_mid_neg2,con_0seas,lab_mid_neg2,i);
    coneff_hig2(i) = ...
        welfaregain(con_hig_neg2,lab_hig_neg2,con_0seas,lab_hig_neg2,i);
end

% Effects of labor
for i = 1:N
    labeff_low1(i) = ...
        welfaregain(con_0seas,lab_low_neg1,con_0seas,lab_0seas,i);
    labeff_mid1(i) = ...
        welfaregain(con_0seas,lab_mid_neg1,con_0seas,lab_0seas,i);
    labeff_hig1(i) = ...
        welfaregain(con_0seas,lab_hig_neg1,con_0seas,lab_0seas,i);
    labeff_low2(i) = ...
        welfaregain(con_0seas,lab_low_neg2,con_0seas,lab_0seas,i);
    labeff_mid2(i) = ...
        welfaregain(con_0seas,lab_mid_neg2,con_0seas,lab_0seas,i);
    labeff_hig2(i) = ...
        welfaregain(con_0seas,lab_hig_neg2,con_0seas,lab_0seas,i);
end

disp('Part (b) Result')
Item = {'Effects of both';'Effects of consumption';'Effects of labor'};
Low1 = [mean(total_low1);mean(coneff_low1);mean(labeff_low1)];
Middle1 = [mean(total_mid1);mean(coneff_mid1);mean(labeff_mid1)];
High1 = [mean(total_hig1);mean(coneff_hig1);mean(labeff_hig1)];
Low2 = [mean(total_low2);mean(coneff_low2);mean(labeff_low2)];
Middle2 = [mean(total_mid2);mean(coneff_mid2);mean(labeff_mid2)];
High2 = [mean(total_hig2);mean(coneff_hig2);mean(labeff_hig2)];
Parta_Result = table(Item,Low1,Middle1,High1,Low2,Middle2,High2)

%% The function to compute welfare gains
function g = welfaregain(a,b,c,d,i)
global beta_all kappa nu
tmp = @(g) abs(sum(beta_all.*(log(a(i,:)*(1+g))-kappa*b(i,:).^(1+1/nu)/...
    (1+1/nu))-beta_all.*(log(c(i,:))-kappa*d(i,:).^(1+1/nu)/(1+1/nu))));
g = fminbnd(tmp,-20,20);
end