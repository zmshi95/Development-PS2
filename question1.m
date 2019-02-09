% Growth and Devleopment 
% Problem set 2 
%This code is written by Zhongming Shi

%Question 1. Praying for Rain: The welfare cost of seasons

%initialize parameters
clear all;
rng(2019)
clc


sigma_u=0.2;
sigma_e=0.2;
n=1000;
beta=.99^(1/12);
T=40;
eta=[1,2,4];

%import data from table 1
table1=importdata('table1.txt');
%middle=table1(:,1);
%high=table1(:,2);
%low=table1(:,3);

%draw z for n's households
z=randn(n,1)*sigma_u; % generate lnU
z=exp(z)*exp(-0.5*(sigma_e^2)); % generate z

%draw epsilon_t for different age 
e=randn(T,1)*sigma_e; % generate lne
e=exp(e)*exp(-0.5*(sigma_e^2)); % generate e

%compute Cm,t for each individual
c(3,n,T,12)=zeros; %C(middle/high/low;individual;age;month)
for i=1:3;
    for j=1:n
    c(i,j,:,:)=z(j)*(e*table1(:,i)');%i=1: middle i=2: high; i=3: low
    end
end

u(length(eta),3,n,T,12)=zeros; % u is 5-D and first dimension is for different eta
for k=1:length(eta)
    ita=eta(k);
    if ita==1
        u(k,:,:,:,:)=log(c);
    else
        u(k,:,:,:,:)=c.^(1-ita)./(1-ita);
    end
end

%compute W(z) for each individual 
w(length(eta),3,n)=zeros;
beta1=beta.^(1:12);
beta2=beta^12.^(1:40);

for k=1:length(eta)
    eta_is=eta(k)
    for s=1:3
        season_is=s
        for i=1:n
             a=reshape(u(k,s,i,:,:),T,12);
             w(k,s,i)=beta2*a*beta1';
        end         
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%remove seasonal component%%%%%%%%%%%%%%%%%%%
c_season(3,n,T,12)=zeros; %C(middle/high/low;individual;age;month)
for i=1:3;
    for j=1:n
    c_season(i,j,:,:)=z(j)*(e*ones(12,1)');%i=1: middle i=2: high; i=3: low
    end
end

u_season(length(eta),3,n,T,12)=zeros; % u is 5-D and first dimension is for different eta
for k=1:length(eta)
    ita=eta(k);
    if ita==1
        u_season(k,:,:,:,:)=log(c_season);
    else
        u_season(k,:,:,:,:)=c_season.^(1-ita)./(1-ita);
    end
end

%compute W(z) for each individual 
w_season(length(eta),3,n)=zeros;
beta1=beta.^(1:12);
beta2=beta^12.^(1:40);

for k=1:length(eta)
    eta_is=eta(k)
    for s=1:3
        season_is=s
        for i=1:n
             a=reshape(u_season(k,s,i,:,:),T,12);
             w_season(k,s,i)=beta2*a*beta1';
        end         
    end
end

%remove nonseasonal component
c_nonseason(3,n,T,12)=zeros; %C(middle/high/low;individual;age;month)
for i=1:3;
    for j=1:n
    c_nonseason(i,j,:,:)=z(j)*(ones(40,1)*table1(:,i)');%i=1: middle i=2: high; i=3: low
    end
end

u_nonseason(length(eta),3,n,T,12)=zeros; % u is 5-D and first dimension is for different eta
for k=1:length(eta)
    ita=eta(k);
    if ita==1
        u_nonseason(k,:,:,:,:)=log(c_nonseason);
    else
        u_nonseason(k,:,:,:,:)=c_nonseason.^(1-ita)./(1-ita);
    end
end

%compute W(z) for each individual 
w_nonseason(length(eta),3,n)=zeros;
beta1=beta.^(1:12);
beta2=beta^12.^(1:40);

for k=1:length(eta)
    eta_is=eta(k)
    for s=1:3
        season_is=s
        for i=1:n
             a=reshape(u_nonseason(k,s,i,:,:),T,12);
             w_nonseason(k,s,i)=beta2*a*beta1';
        end         
    end

end

%%Now display the result
%(a)
v=sum(w,3);
v_season=sum(w_season,3);
v_nonseason=sum(w_nonseason,3);

for k=1:3
    for s=1:3
    disp(['When eta is ', num2str(eta(k)), ', season is (middle=1,high=2,low=3)',num2str(s),', remove seasonal component creates,',num2str(v_season(k,s)-v(k,s)),' aggregate welfare'])
    disp(['When eta is ', num2str(eta(k)), ', season is (middle=1,high=2,low=3)',num2str(s),', remove nonseasonal component creates,',num2str(v_nonseason(k,s)-v(k,s)),' aggregate welfare'])
    end
end

%from the results we can see, when eta=1, remove nonseasonal component will
%decrease utility but remove seasonal one will increase social welfare
%(which also depends on the seasonality).
%As eta increases from 1 to 4, the social welfare increases more from the
%removal.

