%% GET OLG AGGREGATE CONSUMPTION FUNCTION
% Marios Angeletos, Chen Lian, & Christian Wolf
% this version: 06/17/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Dropbox/Research/self_finance/codes/ckw';
path = [local '/ecma_replication/analytical'];
experiment = '/_inputs';

addpath([path '/_auxiliary_functions']);
cd([path experiment]);

%% PARAMETERS

%----------------------------------------------------------------
% Preferences
%----------------------------------------------------------------

global beta omega sigma

beta  = 0.99^(0.25);
omega = 0.75;

sigma = 1;

%----------------------------------------------------------------
% Steady-State Objects
%----------------------------------------------------------------

global tau_y Y_SS C_SS r_SS D_SS

tau_y    = 1/3;
Y_SS     = 1;
C_SS     = Y_SS;
r_SS     = 1/beta - 1;
D_SS     = 1.04;

%% SETTINGS

global T

T    = 500;
step = 1;

exo.y_hat     = zeros(T,1);
exo.i_hat     = zeros(T,1);
exo.pi_hat    = zeros(T,1);
exo.zeta_hat  = zeros(T,1);

r_NPV = zeros(T,1);
for t = 1:T
    r_NPV(t) = (1/(1+r_SS))^(t-1);
end

%% DERIVATIVE MATRICES

%----------------------------------------------------------------
% Auxiliary Matrix
%----------------------------------------------------------------

A = zeros(2*T,2*T); % order (c,d) & (BC, EE)
for t = 1:T
    A(t,t) = 1;
    A(t,T+t) = 1;
    if t > 1
        A(t,T+t-1) = -1/beta;
    end
end
for t = T+1:2*T
    if t == T+1
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-T+1) = - beta * omega;
    elseif t < 2*T
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-T+1) = - beta * omega;
        A(t,t-1) = - (1-beta*omega) * (1-omega) * 1/beta;
    elseif t == 2*T
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-1) = 1;
    end
end

A_inv = A^(-1);

%----------------------------------------------------------------
% Baseline
%----------------------------------------------------------------

[c_base,d_base] = c_olg_fn(exo,A_inv);

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

C_y = NaN(T,T);
D_y = NaN(T,T);

for t = 1:T
    exo.y_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_inv);
    C_y(:,t) = (c_shock-c_base)/step;
    D_y(:,t) = (d_shock-d_base)/step;
    exo.y_hat(t) = 0;
end

%----------------------------------------------------------------
% Nominal Interest Rates
%----------------------------------------------------------------

C_i = NaN(T,T);
D_i = NaN(T,T);

for t = 1:T
    exo.i_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_inv);
    C_i(:,t) = (c_shock-c_base)/step;
    D_i(:,t) = (d_shock-d_base)/step;
    exo.i_hat(t) = 0;
end

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_pi = NaN(T,T);
D_pi = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_inv);
    C_pi(:,t) = (c_shock-c_base)/step;
    D_pi(:,t) = (d_shock-d_base)/step;
    exo.pi_hat(t) = 0;
end

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_d = NaN(T,T);
D_d = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_inv);
    C_d(:,t) = (c_shock-c_base)/step;
    D_d(:,t) = (d_shock-d_base)/step;
    exo.zeta_hat(t) = 0;
end

%% SAVE RESULTS

cd([path experiment '/_results']);

save inputs_olg C_y C_i C_pi C_d D_y D_i D_pi D_d ...
    beta omega sigma ...
    r_SS tau_y Y_SS C_SS D_SS

cd([path experiment]);