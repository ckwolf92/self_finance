%% GET 3-TYPE OLG AGGREGATE CONSUMPTION FUNCTION
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
addpath([path experiment '/_impcs']);
cd([path experiment]);

%% PARAMETERS

%----------------------------------------------------------------
% Calibration Target
%----------------------------------------------------------------

global MPC_target

load fhn_results

MPC_target = cumC1(5:10)';

%----------------------------------------------------------------
% Preferences
%----------------------------------------------------------------

global beta omega omega_1 omega_2 sigma MPC_target chi_1 chi_2

beta  = 0.99^(0.25);

calib_fn = @(param) calib_fn_aux(param);

param_sol = fminsearch(calib_fn,[0.8 0.6 0.3 0.5]);

omega_1 = param_sol(1);
omega_2 = param_sol(2);
chi_1   = param_sol(3);
chi_2   = param_sol(4);

sigma   = 1;

%----------------------------------------------------------------
% Steady-State Objects
%----------------------------------------------------------------

global tau_y Y_SS C_SS r_SS D_SS D_OLG_SS

tau_y    = 1/3;
Y_SS     = 1;
C_SS     = Y_SS;
r_SS     = 1/beta - 1;
D_SS     = 1.04;
D_OLG_SS = 1/(chi_1 + chi_2) * D_SS;

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

%% DERIVATIVE MATRICES: TYPE 1

%----------------------------------------------------------------
% Auxiliary Matrix
%----------------------------------------------------------------

omega = omega_1;

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

[c_base,d_base] = c_hybrid_fn(exo,A_inv);

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

C_y_1 = NaN(T,T);
D_y_1 = NaN(T,T);

for t = 1:T
    exo.y_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_y_1(:,t) = (c_shock-c_base)/step;
    D_y_1(:,t) = (d_shock-d_base)/step;
    exo.y_hat(t) = 0;
end

%----------------------------------------------------------------
% Nominal Interest Rates
%----------------------------------------------------------------

C_i_1 = NaN(T,T);
D_i_1 = NaN(T,T);

for t = 1:T
    exo.i_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_i_1(:,t) = (c_shock-c_base)/step;
    D_i_1(:,t) = (d_shock-d_base)/step;
    exo.i_hat(t) = 0;
end

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_pi_1 = NaN(T,T);
D_pi_1 = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_pi_1(:,t) = (c_shock-c_base)/step;
    D_pi_1(:,t) = (d_shock-d_base)/step;
    exo.pi_hat(t) = 0;
end

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_d_1 = NaN(T,T);
D_d_1 = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_d_1(:,t) = (c_shock-c_base)/step;
    D_d_1(:,t) = (d_shock-d_base)/step;
    exo.zeta_hat(t) = 0;
end

%% DERIVATIVE MATRICES: TYPE 2

%----------------------------------------------------------------
% Auxiliary Matrix
%----------------------------------------------------------------

omega = omega_2;

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

[c_base,d_base] = c_hybrid_fn(exo,A_inv);

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

C_y_2 = NaN(T,T);
D_y_2 = NaN(T,T);

for t = 1:T
    exo.y_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_y_2(:,t) = (c_shock-c_base)/step;
    D_y_2(:,t) = (d_shock-d_base)/step;
    exo.y_hat(t) = 0;
end

%----------------------------------------------------------------
% Nominal Interest Rates
%----------------------------------------------------------------

C_i_2 = NaN(T,T);
D_i_2 = NaN(T,T);

for t = 1:T
    exo.i_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_i_2(:,t) = (c_shock-c_base)/step;
    D_i_2(:,t) = (d_shock-d_base)/step;
    exo.i_hat(t) = 0;
end

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_pi_2 = NaN(T,T);
D_pi_2 = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_pi_2(:,t) = (c_shock-c_base)/step;
    D_pi_2(:,t) = (d_shock-d_base)/step;
    exo.pi_hat(t) = 0;
end

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_d_2 = NaN(T,T);
D_d_2 = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_d_2(:,t) = (c_shock-c_base)/step;
    D_d_2(:,t) = (d_shock-d_base)/step;
    exo.zeta_hat(t) = 0;
end

%% DERIVATIVE MATRICES: AGGREGATE

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

C_y = chi_1 * C_y_1 + chi_2 * C_y_2 + (1 - chi_1 - chi_2) * eye(T);
D_y = chi_1 * D_y_1 + chi_2 * D_y_2 + (1 - chi_1 - chi_2) * zeros(T,T);

%----------------------------------------------------------------
% Nominal Interest Rates
%----------------------------------------------------------------

C_i = chi_1 * C_i_1 + chi_2 * C_i_2 + (1 - chi_1 - chi_2) * zeros(T,T);
D_i = chi_1 * D_i_1 + chi_2 * D_i_2 + (1 - chi_1 - chi_2) * zeros(T,T);

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_pi = chi_1 * C_pi_1 + chi_2 * C_pi_2 + (1 - chi_1 - chi_2) * zeros(T,T);
D_pi = chi_1 * D_pi_1 + chi_2 * D_pi_2 + (1 - chi_1 - chi_2) * zeros(T,T);

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_d = chi_1 * C_d_1 + chi_2 * C_d_2 + (1 - chi_1 - chi_2) * zeros(T,T);
D_d = chi_1 * D_d_1 + chi_2 * D_d_2 + (1 - chi_1 - chi_2) * zeros(T,T);

%% SAVE RESULTS

cd([path experiment '/_results']);

save inputs_3type C_y C_i C_pi C_d D_y D_i D_pi D_d ...
    C_y_1 C_y_2 C_i_1 C_i_2 C_pi_1 C_pi_2 C_d_1 C_d_2 D_y_1 D_y_2 D_i_1 D_i_2 D_pi_1 D_pi_2 D_d_1 D_d_2 ...
    beta omega_1 omega_2 chi_1 chi_2 sigma ...
    r_SS tau_y Y_SS C_SS D_SS

cd([path experiment]);