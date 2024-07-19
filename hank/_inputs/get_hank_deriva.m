%% HANK MODEL: AGGREGATE CONSUMPTION FUNCTION MATRICES
% Marios Angeletos, Chen Lian, & Christian Wolf
% this version: 06/17/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Dropbox/Research/self_finance/codes/ckw';
path = [local '/ecma_replication/hank'];
experiment = '/_inputs';

addpath(genpath([path experiment '/_auxiliary_functions']));
addpath([path experiment '/_income_process'])
cd([path experiment]);

%% PREPARATIONS

%----------------------------------------------------------------
% Global Variables
%----------------------------------------------------------------

% aggregate parameters

global beta beta_hat gamma probdeath wealth_0_pos ...
     y_tax BY_ratio
 
% household parameters

global a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist mu Pi_y Pi_yP Pi_yT annuity_gap

% steady-state quantities

global C_SS Y_SS Trans_SS Pi_SS R_n_SS R_b_SS B_SS ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

% other quantities

global grid_a spliorder states states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

load param_agg
load param_households
load SS
load aux

%----------------------------------------------------------------
% Re-definitions
%----------------------------------------------------------------

r_b_SS = R_b_SS - 1;

%----------------------------------------------------------------
% Time Horizon
%----------------------------------------------------------------

global T

T = 500;

%% GET ALL M MATRICES

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

global step

step = 10^(-3);

%----------------------------------------------------------------
% Prep Run
%----------------------------------------------------------------

global c_opt_vec_SS Qt_big_SS

y_seq     = zeros(T,1);
ib_seq    = zeros(T,1);
pi_seq    = zeros(T,1);
tau_x_seq = zeros(T,1);
tau_e_seq = zeros(T,1);
eps_b_seq = zeros(T,1);
    
solve_hh_problem

c_opt_vec_SS = c_opt_t(:,:,end);

Qt_SS     = Qt;
Qt_big_SS = NaN(n_y*n_a,n_y*n_a);

for i_yT = 1:n_yT
    Qt_big_SS(:,1+(i_yT-1)*(n_a*n_yP):i_yT*(n_a*n_yP)) = repmat(Qt_SS,n_yT,1) * yT_dist(i_yT);
end

%----------------------------------------------------------------
% Output
%----------------------------------------------------------------

vars = zeros(6,1);
vars(1) = step;

[Y_c_y,DD_y] = YD_fn(vars);

[C_y,D_y] = jac_fn(Y_c_y,DD_y);

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

vars = zeros(6,1);
vars(3) = step;

[Y_c_pi,DD_pi] = YD_fn(vars);

[C_pi,D_pi] = jac_fn(Y_c_pi,DD_pi);

%----------------------------------------------------------------
% Nominal Interest Rate
%----------------------------------------------------------------

C_i   = [-C_pi(:,2:T),[0;-C_pi(1:T-1,T)]];
D_i   = [-D_pi(:,2:T),[0;-D_pi(1:T-1,T)]];

%----------------------------------------------------------------
% Transfer
%----------------------------------------------------------------

vars = zeros(6,1);
vars(4) = step;

[Y_c_tau,DD_tau] = YD_fn(vars);

[C_tau,D_tau] = jac_fn(Y_c_tau,DD_tau);

M = C_tau * C_SS/Trans_SS;

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_d = zeros(T,T);
D_d = zeros(T,T);

%% ADJUST NOTATION

D_SS  = B_SS;
tau_y = y_tax;
r_SS  = r_b_SS;

%% SAVE RESULTS

save inputs_hank C_i C_pi C_y C_tau C_d ...
    D_i D_pi D_y D_tau D_d ...
    gamma beta ...
    r_SS tau_y Y_SS C_SS D_SS Trans_SS