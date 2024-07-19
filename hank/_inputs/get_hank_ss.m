%% HANK MODEL: STEADY-STATE COMPUTATION
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

%% ECONOMIC PARAMETERS

%----------------------------------------------------------------
% Household
%----------------------------------------------------------------

% preferences

beta      = [];
gamma     = 1;
probdeath = 0;

% wealth lower bound

a_lb = 0;

% income risk

load logyPgrid.txt
load yPtrans.txt
grid_yP = exp(logyPgrid);
n_yP = length(grid_yP);
Pi_yP = yPtrans;
clear logyPgrid yPtrans

load logyTgrid.txt
load yTdist.txt
grid_yT = exp(logyTgrid);
n_yT = length(grid_yT);
Pi_yT = repmat(yTdist',n_yT,1);
clear logyTgrid yTdist

n_y = n_yT * n_yP;
grid_y = repmat(grid_yP,1,n_yT) .* reshape(repmat(grid_yT,1,n_yP)',1,n_y);
indx_yP = reshape(repmat((1:1:n_yP)',1,n_yT),n_y,1)';
indx_yT = reshape(repmat(1:1:n_yT,n_yP,1),n_y,1)';

Pi_y = NaN(n_y,n_y);
for i_y = 1:n_y
    for i_yy = 1:n_y
        Pi_y(i_y,i_yy) = Pi_yP(indx_yP(i_y),indx_yP(i_yy)) * Pi_yT(indx_yT(i_y),indx_yT(i_yy));
    end
end

y_dist  = ergodicdist(Pi_y)'; 
yT_dist = ergodicdist(Pi_yT)';
yP_dist = ergodicdist(Pi_yP)';

mean_y_raw = sum(grid_y(:).*y_dist(:));
grid_y     = grid_y/mean_y_raw;
grid_yP    = grid_yP/mean_y_raw;
mean_y     = grid_y * y_dist';

%----------------------------------------------------------------
% Government
%----------------------------------------------------------------

BY_ratio     = 1.04;

%% STEADY-STATE CALIBRATION

% interest rate

r_b_SS         = 0.01;
annuity_gap    = probdeath/(1-probdeath);

% total output

Z_SS = 1;
Y_SS = 1;

% total debt

B_SS = BY_ratio * Y_SS;

% total government transfers and tax rate

y_tax    = 1/3;
Trans_SS = 0.06;
G_SS     = y_tax * Y_SS + B_SS - (1 + r_b_SS) * B_SS - Trans_SS;

% other aggregates

Pi_SS  = 1;
R_b_SS = 1 + r_b_SS;
R_n_SS = R_b_SS;

% transfer for MPC distribution

gift = 1/15000 * Y_SS;

%% SOLUTION PARAMETERS

%----------------------------------------------------------------
% Grids
%----------------------------------------------------------------

n_a        = 70; % number of grid points 
a_min      = a_lb;
a_max      = 20 * max(grid_yP) * Y_SS;
spliorder  = [1 1];

%----------------------------------------------------------------
% Beta Loop
%----------------------------------------------------------------

beta_init   = 0.99 * 1/(1 + r_b_SS);
beta_guess  = beta_init;
beta_upd    = 0.4;
beta_tol    = 10^(-4);
beta_it_max = 50;
beta_ub     = 1;
beta_lb     = 0.90;

%----------------------------------------------------------------
% EGP Iteration
%----------------------------------------------------------------

EGP_tol      = 10^(-8);
disp_EGPdist = 0;

%% ASSEMBLE GRID

%----------------------------------------------------------------
% Asset Grid
%----------------------------------------------------------------

gridparam_a = 0.3; %linear = 1, 1 L-shaped = 0;

grid_a = linspace(0,1,n_a);
grid_a = grid_a.^(1./gridparam_a);
grid_a = a_min + (a_max-a_min).*grid_a;

wealth_0_pos = find(grid_a == 0);

%----------------------------------------------------------------
% Splines
%----------------------------------------------------------------

% put productivity and asset grids together

n_s  = n_a * n_yP;

fspace          = fundef({'spli',grid_a,0,spliorder(1)},...
                         {'spli',grid_yP,0,spliorder(2)});
states_grid  = funnode(fspace);
states       = gridmake(states_grid);

states_a   = states(:,1);
states_yP  = states(:,2);

Phi_yP     = splibas(grid_yP,0,spliorder(2),states(:,2));
Phi_A      = splibas(grid_a,0,spliorder(1),states(:,1));
Phi        = dprod(Phi_yP,Phi_A);
Emat_yP    = kron(Pi_yP,speye(n_a));

%----------------------------------------------------------------
% Return Grid
%----------------------------------------------------------------

r_b_grid = (states_a >= 0) .* (r_b_SS + annuity_gap) + (states_a < 0) .* (r_b_SS + annuity_gap);

%% MAIN SOLUTION LOOP

for beta_it = 1:beta_it_max
    
%----------------------------------------------------------------
% Discount Facor
%----------------------------------------------------------------

beta     = beta_guess;
beta_hat = beta * (1 - probdeath);  
 
%----------------------------------------------------------------
% Endogenous Gridpoint Iteration
%----------------------------------------------------------------

% preparations

dist_EGP = 1;
EGP_it   = 0;

% initial guess

cp_opt  = (1-y_tax) * Y_SS * (states_yP * grid_yT') + r_b_grid .* states_a + Trans_SS;
mutilde_opt = Emat_yP * (((1 + r_b_grid) .* cp_opt.^(-gamma)) * yT_dist');

% iteration

while dist_EGP > EGP_tol
    
% one step for EGP

[c_opt,ap_opt,mutilde_upd] = EGP_fun(mutilde_opt,r_b_grid,Y_SS,Trans_SS,0,...
    beta_hat,gamma,y_tax,states,grid_a,Emat_yP,grid_yT,yT_dist);
    
% update

dist_EGP = norm(mutilde_upd - mutilde_opt)/norm(mutilde_opt);

if disp_EGPdist == 1 && mod(EGP_it,100) == 0
    disp(dist_EGP)
end

mutilde_opt = mutilde_upd;
EGP_it      = EGP_it + 1;

end

mutilde_SS = mutilde_opt;

%----------------------------------------------------------------
% Distribution
%----------------------------------------------------------------

ap_opt     = max(min(ap_opt,a_max),a_min);
fspaceerga = fundef({'spli',grid_a,0,1});

QZ_live = kron(Pi_yP,ones(n_a,1));
QA_live = 0;
for i_yT = 1:n_yT
    QA_live = QA_live + yT_dist(i_yT) * funbas(fspaceerga,ap_opt(:,i_yT));
end
Q_live  = dprod(QZ_live,QA_live);

QA_death = sparse(n_s,n_a);
QA_death(:,wealth_0_pos) = 1;
QZ_death = repmat(yP_dist,n_s,1);
Q_death  = dprod(QZ_death,QA_death);

Q = (1-probdeath) * Q_live + probdeath * Q_death;

lambda_SS      = ergodicdist(Q,2);
lambda_vec_SS  = lambda_SS;
lambda_SS      = permute(reshape(lambda_SS,[n_a,n_yP]),[2 1]);
lambdafull_SS  = kron(yT_dist',lambda_SS);

%----------------------------------------------------------------
% Compute Aggregates
%----------------------------------------------------------------

% re-order policy function

c_opt_SS  = permute(reshape(c_opt,[n_a,n_y]),[2 1]);
ap_opt_SS = permute(reshape(ap_opt,[n_a,n_y]),[2 1]);

% compute other aggregates

C_grid = c_opt_SS;
C_SS   = full(sum(sum(C_grid(:,:).*lambdafull_SS(:,:))));

A_grid    = repmat(grid_a,n_y,1);
Omega_SS  = sum(sum(A_grid(:,:).*lambdafull_SS(:,:)));

%----------------------------------------------------------------
% Check Asset Market Clearing
%----------------------------------------------------------------

Wealth_err = Omega_SS - B_SS;

% ----------------------------------------------------------------
% Beta Update
% ----------------------------------------------------------------

if Wealth_err < -beta_tol
    disp(['Discount factor: ' num2str(beta_guess) ', Discount factor too low: ' num2str(Wealth_err) ]);
    beta_lb = beta_guess;
    beta_guess = (1-beta_upd) * beta_guess + beta_upd * beta_ub;
elseif Wealth_err > beta_tol
    disp(['Discount factor: ' num2str(beta_guess) ', Discount factor too high: ' num2str(Wealth_err) ]);
    beta_ub = beta_guess;
    beta_guess = (1-beta_upd) * beta_guess + beta_upd * beta_lb;
elseif abs(Wealth_err) <= beta_tol
    disp(['Steady State Found, Discount factor = ' num2str(beta_guess)]);
    break
end

end

%% SAVE RESULTS

%----------------------------------------------------------------
% Aggregate Parameters
%----------------------------------------------------------------

save param_agg beta beta_hat gamma probdeath wealth_0_pos ...
     y_tax BY_ratio

%----------------------------------------------------------------
% Household Parameters
%----------------------------------------------------------------

save param_households a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap

%----------------------------------------------------------------
% Steady State
%----------------------------------------------------------------

save SS C_SS Y_SS Trans_SS Pi_SS R_n_SS R_b_SS B_SS ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

%----------------------------------------------------------------
% Other Quantities
%----------------------------------------------------------------

save aux grid_a spliorder states states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max