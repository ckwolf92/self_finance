%% LARGEST DEGREE OF SELF-FINANCING WITH TAYLOR RULE
% Marios Angeletos, Chen Lian, & Christian Wolf
% this version: 06/17/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Dropbox/Research/self_finance/codes/ckw';
path = [local '/ecma_replication/analytical'];
experiment = '/tau_activemp';

addpath([path '/_auxiliary_functions']);
addpath([path '/_inputs/_results']);
addpath([path experiment '/_aux_ge']);
cd([path experiment]);

%% PARAMETERS

%----------------------------------------------------------------
% Household Block
%----------------------------------------------------------------

global tau_y beta r_SS D_SS Y_SS D_i D_pi D_y C_i C_pi C_y

load inputs_hybrid

global T

T = size(C_y,1);

%----------------------------------------------------------------
% Price Stickiness
%----------------------------------------------------------------

global kappa

kappa_grid = [0.0062, 2*0.0062, 0.1];
n_kappa    = length(kappa_grid);

%----------------------------------------------------------------
% NPV-Maker
%----------------------------------------------------------------

global r_NPV

r_NPV = NaN(T,1);
for t = 1:T
    r_NPV(t) = (1/(1 + r_SS))^(t-1);
end

%----------------------------------------------------------------
% Policy
%----------------------------------------------------------------

global phi_pi tau_d H tau_x_seq

% monetary rule

phi_pi_grid = [1, 1.25, 1.5];
n_phi_pi    = length(phi_pi_grid);

% fiscal rule

n_tau_d    = 50;
tau_d_list = linspace(0,0.1,n_tau_d);
H          = 1;

% policy shock

tau_x_seq = zeros(T,1);
tau_x_seq(1) = -1;

%% COMPUTE LARGEST SHARE OF SELF-FINANCING

%----------------------------------------------------------------
% Outer Loop
%----------------------------------------------------------------

nu_max_grid = NaN(n_kappa,n_phi_pi);

for i_kappa = 1:n_kappa
    
disp(['I am at step ' num2str(i_kappa) ' of the outer loop.'])
    
for i_phi_pi = 1:n_phi_pi
    
kappa  = kappa_grid(i_kappa);
phi_pi = phi_pi_grid(i_phi_pi);

%----------------------------------------------------------------
% GE Loop
%----------------------------------------------------------------

temp_nu = NaN(n_tau_d,1);

for i_tau_d = 1:n_tau_d

tau_d = tau_d_list(i_tau_d);

% initial wedge

excess_demand_init = excess_demand_fn_tr(zeros(T,1));

% GE updating

step = 1;
A_upd = NaN(T,T);

for i_deriv = 1:T
    guess_seq_deriv = zeros(T,1);
    guess_seq_deriv(i_deriv,1) = guess_seq_deriv(i_deriv,1) + step;
    excess_demand_up = excess_demand_fn_tr(guess_seq_deriv);
    A_upd(:,i_deriv) = (excess_demand_up-excess_demand_init)/step;
end

sol_seq = -A_upd^(-1) * excess_demand_init;

% GE transition

y_seq = sol_seq(1:T);

if abs(y_seq(T)) < 1e-8

    get_aggregates_tr
    
    expenditure = (r_NPV' * tau_x_seq - r_NPV(2:T)' * D_SS * (i_seq(1:T-1) - pi_seq(2:T)));

    nu_base   = (-r_NPV' * Y_SS * tau_y * y_seq) ./ expenditure;
    nu_prices = (- D_SS * pi_seq(1)) ./ expenditure;

    temp_nu(i_tau_d,1) = nu_base + nu_prices;

else

    temp_nu(i_tau_d,1) = NaN;
    
end

end

nu_max_grid(i_kappa,i_phi_pi) = max(temp_nu);

end

end

%% PRINT TABLE

save_filename = 'maxnu_table';

f = fopen(fullfile('_results', strcat(save_filename, '.tex')), 'w'); % open file for writing

fprintf(f, '%s%s%s%s%s%s%s\n', '\begin{tabular}{r|', repmat('c', 1, 1), '', repmat('c', 1, 1), '', repmat('c', 1, 1), '}');
fprintf(f, '%s%d%s%d%s%d%s\n', '$\kappa$ & \multicolumn{', 1, '}{c}{$\psi = 1$} & \multicolumn{', 1, '}{c}{$\psi = 1.25$} & \multicolumn{', 1, '}{c}{$\psi = 1.5$} \\');
fprintf(f, '%s\n%s\n', '\hline');

for i_kappa = 1:n_kappa
    kappa = kappa_grid(i_kappa);
    fprintf(f, '%4.2f', kappa);
    phi_pi_1 = phi_pi_grid(1);
    fprintf(f, '%s%5.3f', ' & ', nu_max_grid(i_kappa,1));
    phi_pi_2 = phi_pi_grid(2);
    fprintf(f, '%s%5.3f', ' & ', nu_max_grid(i_kappa,2));
    phi_pi_3 = phi_pi_grid(3);
    fprintf(f, '%s%5.3f', ' & ', nu_max_grid(i_kappa,3));
    fprintf(f, '%s\n', ' \\');
end

fprintf(f, '%s', '\end{tabular}');

fclose(f);