%% SELF-FINANCING SHARE FOR SMALL MARGINS OF NEAR-PIH HOUSEHOLDS
% Marios Angeletos, Chen Lian, & Christian Wolf
% this version: 06/17/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Dropbox/Research/self_finance/codes/ckw';
path = [local '/ecma_replication/analytical'];
experiment = '/tau_base';

addpath([path '/_auxiliary_functions']);
addpath([path '/_inputs/_results']);
addpath([path experiment '/_aux_ge']);
addpath([path experiment '/_aux_jacs']);
cd([path experiment]);

%% PARAMETERS

%----------------------------------------------------------------
% Household Block
%----------------------------------------------------------------

% preferences

global beta omega omega_1 omega_2 sigma chi

beta    = 0.99^(0.25);
omega_1 = 0.75;
chi     = 0.01;
sigma   = 1;

n_omega_2    = 50;
omega_2_grid = linspace(0.95,1,n_omega_2);

% steady-state objects

global tau_y Y_SS C_SS r_SS D_SS

tau_y    = 1/3;
Y_SS     = 1;
C_SS     = Y_SS;
r_SS     = 1/beta - 1;
D_SS     = 1.04;

% more settings & inputs

global T

T    = 500;
step = 1;

exo.y_hat     = zeros(T,1);
exo.i_hat     = zeros(T,1);
exo.pi_hat    = zeros(T,1);
exo.zeta_hat  = zeros(T,1);

% type 1 matrices

get_omega_1_jacs

% global variables

global D_i D_pi D_y C_i C_pi C_y

%----------------------------------------------------------------
% Price Stickiness
%----------------------------------------------------------------

global kappa

kappa = 0.1;

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

global phi tau_d H tau_x_seq

% monetary rule

phi = 0;

% fiscal rule

tau_d  = [];
H_list = [28,60,100];
n_H    = length(H_list);

% policy shock

tau_x_seq = zeros(T,1);
tau_x_seq(1) = -1;

%% COMPUTE GE TRANSITION

%----------------------------------------------------------------
% Placeholders
%----------------------------------------------------------------

y_seq_all      = NaN(T,n_omega_2,n_H);
d_seq_all      = NaN(T,n_omega_2,n_H);
pi_seq_all     = NaN(T,n_omega_2,n_H);
t_seq_all      = NaN(T,n_omega_2,n_H);
nu_base_all    = NaN(n_omega_2,n_H);
nu_rates_all   = NaN(n_omega_2,n_H);
nu_prices_all  = NaN(n_omega_2,n_H);
nu_all         = NaN(n_omega_2,n_H);

%----------------------------------------------------------------
% Loop
%----------------------------------------------------------------

for i_H = 1:n_H

H = H_list(i_H);

disp(['I am at run ' num2str(i_H)])

for i_omega_2 = 1:n_omega_2

% get relevant matrices

omega_2 = omega_2_grid(i_omega_2);

get_omega_2_jacs

C_y = (1-chi) * C_y_1 + chi * C_y_2;
D_y = (1-chi) * D_y_1 + chi * D_y_2;

C_i = (1-chi) * C_i_1 + chi * C_i_2;
D_i = (1-chi) * D_i_1 + chi * D_i_2;

C_pi = (1-chi) * C_pi_1 + chi * C_pi_2;
D_pi = (1-chi) * D_pi_1 + chi * D_pi_2;

C_d = (1-chi) * C_d_1 + chi * C_d_2;
D_d = (1-chi) * D_d_1 + chi * D_d_2;

% initial wedge

excess_demand_init = excess_demand_fn(zeros(T,1));

% GE updating

step = 1;
A_upd = NaN(T,T);

for i_deriv = 1:T
    guess_seq_deriv = zeros(T,1);
    guess_seq_deriv(i_deriv,1) = guess_seq_deriv(i_deriv,1) + step;
    excess_demand_up = excess_demand_fn(guess_seq_deriv);
    A_upd(:,i_deriv) = (excess_demand_up-excess_demand_init)/step;
end

sol_seq = -A_upd^(-1) * excess_demand_init;

% collect results

y_seq = sol_seq(1:T);

get_aggregates

% save results

y_seq_all(:,i_omega_2,i_H)     = y_seq;
d_seq_all(:,i_omega_2,i_H)     = d_seq;
pi_seq_all(:,i_omega_2,i_H)    = pi_seq;
t_seq_all(:,i_omega_2,i_H)     = t_seq;

expenditure = -(r_NPV' * tau_x_seq);

nu_base_all(i_omega_2,i_H)   = (r_NPV' * Y_SS * tau_y * y_seq) ./ expenditure;
nu_rates_all(i_omega_2,i_H)  = (-r_NPV(2:T)' * D_SS * (i_seq(1:T-1) - pi_seq(2:T))) ./ expenditure;
nu_prices_all(i_omega_2,i_H) = (D_SS * pi_seq(1)) ./ expenditure;
nu_all(i_omega_2,i_H)        = nu_base_all(i_omega_2,i_H) + nu_rates_all(i_omega_2,i_H) + nu_prices_all(i_omega_2,i_H);

end

end

%% PLOT RESULTS

%----------------------------------------------------------------
% Color Preparation
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.lpurple = 0.25 * [160/255 32/255 240/255] + 0.75 * [1 1 1];
settings.colors.lorange = 0.25 * [204/255 102/255 0/255] + 0.75 * [1 1 1];

settings.colors.list = [0 * [0.5 0.5 0.5]; 0.5 * [0.5 0.5 0.5]; 1 * [0.5 0.5 0.5]; 1.5 * [0.5 0.5 0.5]; 1.8 * [0.5 0.5 0.5]];

settings.colors.finance = [0.25 * [196/255 174/255 120/255] + 0.75 * [1 1 1]; ...
                            0.5 * [0 0 0] + 0.5 * [1 1 1]; ...
                            0.25 * [189/255 142/255 131/255] + 0.75 * [1 1 1]];

%----------------------------------------------------------------
% Plot
%----------------------------------------------------------------

cd([path experiment '/_results']);

figure(1)
pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(omega_2_grid,nu_all(:,1),'linewidth',5,'linestyle','-','color',settings.colors.list(2,:))
hold on
plot(omega_2_grid,nu_all(:,2),'linewidth',5,'linestyle','-','color',settings.colors.list(3,:))
hold on
plot(omega_2_grid,nu_all(:,3),'linewidth',5,'linestyle','-','color',settings.colors.list(4,:))
hold on
set(gcf,'color','w')
xlabel('Near-PIH $\bar{\omega}$','interpreter','latex','FontSize',20)
ylabel('$\nu$','interpreter','latex','FontSize',20,'Rotation',0)
xlim([min(omega_2_grid),max(omega_2_grid)])
yticks([0 0.2 0.4 0.6 0.8 1])
legend({'$H = 28$','$H = 60$','$H = 100$'},'Location','Southwest','fontsize',18,'interpreter','latex')
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_b3','-dpng');

cd([path experiment]);