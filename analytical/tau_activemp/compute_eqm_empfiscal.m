%% TRANSFER SHOCK IRFs FOR EMPIRICALLY RELEVANT FISCAL RULES + ACTIVE MP
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

kappa = 0.0062;

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

global phi_pi tau_d H tau_x_seq m_seq

% monetary rule

phi_list = [1, 1.25, 1.5];
n_phi    = length(phi_list);

% fiscal rule

n_taud = 60;

tau_d_list_1 = linspace(1,0,n_taud);
tau_d_list_2 = linspace(1,0,n_taud);
tau_d_list_3 = linspace(1,0,n_taud);

H = 1;

% policy shocks

tau_x_seq = zeros(T,1);
tau_x_seq(1) = -1;

m_seq = zeros(T,1);

%% COMPUTE GE TRANSITION

%----------------------------------------------------------------
% Placeholders
%----------------------------------------------------------------

y_seq_all      = NaN(T,n_taud,n_phi);
d_seq_all      = NaN(T,n_taud,n_phi);
pi_seq_all     = NaN(T,n_taud,n_phi);
t_seq_all      = NaN(T,n_taud,n_phi);
i_seq_all      = NaN(T,n_taud,n_phi);
nu_base_all    = NaN(n_taud,n_phi);
nu_rates_all   = NaN(n_taud,n_phi);
nu_prices_all  = NaN(n_taud,n_phi);
nu_all         = NaN(n_taud,n_phi);

%----------------------------------------------------------------
% Loop
%----------------------------------------------------------------

for i_phi = 1:n_phi

disp(['I am at run ' num2str(i_phi)])

phi_pi = phi_list(i_phi);

for i_taud = 1:n_taud

if i_phi == 1

    tau_d = tau_d_list_1(i_taud);

elseif i_phi == 2

    tau_d = tau_d_list_2(i_taud);

elseif i_phi == 3

    tau_d = tau_d_list_3(i_taud);

end

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

% collect results

y_seq = sol_seq(1:T);

get_aggregates_tr

% save results

y_seq_all(:,i_taud,i_phi)     = y_seq;
d_seq_all(:,i_taud,i_phi)     = d_seq;
pi_seq_all(:,i_taud,i_phi)    = pi_seq;
t_seq_all(:,i_taud,i_phi)     = t_seq;
i_seq_all(:,i_taud,i_phi)     = i_seq;

expenditure = (r_NPV' * tau_x_seq - r_NPV(2:T)' * D_SS * (i_seq(1:T-1) - pi_seq(2:T)));

nu_base_all(i_taud,i_phi)   = (-r_NPV' * Y_SS * tau_y * y_seq) ./ expenditure;
nu_prices_all(i_taud,i_phi) = (- D_SS * pi_seq(1)) ./ expenditure;

nu_all(i_taud,i_phi)        = nu_base_all(i_taud,i_phi) + nu_prices_all(i_taud,i_phi);

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
settings.colors.lred = 0.25 * [1 0 0] + 0.75 * [1 1 1];

settings.colors.list = [0 * [0.5 0.5 0.5]; 0.5 * [0.5 0.5 0.5]; 1 * [0.5 0.5 0.5]; 1.5 * [0.5 0.5 0.5]; 1.8 * [0.5 0.5 0.5]];

settings.colors.finance = [0.25 * [196/255 174/255 120/255] + 0.75 * [1 1 1]; ...
                            0.5 * [0 0 0] + 0.5 * [1 1 1]; ...
                            0.25 * [189/255 142/255 131/255] + 0.75 * [1 1 1]];

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

IRF_plot = 40;

%----------------------------------------------------------------
% Plot
%----------------------------------------------------------------

cd([path experiment '/_results']);

plotwidth = 0.25;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * (gapsize + plotwidth)];

figure(1)

subplot(1,3,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
jbfill(tau_d_list_1,(squeeze(0*nu_prices_all(:,1)))',(squeeze(nu_prices_all(:,1)))',settings.colors.finance(1,:),settings.colors.finance(1,:),0,1);
jbfill(tau_d_list_1,(squeeze(nu_prices_all(:,1)))',(squeeze(nu_prices_all(:,1)+nu_base_all(:,1)))',settings.colors.finance(2,:),settings.colors.finance(2,:),0,1);
hold on
set(gcf,'color','w')
title('$\psi = 1$','interpreter','latex','fontsize',24)
xlabel('$\tau_d$','interpreter','latex','FontSize',20)
ylabel('$\nu$','interpreter','latex','FontSize',20,'rotation',0)
xlim([0 1])
ylim([0 1])
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
jbfill(tau_d_list_2,(squeeze(0*nu_prices_all(:,2)))',(squeeze(nu_prices_all(:,2)))',settings.colors.finance(1,:),settings.colors.finance(1,:),0,1);
jbfill(tau_d_list_2,(squeeze(nu_prices_all(:,2)))',(squeeze(nu_prices_all(:,2)+nu_base_all(:,2)))',settings.colors.finance(2,:),settings.colors.finance(2,:),0,1);
hold on
jbfill([0,tau_d_list_2(end)],[0 0],[1 1],settings.colors.lred,settings.colors.lred,0,1);
hold on
set(gcf,'color','w')
title('$\psi = 1.25$','interpreter','latex','fontsize',24)
xlabel('$\tau_d$','interpreter','latex','FontSize',20)
xlim([0 1])
ylim([0 1])
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
jbfill(tau_d_list_3,(squeeze(0*nu_prices_all(:,3)))',(squeeze(nu_prices_all(:,3)))',settings.colors.finance(1,:),settings.colors.finance(1,:),0,1);
jbfill(tau_d_list_3,(squeeze(nu_prices_all(:,3)))',(squeeze(nu_prices_all(:,3)+nu_base_all(:,3)))',settings.colors.finance(2,:),settings.colors.finance(2,:),0,1);
hold on
jbfill([0,tau_d_list_3(end)],[0 0],[1 1],settings.colors.lred,settings.colors.lred,0,1);
hold on
set(gcf,'color','w')
title('$\psi = 1.5$','interpreter','latex','fontsize',24)
xlabel('$\tau_d$','interpreter','latex','FontSize',20)
xlim([0 1])
ylim([0 1])
% legend({'Date-0 Inflation','Tax Base','No Eq''m'},'Location','Northeast','fontsize',18,'interpreter','latex')
legend({'Date-0 Inflation','Tax Base'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.1*pos(3) 1.15*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_c2','-dpng');

cd([path experiment]);