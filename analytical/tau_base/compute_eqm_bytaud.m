%% TRANSFER SHOCK IRFs FOR DIFFERENT TAU_D
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
cd([path experiment]);

%% PARAMETERS

%----------------------------------------------------------------
% Household Block
%----------------------------------------------------------------

global tau_y beta r_SS D_SS Y_SS D_i D_pi D_y C_i C_pi C_y

load inputs_olg

global T

T = size(C_y,1);

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

n_taud     = 51;
tau_d_list = linspace(1,0,n_taud);

H = 1;

% policy shock

tau_x_seq = zeros(T,1);
tau_x_seq(1) = -1;

%% COMPUTE GE TRANSITION

%----------------------------------------------------------------
% Placeholders
%----------------------------------------------------------------

y_seq_all      = NaN(T,n_taud);
d_seq_all      = NaN(T,n_taud);
pi_seq_all     = NaN(T,n_taud);
t_seq_all      = NaN(T,n_taud);
nu_base_all    = NaN(n_taud,1);
nu_rates_all   = NaN(n_taud,1);
nu_prices_all  = NaN(n_taud,1);
nu_all         = NaN(n_taud,1);

%----------------------------------------------------------------
% Loop
%----------------------------------------------------------------

for i_taud = 1:n_taud

if mod(i_taud,5) == 0
    disp(['I am at run ' num2str(i_taud)])
end

tau_d = tau_d_list(i_taud);

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

y_seq_all(:,i_taud)     = y_seq;
d_seq_all(:,i_taud)     = d_seq;
pi_seq_all(:,i_taud)    = pi_seq;
t_seq_all(:,i_taud)     = t_seq;

expenditure = -(r_NPV' * tau_x_seq);

nu_base_all(i_taud)   = (r_NPV' * Y_SS * tau_y * y_seq) ./ expenditure;
nu_rates_all(i_taud)  = (-r_NPV(2:T)' * D_SS * (i_seq(1:T-1) - pi_seq(2:T))) ./ expenditure;
nu_prices_all(i_taud) = (D_SS * pi_seq(1)) ./ expenditure;
nu_all(i_taud)        = nu_base_all(i_taud) + nu_rates_all(i_taud) + nu_prices_all(i_taud);

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
% Settings
%----------------------------------------------------------------

IRF_plot = 40;

taud_select = [0 25 35 45 50];

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
for iindx_taud = 1:length(taud_select)
    i_taud = taud_select(iindx_taud)+1;
    plot(0:1:IRF_plot,y_seq_all(1:IRF_plot+1,i_taud),'linewidth',5,'linestyle','-','color',settings.colors.list(iindx_taud,:))
    hold on
end
set(gcf,'color','w')
title('Output $y_t$','interpreter','latex','fontsize',24)
xlabel('$t$','interpreter','latex','FontSize',20)
ylabel('\%','interpreter','latex','FontSize',20)
ylim([0 0.5])
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
for iindx_taud = 1:length(taud_select)
    i_taud = taud_select(iindx_taud)+1;
    plot(0:1:IRF_plot,d_seq_all(1:IRF_plot+1,i_taud),'linewidth',5,'linestyle','-','color',settings.colors.list(iindx_taud,:))
    hold on
end
set(gcf,'color','w')
title('Gov''t Debt $d_t$','interpreter','latex','fontsize',24)
xlabel('$t$','interpreter','latex','FontSize',20)
legend({'$\tau_d = 1$','$\tau_d = 0.5$', ...
    '$\tau_d = 0.3$','$\tau_d = 0.1$','$\tau_d = 0$'},'Location','Northeast','fontsize',18,'interpreter','latex')
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
jbfill(tau_d_list,(squeeze(0*nu_prices_all))',(squeeze(nu_prices_all))',settings.colors.finance(1,:),settings.colors.finance(1,:),0,1);
jbfill(tau_d_list,(squeeze(nu_prices_all))',(squeeze(nu_prices_all+nu_base_all))',settings.colors.finance(2,:),settings.colors.finance(2,:),0,1);
hold on
set(gcf,'color','w')
title('Self-Financing Share $\nu$','interpreter','latex','fontsize',24)
xlabel('$\tau_d$','interpreter','latex','FontSize',20)
xlim([0 1])
legend({'Date-0 Inflation','Tax Base'},'Location','Northeast','fontsize',18,'interpreter','latex')
ylim([0 1])
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.1*pos(3) 1.15*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_1_1','-dpng');

cd([path experiment]);