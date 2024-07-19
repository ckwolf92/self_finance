%% TRANSFER SHOCK IRFs FOR EMPIRICALLY RELEVANT FISCAL RULES IN HANK
% Marios Angeletos, Chen Lian, & Christian Wolf
% this version: 06/17/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Dropbox/Research/self_finance/codes/ckw';
path = [local '/ecma_replication/hank'];
experiment = '/tau_base';

addpath([path '/_inputs'])
addpath([path experiment '/_aux_ge'])
addpath([local '/ecma_replication/analytical/_auxiliary_functions'])
cd([path experiment]);

%% PARAMETERS

%----------------------------------------------------------------
% Household Block
%----------------------------------------------------------------

global gamma beta tau_y r_SS D_SS Y_SS Trans_SS D_i D_pi D_y D_tau C_i C_pi C_y C_tau

load inputs_hank

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

global phi tau_d H tau_x_seq

% monetary rule

phi = 0;

% fiscal rule

fiscal_targets = (1 - [0.3 0.1 0.015]).^(1/4);

n_taud     = 51;
tau_d_list = linspace(1,0,n_taud);

tau_d_list = [tau_d_list,1 - fiscal_targets];
n_taud     = length(tau_d_list);

H = 1; % note: alternative assumptions on H not coded up

% policy shock

tau_x_seq    = zeros(T,1);
tau_x_seq(1) = 1/Trans_SS;

%% COMPUTE GE TRANSITION

%----------------------------------------------------------------
% Placeholders
%----------------------------------------------------------------

y_seq_all      = NaN(T,n_taud);
d_seq_all      = NaN(T,n_taud);
pi_seq_all     = NaN(T,n_taud);
tau_seq_all    = NaN(T,n_taud);
nu_base_all    = NaN(n_taud,1);
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
tau_seq_all(:,i_taud)   = tau_seq;

expenditure = (r_NPV' * Trans_SS * tau_x_seq) + (-r_NPV(2:T)' * D_SS * (1+r_SS) * (i_seq(1:T-1) - pi_seq(2:T)));

nu_base_all(i_taud)   = (r_NPV' * Y_SS * tau_y * y_seq) ./ expenditure;
nu_prices_all(i_taud) = (D_SS * (1+r_SS) * pi_seq(1)) ./ expenditure;
nu_all(i_taud)        = nu_base_all(i_taud) + nu_prices_all(i_taud);

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

taud_select = [51 52 53];

%----------------------------------------------------------------
% Plot
%----------------------------------------------------------------

cd([path experiment '/_results']);

plotwidth = 0.25;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * (gapsize + plotwidth)];

indx_taud = length(taud_select);

figure(1)

subplot(1,3,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
for iindx_H = 1:indx_taud
    i_taud = taud_select(iindx_H)+1;
    plot(0:1:IRF_plot,y_seq_all(1:IRF_plot+1,i_taud),'linewidth',5,'linestyle','-','color',settings.colors.list(iindx_H,:))
    hold on
end
set(gcf,'color','w')
title('Output $y_t$','interpreter','latex','fontsize',24)
xlabel('$t$','interpreter','latex','FontSize',20)
ylabel('\%','interpreter','latex','FontSize',22)
ylim([-0.01 0.5])
yticks([0 0.1 0.2 0.3 0.4 0.5])
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
for iindx_H = 1:indx_taud
    i_taud = taud_select(iindx_H)+1;
    plot(0:1:IRF_plot,pi_seq_all(1:IRF_plot+1,i_taud),'linewidth',5,'linestyle','-','color',settings.colors.list(iindx_H,:))
    hold on
end
set(gcf,'color','w')
title('Inflation $\pi_t$','interpreter','latex','fontsize',24)
xlabel('$t$','interpreter','latex','FontSize',20)
legend({'$\tau_d = 0.085$','$\tau_d = 0.026$','$\tau_d = 0.004$'},'Location','Northeast','fontsize',18,'interpreter','latex')
ylim([-0.02*0.1 0.1])
yticks([0 0.025 0.05 0.075 0.1])
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
jbfill(tau_d_list(1:end-3),(squeeze(0*nu_prices_all(1:end-3)))',(squeeze(nu_prices_all(1:end-3)))',settings.colors.finance(1,:),settings.colors.finance(1,:),0,1);
jbfill(tau_d_list(1:end-3),(squeeze(nu_prices_all(1:end-3)))',(squeeze(nu_prices_all(1:end-3)+nu_base_all(1:end-3)))',settings.colors.finance(2,:),settings.colors.finance(2,:),0,1);
hold on
i_taud = taud_select(indx_taud)+1;
plot(tau_d_list(end-2),nu_all(end-2),'o','MarkerSize',10,...
    'MarkerEdgeColor',settings.colors.list(1,:),...
    'MarkerFaceColor',settings.colors.list(1,:))
hold on
plot(tau_d_list(end-1),nu_all(end-1),'o','MarkerSize',10,...
    'MarkerEdgeColor',settings.colors.list(2,:),...
    'MarkerFaceColor',settings.colors.list(2,:))
hold on
plot(tau_d_list(end),nu_all(end),'o','MarkerSize',10,...
    'MarkerEdgeColor',settings.colors.list(3,:),...
    'MarkerFaceColor',settings.colors.list(3,:))
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
print('figure_e2','-dpng');

cd([path experiment]);