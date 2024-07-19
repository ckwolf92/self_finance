%% 2-PERIOD VISUAL ILLUSTRATION
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

kappa = 10^(-5);

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

global phi tau_d H rho_d tau_x_seq

% monetary rule

phi = 0;

% fiscal rule

tau_d  = [];
rho_d  = [];
H      = 70;

% policy shock

tau_x_seq = zeros(T,1);
tau_x_seq(1) = -1;

%% COMPUTE GE TRANSITION

%----------------------------------------------------------------
% Placeholders
%----------------------------------------------------------------

y_seq_all      = NaN(T,1);
y_seq_all_PE   = NaN(T,1);
d_seq_all      = NaN(T,1);
pi_seq_all     = NaN(T,1);
t_seq_all      = NaN(T,1);
t_seq_all_PE   = zeros(T,1);
nu_base_all    = NaN(1,1);
nu_rates_all   = NaN(1,1);
nu_prices_all  = NaN(1,1);
nu_all         = NaN(1,1);

%----------------------------------------------------------------
% PE Wedge
%----------------------------------------------------------------

hor = H+1;

t_seq_all_PE(1,1)   = 1;
t_seq_all_PE(hor,1) = t_seq_all_PE(hor,1)-beta^(-(hor-1));

y_seq_all_PE(:,1) = C_y * t_seq_all_PE(:,1);

clear hor

%----------------------------------------------------------------
% GE Transition
%----------------------------------------------------------------

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

y_seq_all(:,1)     = y_seq;
d_seq_all(:,1)     = d_seq;
pi_seq_all(:,1)    = pi_seq;
t_seq_all(:,1)     = t_seq;

expenditure = -(r_NPV' * tau_x_seq);

nu_base_all(1)   = (r_NPV' * Y_SS * tau_y * y_seq) ./ expenditure;
nu_rates_all(1)  = (-r_NPV(2:T)' * D_SS * (i_seq(1:T-1) - pi_seq(2:T))) ./ expenditure;
nu_prices_all(1) = (D_SS * pi_seq(1)) ./ expenditure;
nu_all(1)        = nu_base_all(1) + nu_rates_all(1) + nu_prices_all(1);

t_e_seq_all = t_seq_all - tau_y * y_seq_all;

%% PLOT RESULTS

%----------------------------------------------------------------
% Color Preparation
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.lpurple = 0.25 * [160/255 32/255 240/255] + 0.75 * [1 1 1];
settings.colors.lorange = 0.25 * [204/255 102/255 0/255] + 0.75 * [1 1 1];
settings.colors.lblue   = 0.25 * [116/255 158/255 178/255] + 0.75 * [1 1 1];
settings.colors.lgrey   = 0.75 * [230/255 230/255 230/255] + 0.25 * [1 1 1];

settings.colors.list = [0 * [0.5 0.5 0.5]; 0.5 * [0.5 0.5 0.5]; 1 * [0.5 0.5 0.5]; 1.5 * [0.5 0.5 0.5]; 1.8 * [0.5 0.5 0.5]];

settings.colors.finance = [0.25 * [196/255 174/255 120/255] + 0.75 * [1 1 1]; ...
                            0.5 * [0 0 0] + 0.5 * [1 1 1]; ...
                            0.25 * [189/255 142/255 131/255] + 0.75 * [1 1 1]];

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

IRF_plot = 100;

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
plot(0,0,'linewidth',5,'linestyle','-','color',settings.colors.blue)
hold on
plot(0,0,'linewidth',5,'linestyle','-','color',settings.colors.list(3,:))
hold on
jbfill(0:1:IRF_plot,(zeros(IRF_plot+1,1))',(y_seq_all_PE(1:IRF_plot+1,1))',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
jbfill(0:1:40,(y_seq_all_PE(1:40+1,1))',(y_seq_all(1:40+1,1))',...
    settings.colors.lgrey,settings.colors.lgrey,0,1);
plot(0:1:IRF_plot,y_seq_all_PE(1:IRF_plot+1,1),'linewidth',5,'linestyle','-','color',settings.colors.blue)
hold on
plot(0:1:IRF_plot,y_seq_all(1:IRF_plot+1,1),'linewidth',5,'linestyle','-','color',settings.colors.list(3,:))
hold on
set(gcf,'color','w')
xlabel('$t$','interpreter','latex','FontSize',20)
ylabel('\%','interpreter','latex','FontSize',20)
ylim([-0.2 0.5])
legend({'PE', 'GE'},'Orientation','horizontal','Location','Southwest','fontsize',18,'interpreter','latex')
text(15,0.45,'"short run"','interpreter','latex','FontSize',20)
text(65,0.45,'"long run"','interpreter','latex','FontSize',20)
annotation('textarrow',[0.225 0.175],[0.5 0.45],'String',' $1/\tau_y$','interpreter','latex','FontSize',20)
annotation('textarrow',[0.16 0.16],[0.28 0.38],'String',' $1$','interpreter','latex','FontSize',20)
annotation('textarrow',[0.6725 0.6725],[0.38 0.28],'String',' $-1$','interpreter','latex','FontSize',20)
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.2*pos(3) 1.2*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_2','-dpng');

cd([path experiment]);