%% EQUILIBRIUM EXISTENCE/UNIQUENESS REGIONS
% Marios Angeletos, Chen Lian, & Christian Wolf
% this version: 06/17/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Dropbox/Research/self_finance/codes/ckw';
path = [local '/ecma_replication/analytical'];
experiment = '/determinacy';

cd([path experiment]);

%% PREPARATIONS

%----------------------------------------------------------------
% Set Grids
%----------------------------------------------------------------

n_phi    = 5000;
phi_grid = linspace(-1,1,n_phi)';

n_taud    = 2500;
taud_grid = linspace(0,1,n_taud);

%----------------------------------------------------------------
% Set Parameters
%----------------------------------------------------------------

beta      = 0.95;
sigma     = 1;
omega_OLG = 0.8;
omega_RA  = 1;
tau_y     = 1/3;
dy_bar    = 1.04;

%% MAIN COMPUTATIONS

%----------------------------------------------------------------
% Cutoffs
%----------------------------------------------------------------

phi_lb_OLG = ((1-beta*omega_OLG)/omega_OLG * (1-omega_OLG) * tau_y)...
    /(sigma * (1-beta) + (1-beta*omega_OLG)/omega_OLG * beta * dy_bar * (1-omega_OLG));

phi_lb_RA  = ((1-beta*omega_RA)/omega_RA * (1-omega_RA) * tau_y)...
    /(sigma * (1-beta) + (1-beta*omega_RA)/omega_RA * beta * dy_bar * (1-omega_RA));

%----------------------------------------------------------------
% Classify Regions
%----------------------------------------------------------------

det_indx_OLG = NaN(n_phi,n_taud); % 1: EU = (0,0), 2: EU = (1,1), 3: EU = (1,0)

for i_phi = 1:n_phi
    for i_taud = 1:n_taud
        phi = phi_grid(i_phi);
        taud = taud_grid(i_taud);
        
        taud_ub = 1 - beta - 1/(sigma*phi) * (1-beta*omega_OLG)/omega_OLG * (tau_y - beta * phi * dy_bar) * (1-omega_OLG);
        
        if phi > phi_lb_OLG
            if taud > taud_ub
                det_indx_OLG(i_phi,i_taud) = 2;
            else
                det_indx_OLG(i_phi,i_taud) = 1;
            end
        elseif phi > 0 && phi <= phi_lb_OLG
            det_indx_OLG(i_phi,i_taud) = 2;
        elseif phi < 0
            if taud > taud_ub
                det_indx_OLG(i_phi,i_taud) = 3;
            else
                det_indx_OLG(i_phi,i_taud) = 2;
            end
        end      
    end
end

det_indx_RA = NaN(n_phi,n_taud); % 1: EU = (0,0), 2: EU = (1,1), 3: EU = (1,0)

for i_phi = 1:n_phi
    for i_taud = 1:n_taud
        phi = phi_grid(i_phi);
        taud = taud_grid(i_taud);
        
        taud_ub = 1 - beta - 1/(sigma*phi) * (1-beta*omega_RA)/omega_RA * (tau_y - beta * phi * dy_bar) * (1-omega_RA);
        
        if phi > phi_lb_RA
            if taud > taud_ub
                det_indx_RA(i_phi,i_taud) = 2;
            else
                det_indx_RA(i_phi,i_taud) = 1;
            end
        elseif phi > 0 && phi <= phi_lb_RA
            det_indx_RA(i_phi,i_taud) = 2;
        elseif phi < 0
            if taud > taud_ub
                det_indx_RA(i_phi,i_taud) = 3;
            else
                det_indx_RA(i_phi,i_taud) = 2;
            end
        end      
    end
end

%% PLOT RESULTS

cd([pwd '/_results']);

plot_legend = {'None','Unique','Multiple'};

n = 200;
clear cmap

cmap(1,:) = 0.9 * [1 1 1];
cmap(2,:) = 0.7 * [1 1 1];
cmap(3,:) = 0.3 * [1 1 1];

[X,Y] = meshgrid([1:3],[1:50]);

cmap = interp2(X([1,25,50],:),Y([1,25,50],:),cmap,X,Y);

figure(1)
imagesc(phi_grid,taud_grid,det_indx_OLG')
colormap(cmap)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
set(gca,'YDir','normal')
xlim([-1 1])
ylim([0 1])
set(gca,'XTick',[-1:0.25:1]);
set(gca,'YTick',[0:0.25:1]);
set(gca,'TickLength',[0 0])
xlabel('$\phi$','interpreter','latex','FontSize',20);
ylabel('$\tau_d$','interpreter','latex','FontSize',20);

hidden_h = [];
hold on;
for K = 1:length(plot_legend)
    hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor', 'none');
end
hold off;
colormap(cmap);
uistack(hidden_h, 'bottom');
legend(hidden_h, plot_legend, 'Location', 'southoutside', 'NumColumns', 3, 'interpreter', 'latex','FontSize',16);

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_b1_2','-dpng');

figure(2)
imagesc(phi_grid,taud_grid,det_indx_RA')
colormap(cmap)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
set(gca,'YDir','normal')
xlim([-1 1])
ylim([0 1])
set(gca,'XTick',[-1:0.25:1]);
set(gca,'YTick',[0:0.25:1]);
set(gca,'TickLength',[0 0])
xlabel('$\phi$','interpreter','latex','FontSize',20);
ylabel('$\tau_d$','interpreter','latex','FontSize',20);

hidden_h = [];
hold on;
for K = 1:length(plot_legend)
    hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor', 'none');
end
hold off;
colormap(cmap);
uistack(hidden_h, 'bottom');
legend(hidden_h, plot_legend, 'Location', 'southoutside', 'NumColumns', 3, 'interpreter', 'latex','FontSize',16);

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_b1_1','-dpng');

cd ..