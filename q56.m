clc;
clear;

%% Step 1: Set the parameters

r_star = 1.04^0.25-1;         % long-run interest rate
beta   = 1/(1+r_star); % discount factor
delta  = 1-(0.9)^0.25;          % depreciation rate
alpha  = 0.32;         % capital income share
sigma  = 2;            % inverse of IES
omega   = 1.455;
phi_k   = 0.028;
psi     = 0.01;
rho_A     = 0.5;
rho_r   = 0.7;
sig_eps = 1;

%% Step 2: Calculate steady states

% Calculate steady states
r_ss = r_star; % steady-state r
kappa = ((r_star + delta)/alpha)^(1/(alpha-1)); % *k/h
h_ss = ((1-alpha)*(kappa^alpha))^(1/(omega-1)); % steady-state h
k_ss = kappa*h_ss; % steady-state k
i_ss = delta*k_ss; % steady-state i
y_ss = (k_ss^alpha)*(h_ss^(1-alpha)); % steady-state y
tb_y = 0.01; % steady-state trade-balance-to-output ratio
d_bar = (tb_y/r_star)*y_ss; % d_bar
d_ss = d_bar; % steady-state external debt
c_ss = -(r_star*d_bar)/(1+r_star) + (kappa^alpha)*h_ss - delta*k_ss; % steady-state c
tb_ss = r_star*d_bar; % steady-state trade balance
ca_ss = 0; % steady-state current account
w_ss = (1-alpha)*(k_ss^alpha)*(h_ss^(-alpha)); % steady-state real wages


% Steady-state matrix
save('ss.mat','r_ss','d_ss','h_ss','k_ss','c_ss','i_ss','tb_ss','ca_ss','y_ss', 'w_ss');

% Parameter matrix
save('param.mat','sigma','delta','r_star','alpha','omega','beta','d_bar','psi','phi_k', 'rho_A', 'rho_r', 'sig_eps');


% Run dynare
addpath('C:\dynare\6.0\matlab')
dynare q56.mod;

% Impulse Responses
T = options_.irf; % number of periods
t = 1:1:T; % time (x-axis for irfs)

figure(1); sgtitle('Impulse Response to a TFP shock');
subplot(3,3,1); plot(t,100*oo_.irfs.y_eps_a,'b','LineWidth',1.5); title('Output');
subplot(3,3,2); plot(t,100*oo_.irfs.c_eps_a,'b','LineWidth',1.5); title('Consumption');
subplot(3,3,3); plot(t,100*oo_.irfs.h_eps_a,'b','LineWidth',1.5); title('Hours');
subplot(3,3,4); plot(t,100*oo_.irfs.i_eps_a,'b','LineWidth',1.5); title('Investment');
subplot(3,3,5); plot(t,100*oo_.irfs.w_eps_a,'b','LineWidth',1.5); title('Wage Rate');
subplot(3,3,6); plot(t,100*oo_.irfs.tb_eps_a,'b','LineWidth',1.5); title('Trade balance');
subplot(3,3,7); plot(t,100*oo_.irfs.tb_y_eps_a,'b','LineWidth',1.5); title('Trade balance/output');
saveas(gcf,'IR_tfp.png');


figure(2); sgtitle('Impulse Response to an interest rate shock')
subplot(3,3,1); plot(t,100*oo_.irfs.y_eps_r,'b','LineWidth',1.5); title('Output');
subplot(3,3,2); plot(t,100*oo_.irfs.c_eps_r,'b','LineWidth',1.5); title('Consumption');
subplot(3,3,3); plot(t,100*oo_.irfs.h_eps_r,'b','LineWidth',1.5); title('Hours');
subplot(3,3,4); plot(t,100*oo_.irfs.i_eps_r,'b','LineWidth',1.5); title('Investment');
subplot(3,3,5); plot(t,100*oo_.irfs.w_eps_r,'b','LineWidth',1.5); title('Wage Rate');
subplot(3,3,6); plot(t,100*oo_.irfs.tb_eps_r,'b','LineWidth',1.5); title('Trade balance');
subplot(3,3,7); plot(t,100*oo_.irfs.tb_y_eps_r,'b','LineWidth',1.5); title('Trade balance/output');
saveas(gcf, 'IR_ir.png')