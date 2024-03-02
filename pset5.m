clc;
clear;

%% Step 1: Set the parameters

sigma  = 2;          % Inverse of intertemp elast of consumption
istar = 0.0316;      % Steady state international ir.
alpha   = 0.75;      % Labor share in NT sector
beta  = 1/(1.0316);  % Quarterly discount factor
theta = 0.7;         % P(no price change in traded)
mu   =  6;           % Elasticity of subs across int goods
tau = 1/mu;
phi  = 1.11;         % Preference parameter
chi   = 1;           % Preference parameter
hbar   = 3;          % Labor endowment
a     = 0.26;        % Share of tradables
xi     = 0.5;        % Elast of subs across final goods
yT = 1;              % Steady state tradable output
psi = 1;
rho_nu = 0;          % Parameter persistence of eta_t
sigma_nu = 0.001;        % Variance of eta_t

%% Step 2: Calculate steady states

% Calculate steady states
cT_ss = yT; %%
D_ss = 0; %%
Dstar_ss = 0; %%
piN_ss = 1; %%
h_ss = 1; %%
lambda_ss = a; %%
p_ss = (1-a)/lambda_ss; %%
epsilon_ss = piN_ss;%%
i_ss = istar; %%
ptilde_ss = ((1-theta*(piN_ss^(mu-1)))/(1-theta))^(1-mu); %%
s_ss = (1-theta)*(ptilde_ss^(-mu/alpha))/(1-theta*piN_ss^(mu/alpha));
yN_ss = (h_ss/s_ss)^alpha; %%
w_ss = phi*(hbar-h_ss)^(-chi)/lambda_ss;
pvmc_ss = ((1-tau)/alpha)*(yN_ss^(1/alpha))*w_ss*ptilde_ss^(-mu/alpha)/(1 -beta*theta*(1/piN_ss)^(1-mu));
pvmr_ss = pvmc_ss; %%((mu - 1)*p_ss)/(mu*(1-beta*theta)); %%
y_ss = (a*(yN_ss)^(1 - 1/xi) + (1-a)*(yT)^(1 - 1/xi))^(1/(1-1/xi)); %%
pi_ss = piN_ss; %%

% Steady-state matrix
save('ss.mat', 'cT_ss', 'D_ss', 'Dstar_ss', 'piN_ss', 'h_ss', 'lambda_ss', 'p_ss', ...
     'epsilon_ss', 'i_ss', 'ptilde_ss', 's_ss', 'yN_ss', 'w_ss', 'pvmc_ss', 'pvmr_ss', ...
     'y_ss', 'pi_ss');

% Parameter matrix
save('param.mat', ...
    'sigma', 'istar', 'alpha', 'beta', 'psi', 'phi', 'chi', 'theta', 'mu', 'a', 'xi', 'tau', ...
    'yT', 'rho_nu', 'hbar');

% Run dynare
addpath('C:\dynare\6.0\matlab')
dynare pset5.mod;


%% Question 5
% Impulse Responses
T = options_.irf; % number of periods
t = 0:1:T-1; % time 
figure(1);
subplot(3,2,1); plot(t,sign(oo_.irfs.i_nu).*((1+abs(oo_.irfs.i_nu)).^4 - 1), 'b','LineWidth',1.5); title('Nominal Interest Rate');
subplot(3,2,2); plot(t,oo_.irfs.big_eps_nu,'b','LineWidth',1.5); title('Nominal Exchange Rate'); 
subplot(3,2,3); plot(t,oo_.irfs.e_nu,'b','LineWidth',1.5); title('Real Exchange Rate');
subplot(3,2,4); plot(t,oo_.irfs.pi_nu,'b','LineWidth',1.5); title('CPI Inflation');
subplot(3,2,5); plot(t,oo_.irfs.piN_nu,'b','LineWidth',1.5); title('Nontraded Goods Inflation' );
subplot(3,2,6); plot(t,oo_.irfs.uid_nu,'b','LineWidth',1.5); title('UID' );
saveas(gcf,'graph1.png');



%% Question 6
clear nom i_m;
for val = 1:3
    psi = val/2;
    save('ss.mat','cT_ss','D_ss','Dstar_ss','piN_ss','h_ss','lambda_ss','p_ss','epsilon_ss','i_ss','ptilde_ss','s_ss','yN_ss','w_ss','pvmc_ss','pvmr_ss', 'y_ss', 'pi_ss');
    save('param.mat','sigma','istar','alpha','beta','psi','phi','chi','theta','mu','a','xi','tau','yT','hbar');
    dynare pset5.mod;
    nom(val,:) = oo_.irfs.big_eps_nu;
    i_m(val,:) = sign(oo_.irfs.i_nu).*((1+abs(oo_.irfs.i_nu)).^4 - 1);
end
%

T = options_.irf; % number of periods
t = 0:1:T-1; % time (x-axis for irfs)
figure(2);
subplot(3,2,1); plot(t, i_m(1,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,2,2); plot(t, nom(1,:),'b','LineWidth',1.5); title('\psi = 0.5');
subplot(3,2,3); plot(t, i_m(2,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,2,4); plot(t, nom(2,:),'b','LineWidth',1.5); title('\psi = 1');
subplot(3,2,5); plot(t, i_m(3,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,2,6); plot(t, nom(3,:),'b','LineWidth',1.5); title('\psi = 1.5');
saveas(gcf,'graph2.png');



%% Question 7
clear nom eta_m i_m y_m pi_m;
sigma_nu = -0.004754;
sig_2 = -0.0061;
sig_5 = -0.034;
sig_6 = 0.02;
sig_9 = 0.00175;
sig_nu_m = [sigma_nu sigma_nu sig_2 sig_2 sig_2 sig_5 sig_6 sig_6 sig_6 sig_9];

for val=1:10
    rho = max(0,(val-1)/10);
    sigma_nu = sig_nu_m(val);
    save('ss.mat','cT_ss','D_ss','Dstar_ss','piN_ss','h_ss','lambda_ss','p_ss', 'epsilon_ss','i_ss','ptilde_ss','s_ss','yN_ss','w_ss','pvmc_ss','pvmr_ss', 'y_ss', 'pi_ss');
    save('param.mat','sigma','istar','alpha','beta','psi','phi','chi','theta','mu','a','xi','tau','yT','hbar');
    if val == 7 | val == 3 | val == 6 | val == 8 | val == 10
        dynare pset5.mod;
    end
    nom(val,:) = oo_.irfs.big_eps_nu;
    eta_m(val,:) = oo_.irfs.eta_nu;
    i_m(val,:) = sign(oo_.irfs.i_nu).*((1+abs(oo_.irfs.i_nu)).^4 - 1);
    y_m(val,:) = oo_.irfs.yN_nu;
    pi_m(val,:) = oo_.irfs.pi_nu;
    end
sigma_nu = -0.004754;
rho = 0;

% Impulse response graphs

T = options_.irf; % number of periods
t = 0:1:T-1; % time (x-axis for irfs)
figure(3);
subplot(3,3,1); plot(t,100*nom(3,:),'b','LineWidth',1.5); title('\epsilon, when \rho = 0.2' );
subplot(3,3,2); plot(t,100*eta_m(3,:),'b','LineWidth',1.5); title('\eta, when \rho = 0.2');
subplot(3,3,3); plot(t,i_m(3,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,3,4); plot(t,nom(7,:),'b','LineWidth',1.5); title('\epsilon, when \rho = 0.6' );
subplot(3,3,5); plot(t,eta_m(7,:),'b','LineWidth',1.5); title('\eta, when \rho = 0.6');
subplot(3,3,6); plot(t,i_m(7,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,3,7); plot(t,nom(10,:),'b','LineWidth',1.5); title('\epsilon, when \rho = 0.9' );
subplot(3,3,8); plot(t,eta_m(10,:),'b','LineWidth',1.5); title('\eta, when \rho = 0.9');
subplot(3,3,9); plot(t,i_m(10,:),'b','LineWidth',1.5); title('Interest Rate');
saveas(gcf,'graph3.png');



%% Question 8(a)

clear nom i_m;
sigma_nu = -0.0029;
for val = 1:3
    psi = val/2;
    % Steady-state matrix
    save('ss.mat','cT_ss','D_ss','Dstar_ss','piN_ss','h_ss','lambda_ss','p_ss', 'epsilon_ss','i_ss','ptilde_ss','s_ss','yN_ss','w_ss','pvmc_ss','pvmr_ss', 'y_ss', 'pi_ss');
    % Parameter matrix
    save('param.mat','sigma','istar','alpha','beta','psi','phi','chi','theta','mu','a','xi','tau','yT','hbar');
    dynare pset5.mod;
    nom(val,:) = oo_.irfs.big_eps_nu;
    i_m(val,:) = sign(oo_.irfs.i_nu).*((1+abs(oo_.irfs.i_nu)).^4 - 1);
end
%%
% Impulse Responses
T = options_.irf; % number of periods
t = 0:1:T-1; % time (x-axis for irfs)
figure(4);
subplot(3,2,1); plot(t,i_m(1,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,2,2); plot(t,nom(1,:),'b','LineWidth',1.5); title('\psi = 0.5');
subplot(3,2,3); plot(t,i_m(2,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,2,4); plot(t,nom(2,:),'b','LineWidth',1.5); title('\psi = 1');
subplot(3,2,5); plot(t,i_m(3,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,2,6); plot(t,nom(3,:),'b','LineWidth',1.5); title('\psi = 1.5');
saveas(gcf,'graph4.png');



%% Question 8(b)
clear nom eta_m i_m;
sigma_nu = -0.0029;
sig_2 = -0.00305;
sig_5 = -0.0379;
sig_6 = -0.0048;
sig_8 = 0.0253;
sig_9 = 0.00303;
sig_nu_m = [sigma_nu sigma_nu sig_2 sig_2 sig_2 sig_5 sig_6 sig_6 sig_8 sig_9];
for val = 1:10
    rho = max(0,(val-1)/10);
    sigma_nu = sig_nu_m(val);
    % Steady-state matrix
    save('ss.mat','cT_ss','D_ss','Dstar_ss','piN_ss','h_ss','lambda_ss','p_ss', 'epsilon_ss','i_ss','ptilde_ss','s_ss','yN_ss','w_ss','pvmc_ss','pvmr_ss', 'y_ss', 'pi_ss');
    % Parameter matrix
    save('param.mat','sigma','istar','alpha','beta','psi','phi','chi','theta','mu','a','xi','tau','yT','hbar');
    
    if val == 7 | val == 3 | val == 6 | val == 9 | val == 10
        dynare pset5.mod;
    end
    
    nom(val,:) = oo_.irfs.big_eps_nu;
    eta_m(val,:) = oo_.irfs.eta_nu;
    i_m(val,:) = sign(oo_.irfs.i_nu).*((1+abs(oo_.irfs.i_nu)).^4 - 1);
end
sigma_nu = -0.004754;
rho = 0;
%%
% Impulse Responses
T = options_.irf; % number of periods
t = 0:1:T-1; % time (x-axis for irfs)

figure(5);
subplot(3,3,1); plot(t,100*nom(3,:),'b','LineWidth',1.5); title('\epsilon, when \rho = 0.2' );
subplot(3,3,2); plot(t,eta_m(3,:),'b','LineWidth',1.5); title('\eta, when \rho = 0.2');
subplot(3,3,3); plot(t,i_m(3,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,3,4); plot(t,nom(9,:),'b','LineWidth',1.5); title('\epsilon, when \rho = 0.8' );
subplot(3,3,5); plot(t,eta_m(9,:),'b','LineWidth',1.5); title('\eta, when \rho = 0.8');
subplot(3,3,6); plot(t,i_m(9,:),'b','LineWidth',1.5); title('Interest Rate');
subplot(3,3,7); plot(t,nom(10,:),'b','LineWidth',1.5); title('\epsilon, when \rho = 0.9' );
subplot(3,3,8); plot(t,eta_m(10,:),'b','LineWidth',1.5); title('\eta, when \rho = 0.9');
subplot(3,3,9); plot(t,i_m(10,:),'b','LineWidth',1.5); title('Interest Rate');
saveas(gcf,'graph5.png');

