/*
 * 
 * 
 */

// Endogenous variables
var c           // log consumption
    h           // log labor
    lambda      // log Lagrange multiplier
    k           // log capital stock
    i           // log investment
    d           // external debt
    r           // interest rate
    w           // real wage
    tb          // trade balance
    ca          // current account
    y           // log output
    tb_y        // trade-balance-to-output ratio
    ca_y        // current-account-to-output ratio    
    a           // total factor productivity
    ;

// Exogenous variables (shocks)
varexo eps_a    // shocks to a
       eps_r    // shocks to r
       ;


// Parameters

parameters sigma
           delta
           r_star
           alpha
           omega
           beta
           d_bar
           psi
           phi_k
           rho_A
           rho_r
           sig_eps
           ;

load param.mat;

set_param_value('sigma',sigma);
set_param_value('delta',delta);
set_param_value('r_star',r_star);
set_param_value('alpha',alpha);
set_param_value('omega',omega);
set_param_value('beta',beta);
set_param_value('d_bar',d_bar);
set_param_value('psi',psi);
set_param_value('phi_k',phi_k);
set_param_value('rho_A', rho_A);
set_param_value('rho_r', rho_r);
set_param_value('sig_eps',sig_eps);

/* Model*/

model;
// Lagrange multiplier
(exp(c) - ((exp(h)^omega)/omega))^(-sigma) = exp(lambda);
// Euler equation for external debt
exp(lambda)*(1-psi*(d/(1+r)-d_bar/(1+r_star))) = beta*(1+r)*exp(lambda(+1));
// Intratemporal condition
exp(h)^(omega-1) = (1-alpha)*exp(a)*(exp(k(-1))^alpha)*(exp(h)^(-alpha));
// Budget constraint
exp(c) + exp(i) + (phi_k/2)*((exp(k)-exp(k(-1)))^2) + (psi/2)*((d/(1+r)-d_bar/(1+r_star))^2) + d(-1) = d/(1+r) + exp(y);
// Euler equation for capital
exp(lambda)*(phi_k*(exp(k)-exp(k(-1))) + 1) = beta*exp(lambda(+1))*(alpha*exp(a(+1))*(exp(k)^(alpha-1))*(exp(h(+1))^(1-alpha)) + phi_k*(exp(k(+1))-exp(k)) + 1- delta);
// Law of motion for capital
exp(k) = (1-delta)*exp(k(-1)) + exp(i);
// Output
exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
// Trade Balance
tb = d(-1) - d/(1+r);
// Current Account
ca = tb - d;
// Trade-balance-to-output ratio
tb_y = tb/exp(y);
// Current-account-to-output ratio
ca_y = ca/exp(y);
// Wages
exp(w) = (1-alpha)*exp(a)*(exp(k(-1))^alpha)*(exp(h)^(-alpha));
// TFP
a = rho_A*a(-1) + sig_eps*eps_a;
//interest rate
log(1+r) - log(1+r_star) = rho_r*(log(1+r(-1))- log(1+r_star)) + sig_eps*eps_r;
end;

load ss.mat;

// Steady State Values
initval;
c      = log(c_ss);
h      = log(h_ss);
lambda = -sigma*log(c_ss - ((h_ss^omega)/omega));
k      = log(k_ss);
i      = log(i_ss);
d      = d_ss;
r      = r_ss;
w      = log(w_ss);
tb     = tb_ss;
ca     = ca_ss;
y      = log(y_ss);
tb_y   = tb/y_ss;
ca_y   = ca/y_ss;
a      = 0; 
end;


// Approximate SS
steady(maxit=10000);
check;


//Manufacturing productivity shocks
shocks;
var eps_a; stderr 0.01;
var eps_r; stderr 0.01;
end;


stoch_simul(order = 1, irf=40, noprint, nodisplay) y h i tb_y c a tb w r;

savefile = 'simulation.mat';
save(savefile, 'oo_');

