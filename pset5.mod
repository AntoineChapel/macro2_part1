/*
 * 
 * 
 */

// Endogenous variables
var cT
    yN
    h
    lambda
    Dstar
    i
    epsilon 
    s
    w
    piN
    ptilde
    pvmr
    pvmc
    p
    y 
    pi

    big_eps
    e
    uid
    big_p

    eta
    ;

// Exogenous variables (shocks)
varexo nu;    // shocks to interest rate
       


// Parameters

parameters sigma
           istar
           alpha
           beta
           theta
           mu  
           tau
           psi
           phi 
           chi 
           hbar 
           a 
           xi 
           yT 
           rho_nu
           sigma_nu
           ;

load param.mat;

set_param_value('sigma', sigma);
set_param_value('istar', istar);
set_param_value('alpha', alpha);
set_param_value('beta', beta);
set_param_value('theta', theta);
set_param_value('mu', mu);
set_param_value('tau', tau);
set_param_value('psi', psi);
set_param_value('phi', phi);
set_param_value('chi', chi);
set_param_value('hbar', hbar);
set_param_value('a', a);
set_param_value('xi', xi);
set_param_value('yT', yT);
set_param_value('rho_nu', rho_nu);
set_param_value('sigma_nu', sigma_nu);



/* Model: TBD*/

model;
//Foreign MC
cT + (1+istar)*log(Dstar(-1)) = yT + log(Dstar) - (psi/2)*log(Dstar)^2;
//Traded Lagrangian
lambda = a*cT^(-1/xi);
//NT Lagrangian
lambda*p = (1-a)*yN^(-1/xi);
//Labor Lagrangian
phi*(hbar - h)^(-chi) = lambda*w;
//Domestic IR
lambda*epsilon(+1) = beta*(i)*lambda(+1);
//Intertemporal
lambda*(1-psi*log(Dstar)) = beta*(1+istar)*lambda(+1);
//Production
yN = (h/s)^alpha;
//intermediate goods dispersion
s = theta*s(-1)*(piN)^(mu/alpha) + (1-theta)*(ptilde)^(-mu/alpha);
//price
1 = theta*(piN)^(mu-1) + (1-theta)*(ptilde)^(1-mu);
//mc=mr
pvmc=pvmr;
pvmc = ((1-tau)/alpha)*yN^(1/alpha)*w*ptilde^(-mu/alpha) + beta*theta*(lambda(+1)/lambda)*(ptilde/ptilde(+1)*(1/(piN(+1))))^(-mu/alpha)*pvmc(+1);
pvmr = ((mu-1)/mu) *yN*p*(ptilde)^(1-mu) + beta*theta*(lambda(+1)/lambda)*((ptilde/ptilde(+1))*(1/piN(+1)))^(1-mu)*pvmr(+1);
//er
p = p(-1)*(piN/epsilon);
//Total output
y = ((1-a)*(yN)^(-1/xi)*yN + a*cT^(-1/xi)*yT)*(a*cT^(1-1/xi) + (1-a)*yN^(1-1/xi))^(1/(1-1/xi)-1);
//Gross inflation
pi = piN*((((yN(-1)))/((yN)))^(-(1/xi)))*(((a*((cT(-1))^(1-(1/xi))) + (1-a)*((yN(-1)))^(1 - (1/xi)))^((1/(1-(1/xi)))-1))/((a*((cT)^(1-(1/xi))) + (1-a)*((yN))^(1 -(1/xi)))^((1/(1-(1/xi)))-1)));



//ETA process
log(eta) = rho_nu*log(eta(-1)) + sigma_nu*nu;
//Taylor rule
log((i)/(1+istar)) = 1.5*log(pi) + 0.125*log(y) + log(eta);



//Real to nominal variables
big_eps = big_eps(-1)*epsilon;
big_p = big_p(-1)*pi;
e = big_eps/big_p;

//uid
uid = ((i)/(1+istar))/(epsilon(+1));

end;

load ss.mat;

// Steady State Values
initval;
cT      = cT_ss; 
yN      = yN_ss; 
h       = h_ss; 
lambda  = lambda_ss; 
Dstar   = exp(Dstar_ss); 
i       = 1+i_ss; 
epsilon = epsilon_ss; 
s       = s_ss; 
w       = w_ss; 
piN     = piN_ss; 
ptilde  = ptilde_ss; 
pvmr    = pvmr_ss; 
pvmc    = pvmc_ss; 
p       = p_ss; 
y       = y_ss; 
pi      = pi_ss;

big_eps = 1;
big_p = 1;
e = 1;

uid = 1;
eta = 1;
end;


// Approximate SS
steady(maxit=10000);
check;


//Manufacturing productivity shocks
shocks;
var nu; stderr 1;
end;


stoch_simul(order = 1, loglinear, irf=40, noprint, nodisplay)  pi piN epsilon big_eps e uid, i eta yN;
savefile = 'simulation.mat';
save(savefile, 'oo_');

