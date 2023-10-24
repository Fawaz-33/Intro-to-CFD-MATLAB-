clear all
close all
clc

np = 101;
dom = 1;
h = dom/(np-1);
x = 0:h:dom;
rho= 1;
u=1;
gamma=1;
Pe = (rho*u*dom)/gamma; % Peclet Number

%% Initial Condition

phi(1) = 10;
phi(np) = 20;

phi_new(1) = 10;
phi_new(np) = 20;


%% Exact Solution
phi_exact = phi(1) + (phi(np)-phi(1)).*((exp(Pe.*x)/dom-1)./(exp(Pe)-1));



%% General Solution
error =1;
error_tol = 1e-07;
iter = 0;

while error > error_tol
    for i= 2:(np-1)
        a_E = (gamma/h) - (rho*u)/2;
        a_W = (gamma/h) + (rho*u)/2;
        a_P = (gamma/h) - (rho*u)/2 + (gamma/h) + (rho*u)/2;
        phi_new(i) = ( a_E*phi(i+1) + a_W*phi(i-1) ) / a_P ;
    end
    iter = iter + 1;
    error = 0;
    for i = 2:(np-1)
        error = error + abs(phi(i)-phi_new(i));
    end
    phi = phi_new;
end

%% Plotting
plot(x,phi_exact, x,phi,'ro')
  
