clear all
close all
clc


%% Problem Initialization
np = 101;
dom = 1;
h = dom/(np-1);
x = 0:h:dom;
y = 0:h:dom;
rho= 1;
u=1;
v=1;
gamma=0.005;
Pe = (rho*u*h)/gamma % Peclet Number


%% Boundary Condition
T = zeros(np,np);
T(1,:) = 1;
T(2,:) = 1;
T(:,1) = 1;
T(:,2) = 1;
T_new = zeros(np,np);
T_new(1,:) = 1;
T_new(2,:) = 1;
T_new(:,1) = 1;
T_new(:,2) = 1;


%% Discretization using CD
error = 1;
iter = 0;
error_tol = 1e-7;

while error > error_tol
    for i = 3: (np-2)
        for j = 3: (np-2)
            a_E = gamma - (3*rho*u*h)/8 ;
    a_W = gamma + (7*rho*u*h)/8 ;
    a_N = gamma - (3*rho*u*h)/8 ;
    a_S = gamma + (7*rho*u*h)/8 ;
    a_SS = (rho*v*h)/8 ;
    a_WW = (rho*u*h)/8 ;
    a_P = (3*rho*u*h)/8 + (3*rho*v*h)/8 + gamma + gamma + gamma + gamma ;
    
    T_new (i,j) = ( T(i+1,j)*a_E + T(i-1,j)*a_W + T(i,j-1)*a_N + T(i,j+1)*a_S - T(i,j+2)*a_SS - T(i-2,j)*a_WW )/a_P ;
        end
    end
    iter = iter +1;
    error = 0;
    for i = 3: (np-2)
        for j = 3: (np-2)
    error = error + abs(T(i,j)-T_new(i,j)) ;
        end
    end
    T = T_new ;
end


%% Plotting
x1 = ((1:np)-1).*h ;
y1 = 1 - ((1:np)-1).*h ;
[X,Y] = meshgrid(x1,y1);
contourf(X,Y,T,30) 
colorbar 


%% centreline
plot(1-y,T(:,(np-1)/2),'ro')
