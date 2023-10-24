clear all
close all
clc

nx = 51;
ny = 301;
domx =1;
domy =6;
hx = domx/(nx-1);
hy = domy/(ny-1);
x = 0:hx:domx;
y = 0:hy:domy;
dt = 0.001;
delta = 4.5;
Re = 100;


%% Initial Condition
% Collocated grid
u_final(nx,ny) = 0;
v_final(nx,ny) = 0;
P_final(nx,ny) =1;
u_final(1,:) = 1;

% Staggered grid
u(nx+1,ny) = 0;
v(nx,ny+1) = 0;
P(nx+1,ny+1) = 1;
u(1,:) = 2;

u_new(nx+1,ny) = 0;
v_new(nx,ny+1) = 0;
P_new(nx+1,ny+1) = 1;
u_new(1,:) = 2;


%% Governing Equations
error = 1;
iter = 0;
error_tol = 1e-6;
while error > error_tol
    % X- Momentum
    for i = 2:nx
      for  j = 2:(ny-1)
          Conv_1 = ((0.5*(u(i,j)+u(i,j+1)))^2 - (0.5*(u(i,j)+u(i,j-1)))^2 )/hx;
          Conv_2 = (0.25*((u(i,j)+u(i-1,j))*(v(i-1,j)+v(i-1,j+1)))- 0.25*((u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))))/hx;
          Press_G = -(P(i,j+1)-P(i,j))/hx;
          Diff = (1/Re)* ((u(i+1,j)-2*u(i,j)+u(i-1,j))/(hx*hx) + (u(i,j+1)-2*u(i,j)+u(i,j-1))/(hx*hx));
          
          u_new(i,j) = u(i,j)+ dt*(Press_G + Diff - Conv_1 - Conv_2);
      end
    end
    
      % BC- X - Momentum
      u_new(1,:) = - u_new(2,:); % Upper
      u_new(nx+1,:) = -u_new(nx,:); % Lower
      u_new(2:nx,1) = 1; 
      u_new(2:nx,ny) = u_new(2:nx,ny-1);
      
      % Y- Momentum
    for i = 2:(nx-1)
      for  j = 2:ny
          Conv_1y = ((0.5*(v(i-1,j)+v(i,j)))^2 - (0.5*(v(i,j)+v(i+1,j)))^2 )/hx;
          Conv_2y = (0.25*((u(i,j)+u(i+1,j))*(v(i,j+1)+v(i,j)))- 0.25*((u(i,j-1)+u(i+1,j-1))*(v(i,j)+v(i,j-1))))/hx;
          Press_Gy = -(P(i,j)-P(i+1,j))/hx;
          Diffy = (1/Re)* ((v(i,j-1)-2*v(i,j)+v(i,j+1))/(hx*hx) + (v(i+1,j)-2*v(i,j)+v(i-1,j))/(hx*hx));
          
          v_new(i,j) = v(i,j)+ dt*(Press_Gy + Diffy - Conv_1y - Conv_2y);
      end
    end
    
      % BC- Y - Momentum
      v_new(:,1) = -v_new(:,2); % Upper
      v_new(:,ny+1) = v_new(:,ny); % Lower
      v_new(1,2:ny) = 0; 
      v_new(nx,2:ny) = 0;
        
    
    
    % Continuity equation
    for i = 2:nx
        for j = 2:ny
            P_new(i,j) = P(i,j) - delta*dt*(u(i,j) - u(i,j-1) + v(i-1,j) - v(i,j))/hx;
        end
    end
    % BC - Continuity equation
    P_new(1,:)= P_new(2,:); % Upper
    P_new(nx+1,:)= P_new(nx,:); % Lower
    P_new(:,1)= P_new(:,2);
    P_new(:,ny+1)= - P_new(:,ny);
    
    % Error Calculation
    error = 0;
    for i = 2:(nx - 1)
        for j =2:(ny - 1)
           error = error + abs((u_new(i,j) - u_new(i,j-1) + v_new(i-1,j) - v_new(i,j))/hx);
        end
    end
    
    % Assigning values
    u = u_new;
    v = v_new;
    P = P_new;
    iter =iter+1;
       
end


%% Mapping staggered grid to collocated grid

    for i=1:nx
        for j= 1:ny
            u_final(i,j) = 0.5*(u(i,j)+u(i+1,j));
            v_final(i,j) = 0.5*(v(i,j)+v(i,j+1));
            P_final(i,j) =0.25*(P(i,j)+P(i+1,j)+P(i,j+1)+P(i+1,j+1));
        end
    end

%% Plotting
x1 = ((1:ny)-1).*hx;
y1 = 1-((1:nx)-1).*hx;
[X,Y] = meshgrid(x1,y1);
contourf(X,Y,u_final,40,'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title('Channel flow')
daspect([1 1 2])
