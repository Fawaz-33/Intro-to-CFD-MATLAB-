clear all
close all
clc

np = 61;
dom =1;
h = dom/(np-1);
x = 0:h:dom;
y = 0:h:dom;
dt = 0.001;
delta = 4.5;
Re = 100;


%% Initial Condition
% Collocated grid
u_final(np,np) = 0;
v_final(np,np) = 0;
P_final(np,np) =1;
u_final(1,:) = 1;

% Staggered grid
u(np+1,np) = 0;
v(np,np+1) = 0;
P(np+1,np+1) = 1;
u(1,:) = 2;

u_new(np+1,np) = 0;
v_new(np,np+1) = 0;
P_new(np+1,np+1) = 1;
u_new(1,:) = 2;


%% Governing Equations
error = 1;
iter = 0;
error_tol = 1e-6;
while error > error_tol
    % X- Momentum
    for i = 2:np
      for  j = 2:(np-1)
          Conv_1 = ((0.5*(u(i,j)+u(i,j+1)))^2 - (0.5*(u(i,j)+u(i,j-1)))^2 )/h;
          Conv_2 = (0.25*((u(i,j)+u(i-1,j))*(v(i-1,j)+v(i-1,j+1)))- 0.25*((u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))))/h;
          Press_G = -(P(i,j+1)-P(i,j))/h;
          Diff = (1/Re)* ((u(i+1,j)-2*u(i,j)+u(i-1,j))/(h*h) + (u(i,j+1)-2*u(i,j)+u(i,j-1))/(h*h));
          
          u_new(i,j) = u(i,j)+ dt*(Press_G + Diff - Conv_1 - Conv_2);
      end
    end
    
      % BC- X - Momentum
      u_new(1,:) = 2 - u_new(2,:); % Upper
      u_new(np+1,:) = -u_new(np,:); % Lower
      u_new(2:np,1) = 0; 
      u_new(2:np,np) = 0;
      
      % Y- Momentum
    for i = 2:(np-1)
      for  j = 2:np
          Conv_1y = ((0.5*(v(i-1,j)+v(i,j)))^2 - (0.5*(v(i,j)+v(i+1,j)))^2 )/h;
          Conv_2y = (0.25*((u(i,j)+u(i+1,j))*(v(i,j+1)+v(i,j)))- 0.25*((u(i,j-1)+u(i+1,j-1))*(v(i,j)+v(i,j-1))))/h;
          Press_Gy = -(P(i,j)-P(i+1,j))/h;
          Diffy = (1/Re)* ((v(i,j-1)-2*v(i,j)+v(i,j+1))/(h*h) + (v(i+1,j)-2*v(i,j)+v(i-1,j))/(h*h));
          
          v_new(i,j) = v(i,j)+ dt*(Press_Gy + Diffy - Conv_1y - Conv_2y);
      end
    end
    
      % BC- Y - Momentum
      v_new(:,1) = -v_new(:,2); % Upper
      v_new(:,np+1) = -v_new(:,np); % Lower
      v_new(1,2:np) = 0; 
      v_new(np,2:np) = 0;
        
    
    
    % Continuity equation
    for i = 2:np
        for j = 2:np
            P_new(i,j) = P(i,j) - delta*dt*(u(i,j) - u(i,j-1) + v(i-1,j) - v(i,j))/h;
        end
    end
    % BC - Continuity equation
    P_new(1,:)= P_new(2,:); % Upper
    P_new(np+1,:)= P_new(np,:); % Lower
    P_new(:,1)= P_new(:,2);
    P_new(:,np+1)= P_new(:,np);
    
    % Error Calculation
    error = 0;
    for i = 2:(np - 1)
        for j =2:(np - 1)
           error = error + abs((u_new(i,j) - u_new(i,j-1) + v_new(i-1,j) - v_new(i,j))/h);
        end
    end
    
    % Assigning values
    u = u_new;
    v = v_new;
    P = P_new;
    iter =iter+1;
       
end


%% Mapping staggered grid to collocated grid

    for i=1:np
        for j= 1:np
            u_final(i,j) = 0.5*(u(i,j)+u(i+1,j));
            v_final(i,j) = 0.5*(v(i,j)+v(i,j+1));
            P_final(i,j) =0.25*(P(i,j)+P(i+1,j)+P(i,j+1)+P(i+1,j+1));
        end
    end

%% Plotting
x1 = ((1:np)-1).*h;
y1 = 1-((1:np)-1).*h;
[X,Y] = meshgrid(x1,y1);
contourf(X,Y,u_final,40,'LineStyle', 'none')
colorbar
colormap('Hot')
xlabel('x')
ylabel('y')
title('Lid Driven Cavity')
