clear all 
close all
clc

np = 31; % Defining no of points
comp_dom = 1;
h= comp_dom/(np-1);
dt = 0.0001; 
alpha = dt/(h*h);


%% Initial Condition

y(np,np) = 0;
y(1,:) = 1;

y_new(np,np) = 0;
y_new(1,:) = 1;

y_transient = 0;

error_mag = 1;
error_tol = 1e-7;
iter = 0;


%% Iteration part

while error_mag > error_tol
    for i = 2:(np-1)
        for j = 2:(np-1)
            y_new(i,j) = y(i,j) + alpha .* (y(i,j-1)+y(i,j+1)+y(i+1,j)+y(i-1,j)-(4*y(i,j)));
            y_transient(iter+1,1:np,1:np) = y_new;
        end
        
    end
    iter = iter+1;
    
    error_mag = 0;
    for i = 2:(np-1)
        for j = 2:(np-1)
            error_mag = error_mag + abs(y(i,j)-y_new(i,j));
        end
        
        
    end
    
    y=y_new ;
end



%% Particular Timestep Calc
timestep = 100;
x = ((1:np)-1).*h ;
y1 = 1 - ((1:np)-1).*h ;
[X,Y] = meshgrid (x,y1);
y_timestep = y_transient(timestep,:,:);
y_timestep = reshape (y_timestep,[np,np]);
contourf(X,Y,y_timestep,12)
colorbar
title(['Time = ' num2str(timestep*dt) 's'])


%% Animating Timestep
t_array = 1:timestep:iter;

for i = 1:length(t_array)
    timestep = t_array(i);
    y_timestep = y_transient(timestep,:,:);
    y_timestep = reshape (y_timestep,[np,np]);
    contourf(X,Y,y_timestep,12)
    colorbar
    title(['Time = ' num2str(timestep*dt) 's'])
    pause (0.25);
end
