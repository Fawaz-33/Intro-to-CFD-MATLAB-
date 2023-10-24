clear all 
close all
clc

np = 51;
comp_dom = 1;
h= comp_dom/(np-1);

%% Initial Condition

y(np,np) = 0;
y(1,:) = 1;

y_new(np,np) = 0;
y_new(1,:) = 1;

%% Main execution

error_mag = 1;
error_tol = 1e-7;
iter = 0;

while error_mag > error_tol
    for i = 2:(np-1)
        for j = 2:(np-1)
            y_new(i,j) = 0.25 .* (y(i,j-1)+y(i,j+1)+y(i+1,j)+y(i-1,j));
        iter = iter +1;
        end
        
    end
    error_mag = 0;
    for i = 2:(np-1)
        for j = 2:(np-1)
            error_mag = error_mag + abs(y(i,j)-y_new(i,j));
        end
        
        
    end
    
    y=y_new ;
end


%% PLOTTING

x = ((1:np)-1).*h ;
y1 = 1 - ((1:np)-1).*h ;
[X,Y] = meshgrid (x,y1);
contourf(X,Y,y,12)
colorbar
