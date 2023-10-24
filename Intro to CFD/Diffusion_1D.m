clear all 
close all
clc

np = 5; % point number
comp_dom = 1; % Domain length
h= comp_dom/(np-1); % Spacing

%% Initial Condition
y(np) = 1;
y(1) = 0;

y_new(np) = 1;
y_new(1) = 0;

%% Main execution

error_mag = 1;
error_tol = 1e-7;
iter = 0;

while error_mag > error_tol
    for i = 2:(np-1)
            y_new(i) = 0.5 .* (y(i-1)+y(i+1));
        iter = iter +1;
        end
        
    
    error_mag = 0;
    for i = 2:(np-1)
            error_mag = error_mag + abs(y(i)-y_new(i));
    end
    y=y_new ;
end        
    
    


%% PLOTTING

x = ((1:np)-1).*h ;
plot(x,y)
