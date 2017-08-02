%% Solving the Heat Equation in 1 Dimension via Finite Difference Method (Backwards Euler)
% dudt function for use in ode45, relies on global A and b
% By Michael Klamkin 2017
function out = dudt( t, u )

    global A b    
    out = (A*u + b);    
    
end

