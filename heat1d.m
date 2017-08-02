%% Solving the Heat Equation in 1 Dimension via Finite Difference Method (Backwards Euler)
% Dependent on the function dudt in dudt.m
% By Michael Klamkin 2017
clear all; close all;
global A b

%% Constants 

n = 100; % number of points on the rod
c = 0.1; % thermal diffusion constant
dx = 1/n; % delta x
omega = c/dx.^2; % just to not type out k/dx.^2 5 times
u_0 = 80; u_end = 20; %initial conditions
%% Initialization 

A = zeros(n); 
for i = 1:n-1
    A(i,i) = -2*omega; %central diag
    A(i, i+1) = omega; %lower diag
    A(i+1, i) = omega; %upper diag
end
A = A(1:n-1, 1:n-1); %upper and lower diags should be n-1, center should be n. saves loop time

b = zeros(n-1, 1);
b(1) = omega * u_0; %initial conditions
b(end) = omega * u_end; %initial conditions

u_0 = b / omega; %initial conditions
%% ode45 and plotting via imagesc

[t, u] = ode45( @dudt, [0,3], u_0); % see help/doc for more info

figure(1);
%    subplot(1, 2, 1);
        pcolor(u'); shading interp; set(gca,'Ydir','reverse'); %pcolor allows for interp shading style
        colorbar; % display colorbar
% 
%     subplot(1, 2, 2);
%     imagesc(u'); colorbar;  %to compare with imagesc