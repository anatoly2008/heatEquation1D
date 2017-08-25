%% Solving the Heat Equation in 1 Dimension via Finite Difference Method
% Dependent on the function dudt in dudt.m
% By Michael Klamkin 2017
clear all; close all;
global A b
ode45tic = tic; % timer for entire script
%% Constants 
n = 500; % number of points on the rod
c = 0.1; % thermal diffusion constant
dx = 1/n; % delta x
omega = c/dx.^2; % just to not type out k/dx.^2 5 times
tstop = 3; %length of simulation
%% Initialization 

A = zeros(n-1); 

for i = 1:n-2
    A(i,i) = -2*omega; %central diag
    A(i, i+1) = omega; %lower diag
    A(i+1, i) = omega; %upper diag
end

A(n-1,n-1) = -2 * omega;

for i = 1:n-1 %init cond
    b(i,1) = omega*(10*sin(pi*i*0.01));
end

u_0 = b / omega; %initial conditions
%% ode45 and plotting via imagesc

[t, u] = ode45( @dudt, [0,tstop], u_0); % see help/doc for more info

figure(1);
     subplot(2, 2, 1);
        fig1 = pcolor(fliplr(u')); shading interp; %colormap jet %set(gca,'Ydir','reverse'); %pcolor allows for interp shading style
        cbar = colorbar;
        cbarpos = get(get(cbar,'YLabel'),'Position');
        set(get(cbar,'YLabel'),'Position',[cbarpos(1) + 3, cbarpos(2) + 40, cbarpos(3)],'Rotation',-90, 'FontSize', 15);
        title('1D Heat Equation - ode45()'); % figure title
        set(gca, 'FontName', 'Times New Roman'); set(gca, 'TitleFontSizeMultiplier', 1.25); set(gca, 'FontSize', 15); % figure styling
        ylabel('Position on Rod'); % y axis label
        xlabel('Time'); % x axis label 
display(toc(ode45tic), 'ode45()'); % final timer
%% Forward Euler For-Loop
% Heat Equation Solution using Forward Euler FEM
% By Michael Klamkin 2017
fwdtic = tic; % timer for entire script

%% Constants
dx = dx * 2; % pos step adjusted
c = c * 1000; % diff. const. adjusted
dt = 0.001; % time step
gm = dt/((c)*(dx^2)); % a const for use in fwd euler

%% Vectors
x = 0:dx:10; % vect of steps of pos on rod
t = 0:dt:7; % vect of steps of time

%% Initial Conditions 
row = length(t); % num time samples
col = length(x); % num pos samples
U = zeros(row,col); % build U
for i = 1:length(x) %init cond
    U(1,i) = 10*sin(pi*x(i));
end

%% for Loop + Plot
for i = 1:row-1
    for j = 2:col-1
        U(i+1,j) = gm*(U(i,j+1) + U(i,j-1)) + (1-2*gm)*(U(i,j));
    end
end

U(:,1) = []; U(:, end) = []; %shave 0-set rows/column

figure(1);
     subplot(2, 2, 2);
     fig2 = pcolor(fliplr(u')); shading interp; %colormap jet %set(gca,'Ydir','reverse'); %pcolor allows for interp shading style
        cbar2 = colorbar;
        cbarpos2 = get(get(cbar2,'YLabel'),'Position');
        set(get(cbar2,'YLabel'),'Position',[cbarpos2(1) + 3, cbarpos2(2) + 40, cbarpos2(3)],'Rotation',-90, 'FontSize', 15);
        title('1D Heat Equation - Forward Euler'); % figure title
        set(gca, 'FontName', 'Times New Roman'); set(gca, 'TitleFontSizeMultiplier', 1.25); set(gca, 'FontSize', 15); % figure styling
        ylabel('Position on Rod'); % y axis label
        xlabel('Time'); % x axis label 
display(toc(fwdfic), 'Forward Euler');

%% Backward Euler For-Loop
% Heat Equation Solution using Forward Euler FEM
% By Michael Klamkin 2017
backtic = tic;
%% Constants
dx = dx * 2; % pos step adjusted
c = c * 1000; % diff. const. adjusted
dt = 0.001; % time step
gm = dt/((c)*(dx^2)); % a const for use in fwd euler

%% Vectors
x = 0:dx:10; % vect of steps of pos on rod
t = 0:dt:7; % vect of steps of time

%% Initial Conditions 
row = length(t); % num time samples
col = length(x); % num pos samples
U = zeros(row,col); % build U
for i = 1:length(x) %init cond
    U(1,i) = 10*sin(pi*x(i));
end

%% for Loop + Plot
for i = 1:row-1
    for j = 2:col-1
        U(i,j) = -gm*(U(i+1,j-1) + U(i+1,j+1)) + (1+2*gm)*(U(i+1,j));
    end
end

U(:,1) = []; U(:, end) = []; %shave 0-set rows/column

figure(1);
     subplot(2, 2, 3);
     fig2 = pcolor(fliplr(u')); shading interp; %colormap jet %set(gca,'Ydir','reverse'); %pcolor allows for interp shading style
        cbar2 = colorbar;
        cbarpos2 = get(get(cbar2,'YLabel'),'Position');
        set(get(cbar2,'YLabel'),'Position',[cbarpos2(1) + 3, cbarpos2(2) + 40, cbarpos2(3)],'Rotation',-90, 'FontSize', 15);
        title('1D Heat Equation - Backward Euler'); % figure title
        set(gca, 'FontName', 'Times New Roman'); set(gca, 'TitleFontSizeMultiplier', 1.25); set(gca, 'FontSize', 15); % figure styling
        ylabel('Position on Rod'); % y axis label
        xlabel('Time'); % x axis label 
        
display(toc(backtic), 'Backward Euler');
